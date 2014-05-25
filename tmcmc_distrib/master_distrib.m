%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%PACKS DISTRIBUTES THE WORK GIVEN TO HIM OVER MULTIPLE
%%SLAVE PROCESSES
%%%%%%%%%%%%%%%%%
function [TMCMC]=master_distrib_cluster(data)
CoefVar(1)=10; %fake first initial CoefVar (arbitrary number)
logselection(1)=0;


%INIT OF MATRICES:
%thetacollector=samples,funceval=Estimations,flc=likelihood,fpc=posterior(p
%rior times likelihood
thetacollector=[];funceval=[];flc=[];fpc=[];
pflag=0;surro_index=[];output={};
%%%RUN KEEPS TRACK OF ALL THE POINTS AT ALL STAGES, AND ALSO THE REAL VS
%%%THE SURROGATES
runinfo.realpoints=[]; %theta values 
runinfo.realevals=[]; %real evaluations of observables (Dim times data.data.Nobs)
runinfo.surrpoints=[]; %surrogate estimates
runinfo.surratio=[];  %Ratio of surrogate to real points per stage/runinfo.gen 
runinfo.R2=[];  %second scaled covariance matrix
runinfo.iR=[]; %inverse of the scaled covarinace matrix
runinfo.currentuniques=[];
runinfo.current_surr=[];
runinfo.combined=[];    %together realpoints and funceval
runinfo.gen=1; %TMCMC STAGE (or generation) 
runinfo.q=[];
runinfo.p(1)=0; % p values p(1)=0;    %initial temparature for annealing at runinfo.gen=0
runinfo.meantheta=[];
runinfo.SS=eye(data.Nth,data.Nth);
runinfo.acceptance(1)=1;
runinfo.havekrig=zeros(data.maxstages,1);
%runinfo.krigpredictor=struct([]);
runinfo.thetacoll_orig=[];
runfin.fpc_orig=[];




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%PREPARE THE STARTING GENERATION IF WE ARE DOING TMCMC%%%%%%%%%%%%%%%%%
%Generation
%--------------------------------------------------------------------------
% GENERATION 1: Generate N_0 samples ?_{0,l}, l=1,...,N_0 from prior pdf ?(?) 
%           and compute likelihood f(?|D) for each sample 
%--------------------------------------------------------------------------

%POINTS THAT HAVE TO BE MONTE CARLO'd IN TMCMC FASHION 
%IF WE ARE DOING A HYPERCUBE SAMPLING IN ORDER TO CONTINUE WITH KRIGING,WE
%RESCALE ALL SAMPLES TO WITHIN THE DATA BOUNDS


if(data.Hypercube==1)
theta.j=lhsamp(data.Num(runinfo.gen),data.Nth); %Latin Hypercube samples (Multi-uniform)
theta.j=bsxfun(@times,diff(data.unifbounds),theta.j);
theta.j=bsxfun(@plus,data.unifbounds(1,:),theta.j);
theta.j1=theta.j; %Copy the points to a holder
else
 
%SAMPLE FROM THE PRIORS!!!!
    
 for i=1:1:data.Nth
    if (length(data.priorparam(i,:))==2)
theta.j(:,i)=random(data.priorf{i},data.priorparam{i,1},data.priorparam{i,2},data.Num(runinfo.gen),1);
    else
theta.j(:,i)=random(data.priorf{i},data.priorparam{i},data.Num(runinfo.gen),1);
    end
end
 
  
     theta.j1=theta.j;
 end





%%%%%%%%%%%%%%%%%NOW WE FOLLOW THE STANDAR FLOWCHART UNTIL TOL COV REACHED%

%while (CoefVar(runinfo.gen)>tolCOV && runinfo.gen<15 && p(runinfo.gen)~=1)
while (runinfo.gen<data.maxstages && pflag<1)

%%%%%%%Pack as equally as possible and ship to the cluster workers    
[output]=package_and_ship_work(theta.j,runinfo,data,thetacollector,flc,funceval,surro_index);

%%%%%%%EXPLICTLY HERE THERE IS A BARRIER.  
%%%%On the cluster it is implemented using LSF's "bsub -K -w { node lists}"


%%%%%%%%AFTER ALL THREADS /CORES OR NODES HAVE FINISHED, WE COLLECT AND
%RESHAPE EVERYTHING BACK
[thetacollector,funceval,flc,fpc,runinfo,TMCMC,surro_index]=collect_work(data,runinfo,logselection,CoefVar,output);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%WEIGHT AND ADAPT AND CALCULATE MODEL EVIDENCE%%%%%%%%%%%%%%%%%%%%%%%%%%
[logselection,CoefVar,pflag,data,runinfo]=calculate_statistics(flc,data,thetacollector,CoefVar,logselection,runinfo);



%%%%%%%%%%%%%%WE CAN NOW PROCEED TO THE NEXT GENERATION
fprintf('MESSAGE FROM MASTER: GENERATION %d IS NOW COMPLETE\n',runinfo.gen);
runinfo.gen=runinfo.gen+1;

%%%%%%%%%%%%%%%%%%%SAVE RESTART INFO!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
namesave=['stage_reg_',num2str(runinfo.gen)]; %saves restart information in case something crashes can be used to restart or refine the simulation
save(namesave);
end
save final
%close all the workers
% if (data.availcpus>1);
% matlabpool close;
% end






function [output]=package_and_ship_work(theta,runinfo,data,thetacollector,flc,funceval,surro_index)

%%%%%THIS IS WHERE ALL THE INPUTS ARE NEEDED!
%experimental input values


if((runinfo.gen==1 && ~data.refine) || (runinfo.gen==1 && data.totalrestart))

 points=theta;  
%NunChains= Number of Chains per worker
%FOR ALL SLAVES
 %If we are Latin Hypercubing (i.e probing the structure model) then divide
%equally the samples along the pool
 
  segmentlength=floor(max(size(points))/data.availcpus) , %segment length (theta's )per cpu for data eval

    for i=1:1:data.availcpus
            thetapercpu(:,:,i)=points((segmentlength*(i-1)+1):(segmentlength*(i)),:);
    end
fprintf('MESSAGE FROM MASTER: FIRST EVALUATION WORK TO %d AVAILABLE PROCESSORS\n',data.availcpus)
    for i=1:1:data.availcpus
      local_work={};
    %WRAP WORK MESSAGE, ALONG WITH PROPOSAL COVARIANCE FOR ACCEPTANCE
    for j=1:1:segmentlength
        
    local_work(j,:)={segmentlength,1, thetapercpu(j,:,i) , flc ,funceval,0};
    end
    
    
         work{i,:}=local_work;
    end

end
%%%%%%%%%%%%%%%%%%IF IT IS A NORMAL CHAIN EVENT NOT THE VERY FIRST ONE
if(runinfo.gen~=1 || (data.refine && ~data.totalrestart))
% data.availcpus = java.lang.Runtime.getRuntime.availableProcessors();

%%%%%%%%%%%%%%%%%%%%%%%%%FIND HOW MANY UNIQUE CHAINS WE HAVE    
          % Repeat Num(runinfo.gen) times. ell now gets the corresponging guys
          % based on their weights.

                     %HOW MANY SLAVES?         
                              Nslaves=data.availcpus();           
                 ell=(mnrnd(1,runinfo.q,data.Num(runinfo.gen)));
              %[ IX ] = reSample(runinfo.q);
              IX = tournamentSelection(data.Num(runinfo.gen),data.tournsize,runinfo.q');
                  for chainn=1:1:data.Num(runinfo.gen)
                 Depth_per_chain(chainn)=length(find(IX==chainn));
                      
                  end
              fprintf('UNIQUE!=%f\n\n',length(unique(IX)));
                   
              
             %   return;

                  
                %Shows how many steps of DRAM one has to do per worker
                %(weight for the load balancing)
              %  Depth_per_chain=sum(ell);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%each index of Depth per chain is how many times chain x got picked . Lets
%match the same chain picks and pack them into correct terms


%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              %%%%%%%%
%%%%CALCULATE EVEN DISTRIBUTION AMONG THE WORKERS  
[sorted_weight,weight_ind]=sort(Depth_per_chain,'descend');
distrib=pack_work(sorted_weight,Nslaves);     

load_balance=[];
for i=1:1:Nslaves
    load_balance(i)=sum(sorted_weight(distrib{i}));
end
fprintf('Done with packing, std between packets is:%f \n',std(load_balance)/mean(load_balance));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% The cell-type distrib now has the work indices. We have to map them back
%To each one of the original ones.
% weight_ind(distrib{i}) are the indices, poies allhsides phre o ka8e
% slave
work=cell(data.availcpus,1);
    for i=1:Nslaves
      local_work={};
    %counter=0; %init
    length(sorted_weight(distrib{i})),
          for chl=1:1:length(distrib{i})   %for as many chains this process got
      %        counter=counter+1; 
                        %PROCESS TAG    %LENGHT OF CHAIN chl   %points
                        %%evals %was it a surrogate?
     %   work(chl,:)={length(sorted_weight(distrib{i})), sorted_weight(distrib{i}(chl)),  thetacollector(qj(find(ell(distrib{i}(chl),:))),:) ,...
      %    flc(qj(find(ell(distrib{i}(chl),:))),:), funceval(qj(find(ell(distrib{i}(chl),:))),:) , surro_index(qj(find(ell(distrib{i}(chl),:))),:)}, %if none is found, this point bye bye
         
           local_work(chl,:)={length(sorted_weight(distrib{i})), sorted_weight(distrib{i}(chl)),  thetacollector(find(ell(distrib{i}(chl),:)),:) ,...
                 flc(find(ell(distrib{i}(chl),:)),:), funceval(find(ell(distrib{i}(chl),:)),:) , surro_index(find(ell(distrib{i}(chl),:)),:)}; %if none is found, this point bye bye
        % thetacollector(find(ell(distrib{i}(chl),:)),:);
          end  %end that's the work of this worker
   


%%%%%%%%%%%%%%%%%%%%%%%%%%COPYING PART%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
    fprintf('MESSAGE FROM MASTER: PREPARING WORK TO PROC %d of  %d AVAILABLE PROCESSORS\n',i,data.availcpus)
%   PROVIDE THE PACKET THE CORRECT NAME (eg : worker_in_1, worker_in_2)
  %    copy_files_toworker(i,data,work,runinfo.gen,q,p,SS,runinfo,flc,funceval,surro_index);

      work{i,:}=local_work;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

end  %of if runinfo.gen~=1
%%%%%%%%%%%%WORK MESSAGE READY, GIVE IT TO WORKER! %%%%%%%%%%%%%%%%%%%%%%%%
 %%USE JAVA RUN TIME TO SPAWN A PROCESS THAT WILL GIVE EACH ELEMENT TO ALL

 %tmcmc_worker(i),
%MULTICORE EXAMPLE SPAWN as many evaluators as you want:
if (~data.refine)
 try
 matlabpool open;
 catch
 end
end
output=cell(data.availcpus,1);

%this way if we have the matlabpool open we do not care to open it again
%use parfor if you want to go distributed!
for i=1:data.availcpus
 output{i,:}=tmcmc_worker(i,data,work{i,:},runinfo);
 %runtime = java.lang.Runtime.getRuntime();
 %pid = runtime.exec(['/Users/pa49/Desktop/mcmc_r0_1/tmcmc_worker.app/Contents/MacOS/applauncher',' ',num2str(i)]);
 %pid.waitFor(); % block until program returns.
end








function [distribution]=pack_work(weights,Nslaves)

w_length=length(weights);
work_distrib=zeros(Nslaves,1);
%distributoin : here we keep all the job indices
distribution=cell(Nslaves,1);
%for all weights (packages)
                for sorted_w_index=1:1:w_length
                %index_lazy=worker with less work allocated till now
                [value,index_lazy]=min(work_distrib);
                work_distrib(index_lazy)=work_distrib(index_lazy)+weights(sorted_w_index);
                 %add the index of packet to the work of the guy only if the
                 %guy is not a COMPLETE FAIL = (weight=0). 
                     if(weights(sorted_w_index))
                    distribution{index_lazy}=[distribution{index_lazy},sorted_w_index];
                     end
                end 


function [thetacollector,funceval,flc,fpc,runinfo,TMCMC,surro_index]=collect_work(data,runinfo,logselection,CoefVar,output)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COLLECT ALL THE DATA FROM EVERYONE!

 

%%%%RESHAPE THE DATA FROM THE CELL COLLECTOR INTO LONG MATRICES THAT WE WANT! 
% WE HAVE COLLECTED ALL INTO OUTPUT
thetacollector=[];funceval=[];flc=[];fpc=[];surro_index=[];
for i=1:1:data.availcpus
dummy=cell2mat(output{i,:});
%size(dummy),
thetacollector=[thetacollector;dummy(:,1:data.Nth)];
%FIX!
funceval=[funceval;dummy(:,((data.Nth+1):(data.Nth+1)))];
flc=[flc;dummy(:,(end-2))];
fpc=[fpc;dummy(:,end-1)];
surro_index=[surro_index;dummy(:,end)];
end

runinfo.thetacoll_orig=thetacollector;
runinfo.fpc_orig=fpc;


dummyfull=[flc thetacollector funceval fpc surro_index];
dummyfull=unique(dummyfull,'rows');
flc=dummyfull(:,1),
thetacollector=dummyfull(:,2:(data.Nth+1));
funceval=dummyfull(:,(2+data.Nth:1+data.Nth+data.Nobs));
fpc=dummyfull(:,(end-1));
surro_index=dummyfull(:,end);
runinfo.currentuniques(runinfo.gen)=length(flc),
pause(4);
runinfo.realpoints=[runinfo.realpoints;thetacollector(find(~surro_index),:)]; %theta values 
runinfo.realevals=[runinfo.realevals;funceval(find(~surro_index),:)]; %real evaluations of observables
runinfo.currentreal=funceval(find(~surro_index),:);
%size(runinfo.realpoints),size(runinfo.realevals),
runinfo.surrpoints=[runinfo.surrpoints;thetacollector(find(surro_index),:)]; %surrogate estimates
runinfo.current_surr=thetacollector(find(surro_index),:);


runinfo.combined=[runinfo.realpoints,runinfo.realevals]; %ola mazi
runinfo.surratio(runinfo.gen)=length(find(surro_index))/runinfo.currentuniques(runinfo.gen); %Ratio of Surrogate to real points per generation
            if(runinfo.gen>1)
            runinfo.acceptance(runinfo.gen)=runinfo.currentuniques(runinfo.gen)/data.Num(runinfo.gen);
            end
            
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PLOT!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (data.iplot)
            HigherDimPlotter(runinfo,data,thetacollector,flc);
            end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%EXPORT ALL TO TMCMC STRUCT%%%%%%%%%%%%%%
 TMCMC.NumSamples= data.Num(runinfo.gen);
      TMCMC.beta2=data.bbeta;
      TMCMC.tolCOV=data.tolCOV;
      TMCMC.unifbounds=data.unifbounds;
      TMCMC.samples=runinfo.thetacoll_orig;
      TMCMC.fj=runinfo.fpc_orig;
      TMCMC.meantheta0=runinfo.meantheta;
      TMCMC.SS=runinfo.SS;
      TMCMC.logEvidence=sum(logselection);
      TMCMC.p=runinfo.p;
      TMCMC.CoefVar=CoefVar;
      TMCMC.surratio=runinfo.surratio;
%save -ascii current_wca.dat runinfo.combined


function [logselection,CoefVar,pflag,data,runinfo]=...
    calculate_statistics(flc,data,thetacollector,CoefVar,logselection,runinfo)
fprintf('I AM CALCuLATING WEIGHTS,COV FOR THE TMCMC AND MAKING EVIDENCE CALCS FOR GEN %d\n',runinfo.gen);
j=runinfo.gen;
theta.j=thetacollector';
pflag=0;
%data.Num(j)=max(size(flc));  %Numbrer of Points
    [runinfo.p(runinfo.gen+1),CoefVar(runinfo.gen+1)]=fmincon(@Objlogp,runinfo.p(runinfo.gen),[],[],[],[],runinfo.p(runinfo.gen),2,...
                                [],data.options,flc,runinfo.p(runinfo.gen),data.tolCOV);
   
   %                         runinfo.p(runinfo.gen+1)=0.000001+2^(runinfo.gen)*runinfo.p(runinfo.gen);
    if (runinfo.p(runinfo.gen+1)>1) pflag=runinfo.p(runinfo.gen); runinfo.p(runinfo.gen+1)=1; 
        data.Num(runinfo.gen+1)=data.LastNum; %Since next generation is going to be the last stage, Set it to the large number
    end
    
    if (runinfo.p(runinfo.gen+1)>data.pstrict)
        data.Rsq=data.Rsqstrict;
    end
    % Compute weights and normalize
    % weight=fj(:,j).^(p(j+1)-p(j));
  
%%%%%%%%%%%%%%
   

fjmax=max(flc),    
weight=exp((flc-fjmax).*(runinfo.p(runinfo.gen+1)-runinfo.p(runinfo.gen)) );      
 %weight=exp((fj(:,j)-fjmax).*(p(j+1)-p(j)) );   
       

    runinfo.q=weight./sum(weight);
logselection(runinfo.gen)=log(sum(weight)/runinfo.currentuniques(runinfo.gen))+fjmax*(runinfo.p(runinfo.gen+1)-runinfo.p(runinfo.gen));
    % Store value of CoefVAr in order to check it against CoefVar=tolCOV
    CoefVar(runinfo.gen+1)=std(runinfo.q)/mean(runinfo.q);


%%%%%%PROPOSAL ADAPTION BASED ON previous runinfo.gen%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%BUILDING COVARIANCE MATRIX
    mean_of_theta=theta.j*runinfo.q;
    runinfo.meantheta(:,j)=mean_of_theta;
    
    meanv=kron( mean_of_theta, ones(1,length(runinfo.q)) );
    for i1=1:data.Nth
        for i2=i1:data.Nth
            runinfo.SS(i1,i2)=sum( runinfo.q'.*(theta.j(i1,:)-meanv(i1,:))...
                        .*(theta.j(i2,:)-meanv(i2,:)) );
            runinfo.SS(i2,i1)=runinfo.SS(i1,i2);
        end
    end
    %%%%%%%%%%%%DR PREPATION UPDATE%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (data.dodr)    
      %R=(R+R')/2;
      % R = bbeta*SS;
    [V,D] = eigs(data.bbeta*runinfo.SS);
    R=V*sqrt(D)*V'; 
    %R   = chol(bbeta*diag(diag(SS)));
          runinfo.R2 = R./data.drscale;     % second proposal for DR try    
          runinfo.iR = inv(R);
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    



  
  
function CoefVar=Objlogp(x,fj,pj,tol)
fjmax=max(fj);
%weight=fj.^(x-pj);
weight=exp( (fj-fjmax).*(x-pj) );
q=weight/sum(weight);
%CoefVar=std(q)/mean(q)-tol;
CoefVar=(std(q)/mean(q)-tol)^2;  




function []=HigherDimPlotter(runinfo,data,thetacollector,flc)

close all;
    fscreen=get(0,'ScreenSize');
        figure('Position',[0 -50 3*fscreen(3)/4 3*fscreen(4)/4]);
%determine how many slices we have
Ns=length(data.plotslices);
Nsf=Ns;
if(mod(Ns,2)~=0)
Nsf=Nsf+1;
end

if(data.do_surrogate)
    Nsf=Nsf*2; %2 plots per input;
end
Nrows=Nsf/2;
Nrows=Nrows+1; %1 extra row (2 columns) for surrratio and acceptace rate
i=0;
j=0;
while(i<Ns)
    j=2*j+1;
     i=i+1;
     length(data.plotslices{i}),
        %%%%%%%%%%%%%%%%%

    if (data.do_surrogate)
     
         if(length(data.plotslices{i})==2)  %A 2d slice
          ind1=data.plotslices{i}(1);
          ind2=data.plotslices{i}(2);
          ax(j)=subplot(Nrows,2,j);plot(ax(j),runinfo.realpoints(:,ind1),runinfo.realpoints(:,ind2),'b.',runinfo.current_surr(:,ind1),runinfo.current_surr(:,ind2),'ro');
          title(sprintf('Accumulated Real Evaluations and Current Surrogates At stage=%d with model %s',runinfo.gen,data.surrogate_type));
          xlabel(sprintf('x_%d',ind1));ylabel(sprintf('x_%d',ind2));grid on;
     
          x=thetacollector(:,ind1);
          y=thetacollector(:,ind2);
          ax(j+1)=subplot(Nrows,2,j+1);scatter(ax(j+1),x,y,30,flc,'filled');colorbar;
          xlabel(sprintf('x_%d',ind1));ylabel(sprintf('x_%d',ind2));grid on;axis square;title('Current Evaluation');
         end
     % linkaxes([ax(3) ax(2) ax(1)],'xy');
         if(length(data.plotslices{i})==3)  %A 3d slice
             ind1=data.plotslices{i}(1);
             ind2=data.plotslices{i}(2);
             ind3=data.plotslices{i}(3);
     ax(j)=subplot(Nrows,2,j);plot3(ax(j),runinfo.realpoints(:,ind1),runinfo.realpoints(:,ind2),runinfo.realpoints(:,ind3),'b.',runinfo.current_surr(:,ind1),runinfo.current_surr(:,ind2),runinfo.current_surr(:,ind3),'ro');
    title(sprintf('Accumulated Real Evaluations and Current Surrogates At stage=%d with model %s',runinfo.gen,data.surrogate_type));
      xlabel(sprintf('xa_%d',ind1));ylabel(sprintf('x_%d',ind2));zlabel(sprintf('x_%d',ind3));grid on;
  
      x=thetacollector(:,ind1);
      y=thetacollector(:,ind2);
      z=thetacollector(:,ind3);  
      ax(j+1)=subplot(Nrows,2,j+1); scatter3(ax(j+1),x,y,z,20,flc,'filled');title(sprintf('Current Evaluation at stage %d',runinfo.gen));colorbar;
    %legend('Accumulated real function evaluation','Current Generation Surrogate evaluations',0);
      xlabel(sprintf('x_%d',ind1));ylabel(sprintf('x_%d',ind2));zlabel(sprintf('x_%d',ind3));grid on;
          end
    
          else    %if NO surrogates are used
              
            if(length(data.plotslices{i})==2)  %A 2d slice
          ind1=data.plotslices{i}(1);
             ind2=data.plotslices{i}(2);
        x=thetacollector(:,ind1);
        y=thetacollector(:,ind2);
        ax(i)=subplot(Nrows,2,i);scatter(ax(i),x,y,30,flc,'filled');colorbar;
       xlabel(sprintf('x_%d',ind1));ylabel(sprintf('x_%d',ind2));grid on;axis square;title('Current Evaluation');
   
     

            end
           if(length(data.plotslices{i})==3)  %A 3d slice
             ind1=data.plotslices{i}(1);
             ind2=data.plotslices{i}(2);
             ind3=data.plotslices{i}(3);
               x=thetacollector(:,ind1);
               y=thetacollector(:,ind2);
               z=thetacollector(:,ind3);  
      ax(i)=subplot(Nrows,2,i); scatter3(ax(i),x,y,z,20,flc,'filled');title(sprintf('Current Evaluation at stage %d',runinfo.gen));colorbar;
     xlabel(sprintf('x_%d',ind1));ylabel(sprintf('x_%d',ind2));zlabel(sprintf('x_%d',ind3));grid on;axis square;
           end
    end

     
end
         
if (data.do_surrogate)
ax(j+2)=subplot(Nrows,2,j+2);plot(ax(j+2),1:runinfo.gen,runinfo.surratio,'*-');
xlabel('Stage/gen');ylabel('Surrogate ratio');ylim([0 1]);grid on;
ax(j+3)=subplot(Nrows,2,j+3);plot(ax(j+3),1:runinfo.gen,runinfo.acceptance,'k*-');
xlabel('Stage/gen');ylabel('Acceptance rates');ylim([0 1]);grid on;
    pause(1);
       if (data.iprint)
        plot_name=['tmcmc_mc_',num2str(runinfo.gen),'.png'];
         print('-dpng',plot_name);
       end
else
%   ax(i+1)=subplot(Nrows,2,i+1);plot(ax(i+1),1:runinfo.gen,runinfo.surratio,'*-');
% xlabel('Stage/gen');ylabel('Surrogate ratio');ylim([0 1]);grid on;
ax(i+1)=subplot(Nrows,2,i+1);plot(ax(i+1),1:runinfo.gen,runinfo.acceptance,'k*-');
xlabel('Stage/gen');ylabel('Acceptance rates');ylim([0 1]);grid on;  
 if (data.iprint)
        plot_name=['tmcmc_mc_',num2str(runinfo.gen),'.png'];
         print('-dpng',plot_name);
       end
end
