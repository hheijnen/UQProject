function [output]=tmcmc_worker(worker_index,data,work,runinfo)

i=worker_index;

worker_name=['worker_in_',num2str(i)];
fprintf('Hello from worker%s\n',worker_name);
v_surr=[];count=0;
%%%%%


%%%TELL WHAT I HAVE TO DO
fprintf('My homework is: i have to do %f chains\n\n',max(size(work)));

    
    
    
        if (data.Hypercube~=1)     %%%THEN WE ARE DOING TMCMC
            
          %identify how many chains do we have to do as a worker.
          NumChains_per_worker=work{1}(1);
           point_counter=1;
        for j=1:1:NumChains_per_worker
           %surr.vsurrcurrent=zeros(1,length(data.exptemps));surr.Rsq=1;
              theta_size=work{j,2}; %theta_size(1) length of chain
                theta=work{j,3};  %FIRST TIME LEADER
                thetaleader=theta;
                
                 if((runinfo.gen==1 && ~data.refine) || (runinfo.gen==1 && data.totalrestart && data.refine))  %at first generation/stage calculate the point likelihood
                
                    if (data.do_surrogate && data.refine)
                         [sur_outl,vsurrcurrent]=surrogate_estimate(runinfo,thetaleader,data);
                         if (sur_outl)
                         did_surrogate(1,j)=1;did_surrogate(2,j)=1;
                         fl(j,1)=feval('uq_tmcmc_surr',thetaleader,vsurrcurrent,data);%fl(j,nc+1),
                         fp(j,1)=feval('Posterior',thetaleader,data,data.prior,fl(j,1)); %posterior
 			 v_surr(1,:)=vsurrcurrent;v_surr(2,:)=vsurrcurrent;
                         else
                         [fl(j,1),vsurrcurrent]=feval(data.likelihood,thetaleader,data);    %Likelihood
                         fp(j,1)=feval('Posterior',thetaleader,data,data.prior,fl(j,1)); %posterior
                         v_surr(1,:)=vsurrcurrent;v_surr(2,:)=vsurrcurrent;did_surrogate(1,j)=0;did_surrogate(2,j)=0;
                         end
                    else
                        sur_outl=0;
                        [fl(j,1),vsurrcurrent]=feval(data.likelihood,thetaleader,data);    %Likelihood
                         fp(j,1)=feval('Posterior',thetaleader,data,data.prior,fl(j,1)); %posterior
                         v_surr(1,:)=vsurrcurrent;v_surr(2,:)=vsurrcurrent;did_surrogate(1,j)=0;did_surrogate(2,j)=0;

                    end
                
                
  
                
                else      %otherwise it is stored in work cell!
                    fl(j,1)=work{j,4};   %get the likelihood value stored for the leader
                    v_surr(1,:)=work{j,5};v_surr(2,:)= v_surr(1,:); %and store the observations of the leader
                    fp(j,1)=feval('Posterior',thetaleader,data,data.prior,fl(j,1)); %posterior of the leader
                    did_surrogate(1,j)=work{j,6}; did_surrogate(2,j)=did_surrogate(1,j);%whatever was the guy before;
                end
                %%%%EVALUATE FEASIBILITY OF SURROGATE ESTIMATE, IF FEASIBLE
                %%%%ESTIMATE IT  
                for nc=1:1:theta_size(1) %perform now x trials for this chain
                did_surrogate(nc+1,j)=did_surrogate(nc,j); %whatever was the guy before
                theta(nc+1,:)=thetaleader; % the next one is still there if not changed in the end
                fl(j,nc+1)=fl(j,nc); %Likelihood
                fp(j,nc+1)=fp(j,nc);
                v_surr(nc+1,:)=v_surr(nc,:);
                if(runinfo.gen>1 || (data.refine && ~data.totalrestart))
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MCMC-MH%%%%
                %PICKED A POINT AROUND thetac,MAKE A CHOICE AROUND IT WITH
                %THE ADAPTED PROPOSAL
                thetac=proprnd(thetaleader,data.bbeta*runinfo.SS,data.unifbounds);
                %%%%%%%%WHAT IS THE CANDIDATE EVALUATION?
              %%%%%%%%%%%%%CHECK ABOUT SURROGATES! IF GEN>2 AND A VARIOUS OF OTHER TESTS
              %INCLUDING IF IT LIES WITHIN THE CONVEX HULL OF ITS
              % SURROUNDING POINTS
                             if (data.do_surrogate)
                                [sur_out,vsurrcurrentc]=surrogate_estimate(runinfo,thetac,data);
                                if (sur_out)
                                did_surrogate(nc+1,j)=1;
                                flc=feval(@uq_tmcmc_surr,thetac,vsurrcurrentc,data);%fl(j,nc+1),
                                fpc=feval('Posterior',thetac,data,data.prior,flc); %posterior
                                else
                                did_surrogate(nc+1,j)=0;
                                [flc,vsurrcurrentc]=feval(data.likelihood,thetac,data); %proposal likelihood
                                fpc=feval('Posterior',thetac,data,data.prior,flc); %proposal posterior
                                end
                             else
                               sur_out=0;
                               did_surrogate(nc+1,j)=0;
                               [flc,vsurrcurrentc]=feval(data.likelihood,thetac,data); %PROPOSED POST
                               fpc=feval('Posterior',thetac,data,data.prior,flc); %posterior
                            end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %DECIDE%%%%%%%%%%%%%%%%%%%%
                    if (data.decidemethod==1)
                 r=exp( (fpc-fp(j,nc+1))*runinfo.p(runinfo.gen));  %New-old
                    else
                        r=exp( (fpc-fp(j,nc)));
            % change for antehr method 
                    end
                 r=min(1,r);
                state=find( mnrnd(1,[r,1-r]) );        
                        if (state==1)
                         %if the new point is better, update its
                         %information of the current one!
                         theta(nc+1,:)=thetac; %will send it to the output 
                         thetaleader=thetac;   %the leader is the new guy
                        %also keep the value of probability 
                         fl(j,nc+1)=flc;
                         fp(j,nc+1)=fpc;
                         count=count+1;
                         v_surr(nc,:)=vsurrcurrentc; 
                         v_surr(nc+1,:)=vsurrcurrentc;
                       
                        
                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%DELAYED REJECTION%%%%%%%%%%%%%%%%%%%%%%%%%
                    elseif (state==2 && data.dodr)
            %%%WE FAILED 1 st time, SO WE CONTINUE WITH DR
            %%%%%DR PIECEs%%%%%%%%%%%%%%%%%%%%%%
            %GET A NEW SAMPLE , with a different distribution made of the
            %covariance matrix (adapted)

            %NOTE : The ideal 1-step delayed rejection involved a sqrt(2)/3
            %scaling of the proposal which is done in the worker (Green,
            %Haario)
                     thetac2=proprnd(thetaleader,data.bbeta*runinfo.R2,data.unifbounds);  
        % compute the likelihood for the new candidate sample
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         We picked a NEW point based on the ORIGINAL leader, lets go all over again....
%         REMEMBER: the proposal is now scaled dy sqrt(2)/data.drscale!
                  if (data.do_surrogate)
                                [sur_out,vsurrcurrentc2]=surrogate_estimate(runinfo,thetac2,data); %NOTE! the search radius is still the big one
                                if (sur_out)
                                did_surrogate(nc+1,j)=1;
                                flc2=feval('uq_tmcmc_surr',thetac2,vsurrcurrentc2,data);
                                fpc2=feval('Posterior',thetac2,data,data.prior,flc2); %posterior
                                else
                                did_surrogate(nc+1,j)=0;
                                [flc2,vsurrcurrentc2]=feval(data.likelihood,thetac2,data); %proposal likelihood
                                fpc2=feval('Posterior',thetac2,data,data.prior,flc2); %proposal posterior
                                end
                  else
                               sur_out=0;
                               did_surrogate(nc+1,j)=0;
                               [flc2,vsurrcurrentc2]=feval(data.likelihood,thetac2,data); %PROPOSED POST
                               fpc2=feval('Posterior',thetac2,data,data.prior,flc2); %posterior
                            end
    %Acceptance/rejection ratio
    %TO PRWTO EINAI TO fj(k,j) 
    
%ADAPTED KOMMATI
                                  r32=min(1,exp( (fpc-fpc2)));    %part with previously drawn canditate
                                  l2     = exp(fpc2-fp(j,nc+1));  %part with leader
                                  q1      = exp(-0.5*norm((thetac2-thetac)*runinfo.iR)^2-norm((thetaleader-thetac)*runinfo.iR)^2); %second proposal | 1st
                                  r13 = l2*q1*(1-r32)/(1-r);
 
                                 state2=find( mnrnd(1,[r13,1-r13]) );       
                                  if(state2==1)    %IF DR gets accepted
                                  theta(nc+1,:)=thetac2; %will send it to the output 
                                  thetaleader=thetac2;   %the leader is the new dr guy
                                 %also keep the value of probability 
                                  fl(j,nc+1)=flc2;
                                  fp(j,nc+1)=fpc2;
                                  count=count+1;
                                  v_surr(nc,:)=vsurrcurrentc2; 
                                  v_surr(nc+1,:)=vsurrcurrentc2;
                                  %fprintf('Another one bites the DRed')
                                  else  %everything failed...even DR!
                                      if (sur_out==1) did_surrogate(nc+1,j)=did_surrogate(nc,j);
                                      end %if the surrogate did not get selected, then the new one tag is going to what was the first one (if it was a real one)
                                  end     %DR FAILED-NO MORE TRIALS
                                     
                        else      %If not DR and we failed at the first step...
                            if (sur_out==1) did_surrogate(nc+1,j)=did_surrogate(nc,j);end %if the surrogate did not get selected, then the new one tag is going to what was the first one (if it was a real one)
                        
                            
                        end  %of trials (of if state==1)
                        

                
                end  %of runinfo.runinfo.gen>1
                           
                           output(point_counter,:)={theta(nc+1,:),v_surr(nc+1,:),fl(j,nc+1),fp(j,nc+1),did_surrogate(nc+1,j)};
                           point_counter=point_counter+1;
                           sur_out=0;sur_outl=0;
            end  %of number of chain length (nc)
        end 
    end
    
    
    
    
    
    
  
    %%%%%%%%%%%WE ARE FINISHED WITH THE SAMPLING OF THE POINTS, LETS SHIP
    %IT BACK TO THE DISTRIB GUY DIRECTORY
   %copy_files_tomaster(i,data,output);







%%%%%%%%%%%%%%%%%%%%%%%%%%HELPER EVALUATING FUNCTIONS&&&&&&&&&&&&&&&&&&&&&&




function f=Posterior(theta,data,priorpdf,LH)

Prior=feval(priorpdf,theta,data);
%f=LH*Prior;
f=LH+log(Prior);



%%%%%%CALL THE MODEL LIKELIHOOD FUNCTION
function [f,output]=likelihood(theta,data)
%uq_tmcmc returns loglikelihood
%returns log-likelihood
[f,output]=feval(data.model,theta,data);


function thetac=proprnd(theta,S,bounds)
% Likelihood for Random Walk problem
id=1;
while id==1
    thetac=mvnrnd(theta,S);
    if (length(theta)==1)
        if (thetac<bounds(1) || thetac>bounds(2))
        id=1;
        else
        id=0;
        end
        
    else
    % exclude samples that fall outside the range of variation of theta
        if (any(thetac<bounds(1,:)) || any(thetac>bounds(2,:)))
            id=1;
        else
        id=0;
        end
    end
end




function value = priorpdf(theta,data)
%
%size(theta),size(data.unifbounds(1,:)),
if(size(theta)~=size(data.unifbounds(1,:)));
    theta=theta';
end

%create first the input vector that we reckon about each dimension
%unifpdf(theta,data.unifbounds(1,:),data.unifbounds(2,:))
for i=1:1:data.Nth
    if (length(data.priorparam(i,:))==2)
pv(i)=pdf(data.priorf{i},theta(i),data.priorparam{i,1},data.priorparam{i,2});
    else
pv(i)=pdf(data.priorf{i},theta(i),data.priorparam{i});
    end
end
value=prod(pv);




%BAD EXAMPLE OF CONSTRUSTION OF LIKELIHOOD
 function [fval,output]=uq_md(x,data)
%Load the observations (n) 


%The covariance matrix of the errors.


if (data.sumofsquare)
Fapp=hbs(x(1),x(2),data.exptemps,data.exppress);  %krisimo shmeio run (couple gromacs)
U=x(end)*eye(data.Nobs,data.Nobs);

gausexp=diag(U).^(-2)'.*((data.expinputs-Fapp)./data.expinputs).^2;
output=sum(gausexp);

fval=(-0.5*(data.Nobs*log(2*pi)+log(det(U)))-0.5*sum(gausexp));
else
    
U=x(end)*eye(data.Nobs,data.Nobs);
Fapp=hbs(x(1),x(2),data.exptemps,data.exppress);  %krisimo shmeio run (couple gromacs)
output=Fapp;

gausexp=diag(U).^(-2)'.*((data.expinputs-Fapp)./data.expinputs).^2;
fval=(-0.5*(data.Nobs*log(2*pi)+log(det(U)))-0.5*sum(gausexp));
end





  
  
%DO NOT LOOK AT ME I BARELY EXIST
%
function [did_surrogate,vsurrcurrent]=surrogate_estimate(runinfo,newpar,data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SURROGATE MODEL 1
%BEFORE we evaluate using the regular function, we check if the chain
%reached a maturity level for surrogate modelling
%this is done by checking how many neighboors a chair has given a fixed
%distance. If its more than a prespecified threshold we proceed to local
%surface fitting.
%Begin by assuming that the Neighboorhood search is going to fail
did_surrogate=0;
vsurrcurrent=zeros(1,data.Nobs);

%REAL INDEX SEARCH (SCREW SURROGATES FOR ESTIMATION)
size(runinfo.combined);
ready=unique(runinfo.combined,'rows');

reshchain=ready(:,1:(data.Nth-data.modelerrordim));
newpar=newpar(1:(data.Nth-data.modelerrordim));
obligatory_nearest=0;  %if the matrix is rank deficient-> pass it to another method (linear or 0th Krig)

%Also get the number of unique values now

%iu is now the index to these unique

% we reshape the chains to be equally balanced (a.k.a parameter space
% [0-1]x...x [0-1] )


%We want now to search around in an ellipsoid way and find the points
%within reach.This means rather than going with the vnorm , we construct an
%ellipse and check for points within them
[V,D]=eig(runinfo.SS(1:(data.Nth-data.modelerrordim),1:(data.Nth-data.modelerrordim)));  %eigendecomposition of the covariance matrix (only the relevant parameters,no model error ones)
 VV = V*D;  %initial scaling 
if data.surrdim==2
t = linspace(0,2*pi,100); 
 e = [cos(t) ; sin(t)];%# unit circle
 %VV = V*sqrt(D); %# scale eigenvectors (why?)
 %FIRST PASS - with original thresholding
 e = bsxfun(@plus, data.threshold*VV*e, newpar'); %#' project circle back to orig space,with center the new paramater (newpar)!    
 %# plot cov and major/minor axes
 %Now points e are the edges of our polygon. Use inpolygon now to get the
 %points within,scaled by a fixed point (data.threshold)
% plot(e(1,:), e(2,:), 'Color','k'); 
x=reshchain(:,1);
y=reshchain(:,2);
in = inpolygon(x,y,e(1,:),e(2,:));

% figure, plot(e(1,:),e(2,:),x(in),y(in),'r.',x(~in),y(~in),'b.');
%    pause(20);
nn=length(find(in)==1); %how many points are in the interior of the current eclipse

else
difs=bsxfun(@minus,newpar,reshchain);
[ineigh]=find(all(bsxfun(@lt,abs(difs),data.threshold*diag(sqrt(D))'),2));
nn=length(ineigh);


%rescale the points
%mS = mean(reshchain); sS = std(reshchain);
%for i=1:1:data.Nth, reshchain(:,i) = (reshchain(:,i) - mS(i))/sS(i); parresh(i)=(newpar(i)-mS(i))/sS(i);end
%dists=vnorm(reshchain(:,1:2)-bsxfun(@times,ones(up,2),parresh(1:2)),2);
%SORT THE NEAREST NEIGHBOORS
end

                if(nn>data.nneighbors)
                    %FIRST PASS -> if we fail here no need to try
                     if(data.surrdim==2)
                         [K] = convhulln([newpar(1:2);[x(in),y(in)]]);
                     else
                         [K] = convhulln([newpar;reshchain(ineigh,:)]);
                     end
                           if (~isempty(find(K==1, 1)))
                             did_surrogate=0;return; 
                           end
                    %We found some and we can continue cause now we have to
                    %find a minimum set around our point
                    %xin=x(in);
                    %yin=y(in);
                    blowfactor=0.1;
                    while (blowfactor<=1)
                        if (data.surrdim==2)
                            e = [cos(t) ; sin(t)];
                            e = bsxfun(@plus, (blowfactor*data.threshold)*VV*e, newpar(1:2)'); %#' project circle back to orig space,with center the new paramater (newpar)!    
                            in = inpolygon(x,y,e(1,:),e(2,:));
                            nnn=length(find(in)==1);
                        else
                            difs=bsxfun(@minus,newpar,reshchain);
                            [ineigh]=find(all(bsxfun(@lt,abs(difs),blowfactor*data.threshold*diag(sqrt(D))'),2));
                            nnn=length(ineigh);
                        end
                              blowfactor=blowfactor+0.2;
                              if(nnn>5)
                                  if(data.surrdim==2)
                                  [K] = convhulln([newpar(1:2);[x(in),y(in)]]);
                                  else
                                  [K] = convhulln([newpar;reshchain(ineigh,:)]);
                                  end
                              end                      
                              if (isempty(find(K==1, 1)) && nnn>=data.nneighbors) blowfactor=99999;end
                      end  %of while
%CONVEX HULL TESTING
                              obligatory_nearest=0;
                           if (nnn<data.nneighbors)
                               
                             did_surrogate=0;return; 
                           end
%%%%%NOW WE HAVE TO CHECK IF ITS A POINT TO AN EDGE! WE GET THE CONVEX HULL
%%%%%OF The CHOSEN ONES AND CHECK OUR LUCK :) 

%Then we proceed to do a surrogate model                
              did_surrogate=1;


              %do it one more time
      %              vnorm(reshchain(ineigh,1:2)-bsxfun(@times,ones(data.nneighbors,2),parresh(1:2)),2),
%[reshchain,iu,ju]=unique(runinfo.realpoints,'rows','first');
                       %vnorm(reshchain(ineigh,1:2)-bsxfun(@times,ones(data.nneighbors,2),newpar(1:2)),2),

                  
                    
                       
                       
             for  ti=1:1:data.Nobs

                  cleanres=ready(:,data.Nth+ti);
                          
%%%%%%%%%%%
%Now we have the clean nice data, we assign each to where it belongs
            
                            
                 if data.surrdim==2         
                   xdata=x(in);
                   ydata=y(in);
                   zdata=cleanres(in);
                  % [xdata,ydata,zdata];
                   indata=[xdata, ydata];
                 else
                    indata=reshchain(ineigh,:);
                        
                    zdata=cleanres(ineigh);
                 end

  
  
    %%%%%%NOW PICK METHOD%%%%%%%%%%%
    
    if (strcmp(data.surrogate_type,'krig'))
        
theta0 = 0.5*ones(1,data.surrdim); % An initial guess of  the correlation lengths.
lob = 0.01*ones(1,data.surrdim); % A lower bound on the correlation lengths.
upb = 1*ones(1,data.surrdim);  % An upper bound on the correlation lengths.
[dmodel,perf] = dacefit(indata,zdata,@regpoly2,@correxp,theta0,lob,upb);
          runinfo.krigpredictor=dmodel;

% Use the DACE model to interpolate the function at a point in the
% parameter space.
  
 
       
            f_dace = predictor(newpar,runinfo.krigpredictor);
            vsurrcurrent(ti)=f_dace;
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    elseif (strcmp(data.surrogate_type,'radialbasisf'))
     %radial Basis function approximation
      op=rbfcreate(indata',zdata','RBFFunction','Thinplate'); %can be changed to gaussian basis or etc etc

      f_rbf = rbfinterp(newpar',op);
                   vsurrcurrent(ti)=f_rbf;
  
    
    elseif (strcmp(data.surrogate_type,'quadratic') && obligatory_nearest==0)
        %  warning off;
     %%%%             QUADRATIC SURFACE FIT TO THE LOCAL POINTS
alpha = 0.05;
[coeffs,Ibeta,res,Ires,stats] = regress(zdata,[ones(size(xdata)) xdata ydata xdata.^2 xdata.*ydata ydata.^2],alpha);
     Rsq=stats(1);
     
 %%%%%%%%IF THE FIT IS BAD STOP AND FAIL THE SURROGATE ESTIMATION
     if (Rsq<data.Rsq)
         did_surrogate=0;
          return;
     end
    vsurrcurrent(ti)=quadratic(coeffs,[newpar(1),newpar(2)]);


        
    elseif (strcmp(data.surrogate_type,'nearestneighbor')  || obligatory_nearest==1)
           
    
    %%%%%MATLAB's NATURAL METHOD USING DELANUY INTERPOLATION using linear
    %%%%%fit

   f_nn=griddatan(indata,zdata,newpar);
    vsurrcurrent(ti)=f_nn;
    if (isnan(f_nn))
        did_surrogate=0;
        return;
    end
%      % Construct the interpolant
%          F = TriScatteredInterp(xdata,ydata,zdata,'nearest');
% 
%     % Evaluate the interpolant at the locations (qx, qy), qz
%     %    is the corresponding value at these locations.
%             f_nn = F(newpar(1),newpar(2));
%                    vsurrcurrent(ti)=f_nn;
    end

                end  %for Nobs
              
            
             else %of number of points with the initial try
did_surrogate=0;  
end            
      
  
  
function y = vnorm(A,varargin)
% VNORM - Return the vector norm along specified dimension of A
%
%   VNORM(A) returns the 2-norm along the first non-singleton
%   dimension of A
%   VNORM(A,dim) return the 2-norm along the dimension 'dim'
%   VNORM(A,dim,normtype) returns the norm specified by normtype
%   along the dimension 'dim'
%   VNORM(A,[],normtype) returns the norm specified by normtype along
%   the first non-singleton dimension of A
% 
%   normtype may be one of {inf,-inf,positive integer}.
%   For a given vector, v, these norms are defined as
%   inf: max(abs(v))
%   -inf: min(abs(v))
%   p (where p is a positive integer): sum(abs(v).^p)^(1/p)
%   
%   Examples:
%       A =  [8 1 6; 3 5 7; 4 -9 2];
%
%       %Columnwise 2-norm (Euclidean norm)
%       vnorm(A,1) = [9.4340 10.3441 9.4340];
%       vnorm(A,[],2) % Same as above (since first non-singleton dimensions
%                     % is columnwise and default norm is 2-norm.
%       vnorm(A,[],[])% Again, same as above
%
%       % Row-wise maximum of absolute values
%       vnorm(A,2,inf) = [8 7 9]'; 
%
%       % Columnwise minimum of absolute values
%       vnorm(A,[],-inf) = [3 1 2];
%       
%       % Error: Use the inf type and not the string 'inf'
%       vnorm(A,[],'inf')   % Wrong
%       vnorm(A,[],inf)     % Correct

dim = [];
ntype = [];

if nargin>1
    dim = varargin{1};
    if isempty(dim)
       idx = find(size(A)~=1);
       dim = idx(1);
    elseif dim~=floor(dim) || dim<1
        error('Dimension must be positive integer');
    end
    if nargin>2
        ntype = varargin{2};
    end
end


if isempty(ntype)
    y = sqrt(sum( abs(A).^2 , dim) );
elseif ntype==1
    y = sum( abs(A) , dim );    
elseif isinf(ntype)
    if ntype > 0
        y=max(abs(A), [], dim);
    else
        y=min(abs(A), [], dim);
    end
elseif ntype~=floor(ntype) || ntype<1
    error(['Norm type must be one of inf,-inf or a positive ' ...
           'integer']);
else 
    y = (sum( abs(A).^ntype , dim) ).^(1/ntype);
end 




function F = quadratic(a, data) 
%DATA TO BE FITTED FOR THE SURROGATE EVALUATION    
x = data(:, 1);
y = data(:, 2);
F = a(1) + a(2).*x + a(3).*y + a(4).*x.^2 + a(5).*x.*y + a(6).*y.^2;
