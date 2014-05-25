%Code Strategy:
% UQ_2014 Students: you can play with the parameters here to adjust them to
% your problem
%FIX THIS! This points to the directory you have the code
%data.masterpath='/Users/pa49/Desktop/mcmc_r0_1/tmcmc_distrib/'
data.materpath='/Users/stefan/Dropbox/ETHZDropbox/7Sem/Uncertainty Quantification/UQProject/tmcmc_distrib/'

%%%%DACE -package (distributed with the code).
%addpath('/Users/pa49/Desktop/mcmc_r0_1/tmcmc_distrib/dace/');


%What is the dimension of the parameter space?
data.Nth=4;
%%HOW MANY workers to spawn (set as many as cores in a multicore machine
%%e.g single core machine=1, mac laptop= 2, i7=8,brutus old node = 16, new node=48 and ppl screaming at you for their licences)
data.availcpus=2;
%%%MAX STAGES (Ignore)
data.maxstages=5000;
%How many Samples per stage are we going to use? Play with this!
data.Num=100*ones(data.maxstages,1); 
%Last stage, Generate this many samples
data.LastNum=500;
%Covariance tolerance for weights normalization  (Set to values between 0.1
%and 1). 0.1= slow transition , 1= fast transition
data.tolCOV=0.5;
%Scaling of the Covariance Matrix good for the proposal (Ching proposes 0.2)
data.bbeta=0.2;

%FEASIBLE PARAMETER SPACE
data.unifbounds=[0*ones(1,data.Nth);3*ones(1,data.Nth) ];



data.model='uq14_model'; % GOT TO THE FILE AND ADAPTH IT!!!!! your fucntion (likelihood in the uq13_model)


data.tournsize=8;



%%%FIX THE PRIORS! EACH LINE IN THE CELL CORRESPONDS TO ONE OF THE VARS
data.priorf{1}='Uniform';data.priorparam{1,1}=0;data.priorparam{1,2}=0.4; %i.e Gaussian with mean 0.2 std 10.7
data.priorf{2}='Uniform';data.priorparam{2,1}=0;data.priorparam{2,2}=3; %i.e Gaussian with mean -0.2 std 20.7
data.priorf{3}='Uniform';data.priorparam{3,1}=0;data.priorparam{3,2}=4e-6; %i.e Gaussian with mean 0.2 std 10.7
data.priorf{4}='Uniform';data.priorparam{4,1}=0;data.priorparam{4,2}=1; %i.e Gaussian with mean 0.2 std 10.7
%data.priorf{5}='Uniform';data.priorparam{5,1}=0;data.priorparam{5,2}=10000; %i.e Gaussian with mean 0.2 std 10.7




%IMPORTANT! 
%READ AND PASS HERE THE DATA TO PASS TO YOUR LIKELIHOOD!
%Values enter via the structure data
%e.g  : 12 experimental observations with their corresponding uncertainty
% data.expinputs=[0.0619    0.1120    0.0320    0.0484    0.0337    0.0358
% 0.0456    0.08300    0.1560    0.1780    0.2120    0.2490];
% data.expvar=[0.000183 0.00066 0.000096 0.00014 0.00011 0.00010 0.0001 0.001 0.002 0.002  0.002 0.003];
load seasonaldata.mat
data.expvar = season1;


%NOT IMPORTANT BUT CAN BE INSTRUCTIVE
%PLOTTING FOR ND>3 DATA:
%GUIDE: PUT THE SLICES IN THIS FORMAT : data.plotslices={[1 3] [1 2 3]}  %means
% make a 2-window subplot where we have a 2-D slice (x_1 with x_3) and a 3d
% slice x_1 with x_2 and x_3  .Can be any number of slices;
data.plotslices={ [1 2] [1 2 3] [2 3] [1 2 4] [2 3 4] };
data.Nobs=1;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               DO NOT TOUCH!                                             %
%%%%%%%%%%%%%%%%%%%%%%UQ-SPRING_13%%%%IGNORE_FROM_NOW_ON%%%%%%%%%%%%%%%%%%%
%Likelihood (Stupid Nomenclature-DO NOT TOUCH!!!!!!!!!)
data.likelihood='likelihood';
%Prior=
data.prior='priorpdf';
data.throw_surrmax=1;
%METHOD PICK: 1 for NORMAL TMCMC, 2: For DRAM-TMCMC
data.decidemethod=1;
%%%%%%%%%%%%RESTART WITH ALL THE SAME!%%%%%
data.loaddata=0;
data.restart=0;
%%%%%%%%%%%%%APPEND MORE INFO%%%%%%%%%%%%%%%%
%%%%%%%%REFINE%%%%%%%%%%%%%%
data.refine=0;
%if totalrestart together with refine, then we start from the beggining,
%giving a matrix containing all the last inputs.
data.totalrestart=0;
%Do we do Delayed rejection?
data.dodr=0;
%What is the prposal scale for delayed rejection?
data.drscale=3;  
%Do we use surrogate estimates?
data.do_surrogate=0;
%What type of surrogate estimation do we do?
%FOR DIM>3 DO NOT USE quadratic 
data.surrdim=2;
data.expcor=0;
data.surrogate_type='radialbasisf';   %Alternative - 'pseudospectral'
data.surrogate_type='nearestneighbor'
data.surrogate_type='quadratic';
data.surrogate_type='radialbasisf';   %Alternative - 'pseudospectral'
data.surrogate_type='krig';   %Alternative - 'pseudospectal'

%Number of Neighboors used to decide if we are going to do a Surrogate
data.nneighbors=12; %51;/?
%Normalized Distance threshold for neighbor identification in surrogate.
%This is the amplification factor that has the shape of the contour of
%covariance matrix, times this value
data.threshold=10;

%%%%%%IF WE USE RBF, TAKE A LOT
if (strcmp(data.surrogate_type,'radialbasisf'))
   data.nneighbors=12;
   data.threshold=60;
end

%%%%%%IF WE USE KRIG, TAKE ALMOST EVERYTHING
if (strcmp(data.surrogate_type,'krig'))
   data.nneighbors=15;
   data.threshold=9;
end


%%%%%%IF WE USE NEAREST NEIGHBOR, TAKE A LOT CAUSE IT WILL AUTOMATICALLY
%%%%%%STOP FROM GETTING A NON_CONVEX HULL INTERPOLANT
if (strcmp(data.surrogate_type,'nearestneighbor'))
   data.nneighbors=28;
   data.threshold=60;
end
%Relaxed Balue of R.sq for initial fitting!
data.Rsqrelaxed=0.90;
%Original one-Start with relaxed
data.Rsq=data.Rsqrelaxed;
%R-squared hard limit to kill the surrogate evaluation
data.Rsqstrict=0.95;
%At which p do we change to the strict value?
data.pstrict=0.2;
%Are we going to create only a Latin Hypercube of points?
data.Hypercube=0;

if (strcmp(data.surrogate_type,'krig'))
data.Rsqrelaxed=1e-3;
%Original one-Start with relaxed
data.Rsq=data.Rsqrelaxed;
%R-squared hard limit to kill the surrogate evaluation
data.Rsqstrict=1e-4;
end
%Level of suroogates fit
data.sumofsquare=0;
%%%NOBS%%%%%
%data.Nobs=113;   %Number of observations for surrogates
if (data.sumofsquare==1)
    data.Nobs=1;
end
%%%PLOT OUTPUT OR NOT%%%%%
data.iplot=1;
data.iprint=0;
%%%%%%%%
%%%Function names

%data.model='uq_rdf';
%data.model='testing_grid';

%How many model error dimensions do we have? (This is needed for the
%surrogate estimate, so to avoid interpolating there.
data.modelerrordim=1;
data.group={1:1:data.Nobs};
if (data.modelerrordim>=2)
    data.model='uq_mdmultigroup';
    data.group={1:1:7 8:12} ;
end


%Sets for optimization for scaling p's finding
data.options=optimset('Display','iter','TolX',0.00001,'TolFun',0.00001,'Algorithm','interior-point');
