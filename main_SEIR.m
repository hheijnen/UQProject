%% Flu epidemic

close all 
clear all

load seasonaldata.mat;

data = season1;

burnsamples = 1000;
runsamples = 1000;

%% Prior calculation
% optimal value
alpha = 0.1523;
tau = 1.7357;
gamma = 1;
m = 2.7*1e-6;
S0 = 0.2743;
E0 = 0.005;
beta = 50;

prior(1) = alpha;
prior(2) = tau;
prior(3) = gamma;
prior(4) = m;
prior(5) = S0;
prior(6) = E0;

opt_prior = optimize_prior_SEIR(data,prior);

%% Sampling

postSEIR(data,prior,burnsamples,runsamples);

%% Most probable values



%% Calculate Hessian

