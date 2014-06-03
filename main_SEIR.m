%% Flu epidemic

close all 
clear all

load seasonaldata.mat;
load prior.mat;

data = season1;

burnsamples = 1000;
runsamples = 1000;

%% Prior calculation
%optimal value
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

% upperbounds(1) = 0.7;
% upperbounds(2) = 5;
% upperbounds(3) = 2;
% upperbounds(4) = 1e-4;
% upperbounds(5) = .9;
% upperbounds(6) = .9;

parameters = {
    {'alpha',   opt_prior(1),    0,          }
    {'tau',     opt_prior(2),    0,          }
    {'gamma',   opt_prior(3),    0,          }
    {'m',       opt_prior(4),    0,          }
    {'S0',      opt_prior(5),    0,          }
    {'E0',      opt_prior(6),    0,          }
};

%% Sampling

postSEIR(data,parameters,burnsamples,runsamples);

%% Most probable values



%% Calculate Hessian

