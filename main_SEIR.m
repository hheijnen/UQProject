%% Flu epidemic

close all 
clear all

load seasonaldata.mat;

data = season1;

burnsamples = 1000;
runsamples = 1000;

%% Prior calculation
% optimal value
alpha = 0.13;
tau = 1.8;
gamma = 1;
m = 2.7*1e-6;
S0 = 0.3;
E0 = 0.05;
beta = 50;

prior(1) = alpha;
prior(2) = tau;
prior(3) = gamma;
prior(4) = m;
prior(5) = S0;
prior(6) = E0;

opt_prior = optimize_prior_SEIR(data,prior);

upperbounds(1) = 0.7;
upperbounds(2) = 5;
upperbounds(3) = 2;
upperbounds(4) = 1e-4;
upperbounds(5) = .9;
upperbounds(6) = .9;

for i = 1:6
    if (opt_prior(i) <= 0) || (opt_prior(i) >= upperbounds(i))
        opt_prior(i) = prior(i);
    end
end

parameters = {
    {'alpha',   opt_prior(1),    0,   upperbounds(1)       }
    {'tau',     opt_prior(2),    0,   upperbounds(2)       }
    {'gamma',   opt_prior(3),    0,   upperbounds(3)       }
    {'m',       opt_prior(4),    0,   upperbounds(4)       }
    {'S0',      opt_prior(5),    0,   upperbounds(5)       }
    {'E0',      opt_prior(6),    0,   upperbounds(6)       }
};

%% Sampling

postSEIR(data,parameters,burnsamples,runsamples);

%% Most probable values



%% Calculate Hessian

