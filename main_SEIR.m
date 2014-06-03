%% Flu epidemic

close all 
clear all

load seasonaldata.mat;
load prior.mat;

seasonindex = 3;

switch seasonindex
    case 1
        data = season1;
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
        
        parameters = {
            {'alpha',   opt_prior(1),    0,          }
            {'tau',     opt_prior(2),    0,          }
            {'gamma',   opt_prior(3),    0,          }
            {'m',       opt_prior(4),    0,          }
            {'S0',      opt_prior(5),    0,          }
            {'E0',      opt_prior(6),    0,          }
            };
        
    case 2
        data = season2;
        
        parameters = {
            {'alpha',   parF(1),    0,          }
            {'tau',     parF(2),    0,          }
            {'gamma',   parF(3),    0,          }
            {'m',       parF(4),    0,          }
            {'S0',      parF(5),    0,          }
            {'E0',      parF(6),    0,          }
            };
        
    case 3
        data = season3;
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
        
        parameters = {
            {'alpha',   opt_prior(1),    0,          }
            {'tau',     opt_prior(2),    0,          }
            {'gamma',   opt_prior(3),    0,          }
            {'m',       opt_prior(4),    0,          }
            {'S0',      opt_prior(5),    0,          }
            {'E0',      opt_prior(6),    0,          }
            };

end


burnsamples = 1000;
runsamples = 1000;

%% Sampling

postSEIR(data,parameters,burnsamples,runsamples);

%% Most probable values



%% Calculate Hessian

