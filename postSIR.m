function postSIR(data)
clear model params options
close all

%% Model

time = [1:length(data)];
y = data;


model.ssfun = @SIRss;

%% Parameters
% PARAMS  theta structure
%   {  {'par1',initial, min, max, pri_mu, pri_sig, targetflag, localflag}
%      {'par2',initial, min, max, pri_mu, pri_sig, targetflag, localflag}
%      ... }
parameters = {
    {'alpha',   0.2,    0, 0.4      }
    {'tau',     0.1,    0, 2        }
    {'m',       0.2e-6, 0, 0.4e-6   }
    {'S0',      0.4,    0, 1        }
%     {'B0',       26,    1, 52       }
%     {'B1',       20,    0, 26       }
%     {'B2',        1,    0, 1        }
    };
%% Run
%

% burn 2000
options.nsimu = 2000;
[results, chain, s2chain]= mcmcrun(model,data,parameters,options);
%%
% Then re-run starting from the results of the previous run.

options.nsimu = 2000;
[results, chain, s2chain] = mcmcrun(model,data,parameters,options, results);

%% Plot

figure
mcmcplot(chain,[],results,'pairs');
figure
mcmcplot(chain,[],results,'denspanel',2);