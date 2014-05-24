function postSIR(data)

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
    {'tau',     1.5,    1, 2        }
    {'m',       0.2e-6, 0, 0.4e-6   }
    {'S0',      0.5,    0, 1        }
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