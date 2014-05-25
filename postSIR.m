function postSIR(data)
 clear model parameters options results chain s2chain ssmat thetamat
close all

%% Model

time = [1:length(data)];
y = data;

load debugdata.mat

ssmat = [];
thetamat = [];

save debugdata.mat ssmat thetamat

model.ssfun = @SIRss;

%% Parameters
% PARAMS  theta structure
%   {  {'par1',initial, min, max, pri_mu, pri_sig, targetflag, localflag}
%      {'par2',initial, min, max, pri_mu, pri_sig, targetflag, localflag}
%      ... }
parameters = {
    {'alpha',   0.2,    0, 0.4      }
    {'tau',     1.5,    0, 3        }
    {'m',       2e-06, 0,   4e-6    }
    {'S0',      0.2,    0, 1       }
%    {'thres',   9   ,    0, 10      }
%     {'B0',       26,    1, 52       }
%     {'B1',       20,    0, 26       }
%     {'B2',        1,    0, 1        }
    };
%% Run
%

% burn 2000
options.nsimu = 1000;
[results, chain, s2chain]= mcmcrun(model,data,parameters,options);
%%
% Then re-run starting from the results of the previous run.
% load debugdata.mat
% 
% parameters = {
%     {'alpha',   thetamat(1,find(min(ssmat) == ssmat)),    0, 0.4      }
%     {'tau',     thetamat(2,find(min(ssmat) == ssmat)),    0, 3        }
%     {'m',       thetamat(3,find(min(ssmat) == ssmat)),    0,   4e-6    }
%     {'S0',      thetamat(4,find(min(ssmat) == ssmat)),    0, 1       }
% }

options.nsimu = 1000;
[results, chain, s2chain] = mcmcrun(model,data,parameters,options,results);

%% Posterior data
chainstats(chain,results);

%% Plot

figure
mcmcplot(chain,[],results,'denspanel');
% figure
% mcmcplot(chain,[],results,'denspanel',2);

%% Pred plot
% We sample 500 parameter realizations from |chain| and |s2chain|
% and calculate the predictive plots.
I0 = mean(data(1:3,1));

modelfun = @(data,theta) SIRfun(time,theta,[theta(4) mean(data(1:3,1)) 1-theta(4)-mean(data(1:3,1))]');

nsample = 500;
out = mcmcpred(results,chain,s2chain,time,modelfun,nsample);
figure
mcmcpredplot(out);
% add the 'y' observations to the plot
  hold on
  plot(time,data(:,2),'s'); 
%   ylabel(''); title(data.ylabels(i+1));
  hold off
xlabel('weeks');
