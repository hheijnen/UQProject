function [results, chain, s2chain] = postSEIR(data,parameters,burnsamples,runsamples)

%% Model

time = [1:length(data)];
y = data;

load debugdata.mat

ssmat = [];
thetamat = [];

save debugdata.mat ssmat thetamat
clear ssmat thetamat

model.ssfun = @SEIRss;

%% Parameters
% PARAMS  theta structure
%   {  {'par1',initial, min, max, pri_mu, pri_sig, targetflag, localflag}
%      {'par2',initial, min, max, pri_mu, pri_sig, targetflag, localflag}
%      ... }

% parameters = {
%     {'alpha',   1,    0,          }
%     {'tau',     1,    0,          }
%     {'gamma',   1,    0,          }
%     {'m',       1,    0,          }
%     {'S0',      1,    0,          }
%     {'E0',      1,    0,          }
% };
%% Run
%

%  for(i = 1:100)
options.nsimu = burnsamples;
  
[results, chain, s2chain]= mcmcrun(model,data,parameters,options);

options.nsimu = runsamples;

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
% I0 = mean(data(1:3,1));
% 
% modelfun = @(data,theta) SIRfun(time,theta,[theta(4) mean(data(1:3,1)) 1-theta(4)-mean(data(1:3,1))]');
% 
% nsample = 500;
% out = mcmcpred(results,chain,s2chain,time,modelfun,nsample);
% figure
% mcmcpredplot(out);
% % add the 'y' observations to the plot
%   hold on
%   plot(time,data(:,2),'s'); 
% %   ylabel(''); title(data.ylabels(i+1));
%   hold off
% xlabel('weeks');
