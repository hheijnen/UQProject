function [results, chain, s2chain] = postSIR(data,parameters,burnsamples,runsamples,SSFunc)

%% Model

time = 1:length(data);

load debugdata.mat

ssmat = [];
thetamat = [];

save debugdata.mat ssmat thetamat
clear ssmat thetamat

model.ssfun = SSFunc;

%% Parameters
% PARAMS  theta structure
%   {  {'par1',initial, min, max, pri_mu, pri_sig, targetflag, localflag}
%      {'par2',initial, min, max, pri_mu, pri_sig, targetflag, localflag}
%      ... }

%    {'thres',   9   ,    0, 10      }
%     {'B0',       26,    1, 52       }
%     {'B1',       20,    0, 26       }
%     {'B2',        1,    0, 1        }
%% Run
%

%  for(i = 1:100)
options.nsimu = burnsamples;



  
[results, chain, s2chain]= mcmcrun(model,data,parameters,options);

options.nsimu = runsamples;

[results, chain, s2chain] = mcmcrun(model,data,parameters,options,results);

%% Posterior data
%chainstats(chain,results);

%% Plot

%figure
%mcmcplot(chain,[],results,'denspanel');
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
