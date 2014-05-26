function postSEIR(data)
 clear model parameters options results chain s2chain ssmat thetamat
close all

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

%    {'thres',   9   ,    0, 10      }
%     {'B0',       26,    1, 52       }
%     {'B1',       20,    0, 26       }
%     {'B2',        1,    0, 1        }
%% Run
%

%  for(i = 1:100)
options.nsimu = 100;
%%
% Then re-run starting from the results of the previous run.
% load debugdata.mat
% 
% parameters = {
%     {'alpha',   rand*0.4,    0, 0.4      }
%     {'tau',     rand*2,    0, 2       }
%     {'gamma',   rand,    0,   1    }
%     {'m',       rand*4e-06, 0,   6e-6    }
%     {'S0',      0.8*rand,    0, 1       }
%     {'E0',      0.1*rand,    0, 1       }
% };

parameters = {
    {'alpha',   0.4,    0, 0.7      }
    {'tau',     4.34,    0, 5       }
    {'gamma',   1,    0,   2   }
    {'m',       8e-06, 0,   1e-5    }
    {'S0',      0.343,    0, 1       }
    {'E0',      0.0056,    0, 1       }
};
  
[results, chain, s2chain]= mcmcrun(model,data,parameters,options);


% options.nsimu = 2000;
% 
% load debugdata.mat
% 
% ssmat(ssmat < 0 ) = 100;
% 
% alphai =    thetamat(1,min(find(min(ssmat) == ssmat)))
% taui =      thetamat(2,min(find(min(ssmat) == ssmat)))
% gammai =    thetamat(3,min(find(min(ssmat) == ssmat)))
% mi =        thetamat(4,min(find(min(ssmat) == ssmat)))
% S0i =       thetamat(5,min(find(min(ssmat) == ssmat)))
% E0i =       thetamat(6,min(find(min(ssmat) == ssmat)))
% 
% parameters = {
%     {'alpha',   alphai,    0, 0.4      }
%     {'tau',     taui,    0, 2       }
%     {'gamma',   gammai,    0,   1    }
%     {'m',       mi, 0,   4e-6    }
%     {'S0',      S0i,    0, 1       }
%     {'E0',      E0i,    0, 1       }
% };

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
