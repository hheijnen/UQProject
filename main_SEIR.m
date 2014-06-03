%% Flu epidemic

close all 
clear all

load seasonaldata.mat;
load prior.mat;

seasonindex = 1;

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
            {'alpha',   parF2(1),    0,          }
            {'tau',     parF2(2),    0,          }
            {'gamma',   parF2(3),    0,          }
            {'m',       parF2(4),    0,          }
            {'S0',      parF2(5),    0,          }
            {'E0',      parF2(6),    0,          }
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


burnsamples = 500;
runsamples = 500;


%% Sampling

[result, chain, s2chain]  = postSEIR(data,parameters,burnsamples,runsamples);

%% Plot the posterior
figure
mcmcplot(chain,[],result,'denspanel',2)

%% Calculate the modes

findmodesSEIR

%% mean and std of the params
chainstats(chain,result);

%% optimal value
load debugdata.mat
[dummy, index] = min(ssmat);
index
format long;
optpars = thetamat(:,index)

time = 1:length(data);
S0 = optpars(5);
E0 = optpars(6);
I0 = data(1);
R0 = 1-S0-E0-I0;
y00 = [S0 E0 I0 R0];

figure
YModelOpt = SEIRfun(time,optpars,y00);
YModelMean = SEIRfun(time,mean(chain),y00);
plot(time,[YModelOpt(:,3),data]);
legend('Model','Data');
xlabel('Time(weeks)');
ylabel('I');
title('Comparison of the optimal fit for season 10/11');

%% bootstrap confidence intervals


% I0 = data(1);
% I = @(Results) Results(:,2); % wrap the SIRFun to return only I
% bootfun = @(params) I(SIRfun(time,params,[params(4) I0 1-params(4)-I0]));
% 
% nSamples = 10;
% 
% ci = bootci(nSamples,{bootfun,chain},'alpha',0.9);
% 
% errorbar(time,bootfun(optpars),ci(1,:),ci(2,:))
% hold on;
% plot(time,data,'g');


%% conf interval:  randomly choose 1e3 samples
N = 1000;
II = zeros(N,length(time));

for i=1:N
    idx = unidrnd(runsamples);
    
    S0 = chain(idx,5);
    E0 = chain(idx,6);
    I0 = data(1);
    R0 = 1-S0-E0-I0;
    y0 = [S0 E0 I0 R0];
    
   % y0 = [chain(idx,4) data(1) 1-data(1)-chain(idx,4)];
    Y = SEIRfun(time,chain(idx,:),y0);
    II(i,:) = Y(:,3)';
    %plot(time,Y(:,2));
    %hold on;
end
Q95 = quantile(II,0.975);
Q05 = quantile(II,0.025);


%% plot the stuff
figure
XX= [time,fliplr(time)]';
YY = [Q05,fliplr(Q95)];
h = fill(XX,YY,[0.9 0.9 0.9]);
set(h,'EdgeColor','None');
hold on;
plot(time,YModelOpt(:,3),'b','linewidth',2);
plot(time,data,'kd','linewidth',2)
xlabel('Time (weeks)');
ylabel('I');
legend('95% Confidence Interval','Optimal fit','Data');
title('Flu-season 11/12, SEIR')
%% Calculate Hessian

[parF2, fval2, exitflag, output, grad, hessian] = fminunc(@(par) diff_sqr(par), optpars);
hessian

%% plot the data for the different seasons ('FluData.eps')
plot(seasonaldata,'d','linewidth',2);
xlabel('Time (Weeks)');
ylabel('I');
legend('Winter 10/11','Winter 11/12','Winter 12/13');
title('Google Flu Trends Belgium');

