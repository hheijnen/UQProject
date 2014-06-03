%% Flu epidemic
clc;
close all 
clear all

load seasonaldata.mat;

data = season1;

burnsamples = 2000;
runsamples =  2000;

%% Prior calculation
% optimal value
alpha = 0.13;
tau = 1.4;
m   = 2e-6;
S0  = 0.7;

prior(1) = alpha;
prior(2) = tau;
prior(3) = m;
prior(4) = S0;


opt_prior = optimize_prior_SIR(data,prior);

upperbounds(1) = 0.7;
upperbounds(2) = 5;
upperbounds(3) = 1e-4;
upperbounds(4) = .9;

for i = 1:4
    if (opt_prior(i) <= 0) || (opt_prior(i) >= upperbounds(i))
        opt_prior(i) = prior(i);
    end
end

parameters = {
    {'alpha',   opt_prior(1),    0,   upperbounds(1)       }
    {'tau',     opt_prior(2),    0,   upperbounds(2)       }
    {'m',       opt_prior(3),    0,   upperbounds(3)       }
    {'S0',      opt_prior(4),    0,   upperbounds(4)       }
};

%% Sampling

[result, chain, s2chain] = postSIR(data,parameters,burnsamples,runsamples,@SIRss);

%% Plot the posterior
figure
mcmcplot(chain,[],result,'denspanel',2)

%% mean and std of the params
chainstats(chain,result);

%% optimal value
load debugdata.mat
[dummy, index] = min(ssmat);
index
format long;
optpars = thetamat(:,index)

time = 1:length(data);

YModelOpt = SIRfun(time,optpars,[optpars(4) data(1) 1-data(1)-optpars(4)]);
YModelMean = SIRfun(time,mean(chain),[optpars(4) data(1) 1-data(1)-optpars(4)]);
plot(time,[YModelOpt(:,2),data]);
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
    y0 = [chain(idx,4) data(1) 1-data(1)-chain(idx,4)];
    Y = SIRfun(time,chain(idx,:),y0);
    II(i,:) = Y(:,2)';
    %plot(time,Y(:,2));
    %hold on;
end
Q95 = quantile(II,0.975);
Q05 = quantile(II,0.025);


%% plot the stuff
XX= [time,fliplr(time)]';
YY = [Q05,fliplr(Q95)];
h = fill(XX,YY,[0.9 0.9 0.9]);
set(h,'EdgeColor','None');
hold on;
plot(time,YModelOpt(:,2),'b','linewidth',2);
plot(time,data,'kd','linewidth',2)
xlabel('Time (weeks)');
ylabel('I');
legend('95% Confidence Interval','Optimal fit','Data');
title('Flu-season 12/13, SIR')
%% Calculate Hessian

[parF2, fval2, exitflag, output, grad, hessian] = fminunc(@(par) diff_sqr(par), optpars);
hessian

%% plot the data for the different seasons ('FluData.eps')
plot(seasonaldata,'d','linewidth',2);
xlabel('Time (Weeks)');
ylabel('I');
legend('Winter 10/11','Winter 11/12','Winter 12/13');
title('Google Flu Trends Belgium');

