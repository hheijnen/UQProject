%% script which tries to plot the uncertainty around the estimated I
% by using the bootstrap method

%use the chain data from the mcmcrun

%use the  corresponding seasonal data, needed to compute I0,R0
load seasonaldata;
season = season1;

%bootstrap function
%time = 1:length(season);
I0 = season(1);
I = @(Results) Results(:,2); % wrap the SIRFun to return only I
bootfun = @(params) I(SIRfun(time,params,[params(4) I0 1-params(4)-I0]));

nSamples = 10;


%%
ci = bootci(nSamples,{bootfun,chain},'alpha',0.9);

%%
errorbar(time,bootfun(mean(chain)),ci(1,:),ci(2,:))
hold on;
plot(time,season,'g');


%% plot the best fit (mean of the params)
plot(time,bootfun(mean(chain)),time,season);

%% plot the histogram of the stderr
 %cant do bc s2chain is -inf...
%% predictive plot
out = mcmcpred(results,chain,[],season,SIRfunforPred,100);