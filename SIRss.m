function ss = SIRss(theta,ydata)

load debugdata.mat

time = [1:length(ydata)];

I0 = mean(ydata(1,1));

% 3 last parameters are the initial states
y0 = [theta(4) I0 1-theta(4)-I0]';

ymodel = SIRfun(time,theta,y0);
ss = sum(log(ymodel(:,2) - ydata)*2);

ssmat = [ssmat ss];
thetamat = [thetamat theta'];

hold on
plot(time,[ymodel(:,2) ydata])
ylim([0 max(ydata)])

save debugdata.mat ssmat thetamat

end

