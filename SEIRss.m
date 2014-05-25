function ss = SEIRss(theta,ydata)

load debugdata.mat

time = [1:length(ydata)];

I0 = mean(ydata(1,1));

% 3 last parameters are the initial states
y0 = [theta(5) theta(6) I0 1-theta(5)-2*I0]';

ymodel = SEIRfun(time,theta,y0);
ss = sum((ymodel(:,3) - ydata).^4);

ssmat = [ssmat ss];
thetamat = [thetamat theta'];

hold on
plot(time,[ymodel(:,3) ydata])
ylim([0 max(ydata)])

save debugdata.mat ssmat thetamat

end

