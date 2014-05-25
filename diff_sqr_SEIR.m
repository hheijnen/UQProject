function [difference_sqr] = diff_sqr_SEIR(theta)

load('seasonaldata.mat');

alpha = theta(1);
tau = theta(2);
gamma = theta(3);
m = theta(4);
S0 = theta(5);
E0 = theta(6);

par = [alpha, tau, gamma, m, S0, E0];
ode = @(t, y) SEIR(t, y, par); 
I0 = season1(1,1);
R0 = 1-S0-E0-I0;
xinit=[S0 E0 I0 R0];

for i=1:length(season1)
    tspan(i) = i;
end
[T,Y]=ode15s(ode,tspan,xinit);

I = Y(:,3);
difference_sqr = sum((I-season1(:,1)).^2);

plot(I)
hold on
plot(seasonaldata(:,1),'r')

end

