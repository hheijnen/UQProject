function [difference_sqr] = diff_sqr_SEIR(theta,data)

load('seasonaldata.mat');
season = data;


alpha = theta(1);
tau = theta(2);
gamma = theta(3);
m = theta(4);
S0 = theta(5);
E0 = theta(6);

par = [alpha, tau, gamma, m, S0, E0];
ode = @(t, y) SEIR(t, y, par); 
I0 = season(1,1);
R0 = 1-S0-E0-I0;
xinit=[S0 E0 I0 R0];

for i=1:length(season)
    tspan(i) = i;
end
[T,Y]=ode15s(ode,tspan,xinit);

I = Y(:,3);
difference_sqr = sum((I-season).^2);

plot(I)
hold on
plot(season,'r')

end

