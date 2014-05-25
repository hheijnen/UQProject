function [difference_sqr] = diff_sqr(par)

load('seasonaldata');

alpha = par(1);
tau = par(2);
m = par(3);
S0 = par(4);

par = [alpha, tau, m, S0];
ode = @(t, y) sir(t, y, par); 
I0 = mean(seasonaldata(1:3,1));
R0 = 1-S0-I0;
xinit=[S0 I0 R0];
for i=1:30
    tspan(i) = i;
end
[T,Y]=ode15s(ode,tspan,xinit);

I = Y(:,2);
difference_sqr = sum((I-seasonaldata(:,1)).^2);

plot(I)
hold on
plot(seasonaldata(:,1),'r')

end

