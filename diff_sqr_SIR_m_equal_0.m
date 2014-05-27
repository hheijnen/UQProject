function [difference_sqr] = diff_sqr_SIR_m_equal_0(theta,data)
% sum of squared differences of the SIR model with fixed tau

alpha = theta(1);
tau = 1.4;
m = 0;
S0 = theta(2);


par = [alpha, tau, m, S0];
par
ode = @(t, y) sir(t, y, par); 
I0 = data(1);
R0 = 1-S0-I0;
xinit=[S0 I0 R0];
tspan = 1:length(data);
[T,Y]=ode15s(ode,tspan,xinit);

I = Y(:,2);

difference_sqr = sum((I-data).^2);



end

