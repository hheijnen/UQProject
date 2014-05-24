function [I]=eval_sir(input_pars)

load('seasonaldata');

if nargin==0 %if no input given , use the standard params
alpha = 0.2*rand;
tau = 1.5;
m = 4*1e-6*rand;
S0 = 0.5;
Re = 1.1;
par=[alpha, tau, m, S0];
else
 par=input_pars;
end

% beta = Re*tau/(alpha*S0);
% options = odeset('MaxStep',0.01);
ode = @(t, y) sir(t, y, par); 

I0 = seasonaldata(1,1)/max(seasonaldata(:,1));
R0 = 1-S0-I0;
xinit=[S0 I0 R0];
for i=1:52
    tspan(i) = i;
end
[T,Y]=ode15s(ode,tspan,xinit);  

S = Y(:,1);
I = Y(:,2);
R = Y(:,3);
plot(I)
hold on
max_data = max(seasonaldata(:,1));
plot(seasonaldata(:,1)/max_data,'r')
% hold on
% plot(S,'r')
% plot(R,'g')

%change to if(0) or if(1) regarding if you want plots required
% if (0)
% figure(1)
% plot(T,-log(Y(:,3)));
% xlabel('t');
% ylabel('-ln(I)');
% 
% figure(2);
% plot(-log(Y(:,1)),-log(Y(:,3)));
% xlabel('-ln(S)');
% ylabel('-ln(I)');

end

