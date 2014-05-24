function dy = sir(t, y,theta)

alpha    = theta(1);
tau      = theta(2);
m        = theta(3);

S = y(1);
I = y(2);
R = y(3);
%tt = y(4);  % whats this?

%beta should depend on time: almost zero in summer, big enough in winter...
beta = -4.5*cos(pi/26)+5.5;
%% New
lambda = beta*(alpha*I + m);
dS = -lambda*S;
dI = lambda*S - tau*I;
dR = tau*I;


dy = [dS;dI;dR];