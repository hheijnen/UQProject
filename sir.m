function dy = sir(t, y,theta)

alpha    = theta(1);
tau      = theta(2);
m        = theta(3);
% thres     = theta(5);

% B0       = theta(5);
% B1       = theta(6);
% B2       = theta(7);

S = y(1);
I = y(2);
R = y(3);

%beta should depend on time: almost zero in summer, big enough in winter...
% beta = B2*normpdf(t,B0,B1);
% if (t > thres) 
    beta = 50;
% else
%     beta = 1;
% end

%% New
lambda = beta*(alpha*I + m);
dS = -lambda*S;
dI = lambda*S - tau*I;
dR = tau*I;


dy = [dS;dI;dR];