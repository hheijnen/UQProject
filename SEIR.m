function dy = SEIR(t, y,theta)

alpha    = theta(1);
tau      = theta(2);
gamma    = theta(3);
m        = theta(4);
% thres     = theta(5);

% B0       = theta(5);
% B1       = theta(6);
% B2       = theta(7);

S = y(1);
E = y(2);
I = y(3);
R = y(4);

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
dE = lambda*S - gamma*E;
dI = gamma*E - tau*I;
dR = tau*I;


dy = [dS;dE;dI;dR];