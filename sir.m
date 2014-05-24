function dy = sir(t, y, par)

%% Old, needs checking
alpha    = par(1);
tau      = par(2);
m        = par(3);
S0       = par(4);


beta0 = 1e6;
beta1 = 1e7;

S = y(1);
I = y(2);
R = y(3);
% tt = y(4);  % whats this?

%beta_SEIR depends on time! it is periodic and has to do with the contact
%rate!
beta_SIR = @(beta0,beta1,t) beta0*(1+beta1*cos(2*pi*t));

%beta should depend on time: almost zero in summer, big enough in winter...
% if (t <= 8 || t > 48)
% if (t <= 8 || t > 48)
%     beta_SIR = 1;
% else
%     beta_SIR = 1+beta;
% end

% beta_SIR = @(t,beta0) beta0*cos(2*pi*t);

%% New
lambda = beta_SIR(t,beta0,beta1)*(alpha*I + m);

dS = -lambda*S;
dI = lambda*S - tau*I;
dR = tau*I;


dy = [dS;dI;dR];