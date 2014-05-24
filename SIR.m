function dy = sir(t, y,par)

%% Old, needs checking
alpha    = par(1);
beta     = par(2);
tau      = par(3);
m        = par(4);
S0       = par(5);

S = y(1);
I = y(2);
R = y(3);
tt = y(4);  % whats this?

%beta_SEIR depends on time! it is periodic and has to do with the contact
%rate!
%beta_SEIR = @(beta0,beta1,t) beta0*(1+beta1*cos(2*pi*t));

%beta should depend on time: almost zero in summer, big enough in winter...
beta_SIR = @(beta,t) 1 + 10*t^2 + 28;
%% New
lambda = beta_SIR(beta,t)*(alpha*I + m);

dS = -lambda*S;
dI = lambda*S - tau*I;
dR = tau*I;


dy = [dS;dI;dR];