function dy = sir(t, y,par)

%% Old, needs checking
mu    = par(1);
alpha = par(2);
gamma = par(3);
beta0 = par(4);
beta1 = par(5);

S = y(1);
E = y(2);
I = y(3);
tt = y(4);

%beta_SEIR depends on time! it is periodic and has to do with the contact
%rate!
beta_SEIR = @(beta0,beta1,t) beta0*(1+beta1*cos(2*pi*t));

%% New
lambda = beta_SEIR(beta0,beta1,tt)*(alpha*I + m);

dS = -lambda*S;
dI = lambda*S - tau*I;
dR = tau*I;


dy = [dS;dI;dR];