% seir.m
%
% Imlements a SEIR infection model
%
%
%http://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology#The_SEIR_model
%
%Also : Seasonality and period-doubling bifurcations in an epidemic model.
%J Theor Biol. 1984 Oct 21;110(4):665-79. Aron JL, Schwartz IB.
% 
%
%S is the fraction of susceptible individuals (those able to contract the disease),
%E is the fraction of exposed individuals (those who have been infected but are not yet infectious),
%I is the fraction of infective individuals (those capable of transmitting the disease),
%R is the fraction of recovered individuals (those who have become immune).
%Note that the variables give the fraction of individuals - that is, they are 
%S+E+I+R = 1.
%
% Since one of the params (check cited paper) is time varying , we have a non-autonomous
% system.
% However, we can convert it to a four-dimensional autonomous system
% Inputs:
%   t - Time variable:
%   y - Independent variable: populations fractions (1-3), 4 is time
% par - parameters regarding infection and migration rates
% Output:
%   dy - First derivative: the rate of change 
% 
%UQ class 2014 Spring Semester ETH-Z

function dy = seir(t, y,par)
%% Old, needs checking

mu  = par(1);
alpha=par(2);
gamma=par(3);
beta0=par(4);
beta1=par(5);

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
dE = lambda*S - sigma*E;
dI = sigma*E - tau*I;
dR = tau*I;

dy = [dS;dE;dI;dR];