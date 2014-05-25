% Parametric uncertainty calculation

clear all;
close all;

load('seasonaldata.mat');

% optimal value
alpha = 0.1523;
tau = 1.7357;
gamma = 0.5;
m = 2.7*1e-6;
S0 = 0.2743;
E0 = 0.05;
beta = 50;

theta(1) = alpha;
theta(2) = tau;
theta(3) = gamma;
theta(4) = m;
theta(5) = S0;
theta(6) = E0;
par0 = theta;

[parF,fval] = fminsearch(@(par) diff_sqr_SEIR(par),par0);

theta_opt = parF

options = optimset('TolFun', 1e-6);
options = optimset('TolX', 1e-6);
[parF2,fval2,exitflag,output,grad,hessian] = fminunc(@(par) diff_sqr(par), parF);
% [hess, err] = hessian(@(par) diff_sqr(par), parF);

hessian



