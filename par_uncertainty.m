% Parametric uncertainty calculation

clear all;
close all;

load('seasonaldata');

% optimal value
alpha = 0.1523;
tau = 1.7357;
m = 2.7*1e-6;
S0 = 0.2743;
beta = 50;

theta(1) = alpha;
theta(2) = tau;
theta(3) = m;
theta(4) = S0;
par0 = theta;

[parF,fval] = fminsearch(@(par) diff_sqr(par),par0);

alpha_opt = parF(1)
tau_opt = parF(2)
m_opt = parF(3)
S0_opt = parF(4)

options = optimset('TolFun', 1e-6);
options = optimset('TolX', 1e-6);
[parF2,fval2,exitflag,output,grad,hessian] = fminunc(@(par) diff_sqr(par), parF);
% [hess, err] = hessian(@(par) diff_sqr(par), parF);

hessian



