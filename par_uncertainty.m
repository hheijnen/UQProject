% Parametric uncertainty for SIR model

clear all;
close all;

load('seasonaldata.mat');

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

% options =  optimset('MaxIter',1000,'MaxFunEvals',2000);

% Calculation of the optimal parameters

[parF,fval] = fminsearch(@(par) diff_sqr(par),par0);
parF

% Calculation of the hessian
[parF2, fval2, exitflag, output, grad, hessian] = fminunc(@(par) diff_sqr(par), parF);

N = length(seasonaldata(:,1));
sigma = sqrt( 1/N * diff_sqr(parF));
H = 1/(2*sigma^2) * hessian;
uncertainty = inv(H)


% % Hessian
% N = length(seasonaldata(:,1));
% sigma = sqrt( 1/(N-1) * diff_sqr(parF));
% [hess, err] = hessian(@(par) diff_sqr(par), parF);
% H = 1/(2*sigma^2) * hess


% % [parF2,fval2,exitflag,output,grad,hessian] = fminunc(@(par) diff_sqr(par), parF);
% % hessian



