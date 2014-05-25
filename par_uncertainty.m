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

% Likelihood function

% diff_sqr = sum((I - seasonaldata(:,1)).^2);
% N = length(seasonaldata);
% sigma = sqrt(1/(N-1) * diff_sqr);
% Li = @(alpha, tau, m, beta) 1/(sqrt(2*pi)^N*sigma^N) * exp(-1/(2*sigma^2) * diff_sqr);

% (-1)*log posterior

% prior = 1; % set it equal to 1, doesnt play a role for the calculation of the Hessian
% L = -log(Li*prior);



