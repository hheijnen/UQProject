function [parF] = optimize_prior_SEIR(data,prior)

% Parametric uncertainty for SIR model

clear all;
close all;

load('seasonaldata.mat');



par0 = prior;

% options =  optimset('MaxIter',1000,'MaxFunEvals',2000);

% Calculation of the optimal parameters

[parF,fval] = fminsearch(@(par) diff_sqr_SEIR(par,data),par0);

% Calculation of the hessian
%[parF2, fval2, exitflag, output, grad, hessian] = fminunc(@(par) diff_sqr_SEIR(par), parF);

% N = length(seasonaldata(:,1));
% sigma = sqrt( 1/N * diff_sqr_SEIR(parF));
%H = 1/(2*sigma^2) * hessian;
%uncertainty = inv(H)


% % Hessian
% N = length(seasonaldata(:,1));
% sigma = sqrt( 1/(N-1) * diff_sqr(parF));
% [hess, err] = hessian(@(par) diff_sqr(par), parF);
% H = 1/(2*sigma^2) * hess


% % [parF2,fval2,exitflag,output,grad,hessian] = fminunc(@(par) diff_sqr(par), parF);
% % hessian



