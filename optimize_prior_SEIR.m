function [parF] = optimize_prior_SEIR(data,theta)

% Parametric uncertainty for SIR model

load('seasonaldata.mat');



par0 = theta;

% figure
[parF,fval] = fminsearch(@(par) diff_sqr_SEIR(par,data),par0);
% hold off




