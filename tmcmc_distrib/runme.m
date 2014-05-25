% Main program that integrates model parameterization and 
% Bayesian parameter estimation based on stochastic simulation algorithms
% for model updating 
% Transitional Markov Chain Monte Carlo (Ching, J. and Chen, Y. (2007).
% ?Transitional Markov Chain Monte Carlo Method for Bayesian Model Updating, 
% Model Class Selection, and Model Averaging.? J. Eng. Mech., 133(7), 816?832.)
%
% If using this code in a publication pls cite J Chem Phys. 2012 Oct 14;137(14):144103. doi: 10.1063/1.4757266.
%
% Matlab Code distribution For the purpose of UQ_14 class
% shoot as is to sample from a mutlimodal function.
%Go and edit the userinput file and uq14_model to adapt it to your
%likelihoods and priors!
close all;  clear all;  clc;

%Multicore Version-Brutus version with scripting on request
try
    matlabpool close force local;   %in case it was left open
catch
end

%User parameters to be changed
userinput;



%--------------------------------------------------------------------------
% TMCMC run
%--------------------------------------------------------------------------
tic;
[TMCMC]=master_distrib(data);
toc;
%
%
%EXITING
fprintf('Finished with Run, all output are in the field TMCMC\n');
fprintf('The Model Evidence is also on TMCMC.logEvidence!\n');

