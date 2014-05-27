function [ parF ] = optimize_prior_SIR_m_equal_0(data,theta)
%TOptimizes the params for SIR using fminsearch
%   given the start values in theta (tau gets fixed to 1.4, m to 0)

par0 = theta;

%leave out tau in the optimization
par0 = [par0(1) par0(4)];
[parF,fval] = fminsearch(@(par) diff_sqr_SIR_m_equal_0(par,data),par0);

%put tau back in
parF = [parF(1) 1.4 0 parF(2)];
end

