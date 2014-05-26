%% script which plots the optimal fits found by fminsearch
%  for uniformly distributed random initial guesses for the params

clc;
clear all;
close all;

load seasonaldata.mat

%% SIR optimization
%try fminsearch for different starting values and plot the results.

% optimal values
%alpha = 0.1523;
%tau = 1.7357;
%m = 2.7*1e-6;
%S0 = 0.2743;

%uniform bounds for all the params
%alpha [0,0.4]
%tau [0 3]
%m [0 1e-5]
%S0 [0 1]

%choose the season
season = season3;
tt= 1:length(season);
for i=1:4
    par0 = [0.4*rand(1,1) 3*rand(1,1) 1e-5*rand(1,1) rand(1,1)]'
    [parF,fval] = fminsearch(@(par) diff_sqr(par),par0);
    parF
    
    %solve the system with the optimal params
    y0 = [parF(4) season(1,1) 1-parF(4)-season(1,1)]'
    res = SIRfun(tt,parF,y0);
    
    %solve the system with the initial params
    y0init = [par0(4) season(1,1) 1-par0(4)-season(1,1)]';
    initres = SIRfun(tt,par0,y0init);
    
    
    % plotting of the results
    figure;
    alpha0 = par0(1);
    tau0   = par0(2);
    m0     = par0(3);
    S00    = par0(4);
    plot(tt,season,tt,res(:,2),...
        tt,initres(:,2));
    ylim([0 max(max(res(:,2)),max(season))]);
    xlabel('Time (Weeks)');
    ylabel('I');
    legend('data','optfit','initfit');
    
    txstr(1) = {'initial values'};
    txstr(2) = {['\alpha_0 ', num2str(alpha0)]};
    txstr(3) = {['\tau_0 ',num2str(tau0)]};
    txstr(4) = {['m_0 ',num2str(m0)]};
    txstr(5) = {['S0_0 ',num2str(S00)]};
    text(12,0.015,txstr,'HorizontalAlignment','left');
    
    
    txstr2(1) = {'optimal values'};
    txstr2(2) = {['\alpha ', num2str(parF(1))]};
    txstr2(3) = {['\tau ',num2str(parF(2))]};
    txstr2(4) = {['m ',num2str(parF(3))]};
    txstr2(5) = {['S0 ',num2str(parF(4))]};
    text(12,0.01,txstr2,'HorizontalAlignment','left');
    
end

%% and plot 'our' optimal values as wel
alpha = 0.1523;
tau = 1.7357;
m = 2.7*1e-6;
S0 = 0.2743;
parF = [alpha tau m S0]';
y0 = [parF(4) season(1,1) 1-parF(4)-season(1,1)]';
res = SIRfun(tt,parF,y0);

figure;
plot(tt,season,tt,res(:,2));
legend('data','optfit');

txstr2(1) = {'optimal values'};
txstr2(2) = {['\alpha ', num2str(parF(1))]};
txstr2(3) = {['\tau ',num2str(parF(2))]};
txstr2(4) = {['m ',num2str(parF(3))]};
txstr2(5) = {['S0 ',num2str(parF(4))]};
text(12,0.01,txstr2,'HorizontalAlignment','left');


%% SEIR optimization