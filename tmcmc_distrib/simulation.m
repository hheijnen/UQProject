% Uncertainty Quantification 2014 spring - Project 2
% This script is devoted to simulating state variables (x1, x2,.. xN)
% given parameters (\theta)

function [ X ] = simulation( theta )

% initial values
x0 = [1,1,1,1,1,1,1,1]';

% time points
time_pt = [0:1:30];

% define kinetic rate functions
rates =@(c,x)[
           10^(c(1))*x(1)*x(7);
           10^(c(2))*x(3);
           10^(c(3))*x(2)*x(7);
           10^(c(4))*x(4);
           10^(c(5))*x(1);
           10^(c(6))*x(5);
           10^(c(7))*x(2);
           10^(c(8))*x(6);
           10^(c(9))*x(5);
           10^(c(10))*x(6);
           10^(c(11))*x(7);
           10^(c(12))*x(8)];
       
% define stoichiometry matrix
stoich = [-1 1 0 0 0 0 0 0 0 0 0 0;
          0 0 -1 1 0 0 0 0 0 0 0 0;
          1 -1 0 0 0 0 0 0 0 0 0 0;
          0 0 1 -1 0 0 0 0 0 0 0 0;
          0 0 0 0 1 0 0 0 -1 0 0 0;
          0 0 0 0 0 0 1 0 0 -1 0 0;
          -1 1 -1 1 0 1 0 0 0 0 -1 0;
          0 0 0 0 0 0 0 1 0 0 0 -1];

% define the ODE
dydt = @(c,x) stoich*rates(c,x);
odefun =@(t,x) dydt(theta,x);

% integrator options
options = odeset('Stats','on','AbsTol',1e-9);

% perform integration
[T,X] = ode15s(odefun,time_pt,x0,options);

end