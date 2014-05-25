function y = SEIRfun(time, theta, y0)

[t,y] = ode15s(@SEIR,time,y0,[],theta);


end