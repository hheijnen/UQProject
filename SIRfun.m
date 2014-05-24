function y = SIRfun(time, theta, y0)

[t,y] = ode15s(@SIR,time,y0,[],theta);


end

