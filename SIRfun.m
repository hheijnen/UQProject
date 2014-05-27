function y = SIRfun(time, theta, y0)

[t,y] = ode15s(@sir,time,y0,[],theta);


end

