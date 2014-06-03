clear
close all;

load('seasonaldata.mat');

data = season2;

amidone = 1;

while amidone
    alpha = rand*0.5;
    tau = rand*3 + 1;
    gamma = rand+0.5;
    m = 2.7*1e-5*rand;
    S0 = 0.2+0.5*rand;
    E0 = 0.01*rand;
    beta = 50;
    
    
    prior2(1) = alpha;
    prior2(2) = tau;
    prior2(3) = gamma;
    prior2(4) = m;
    prior2(5) = S0;
    prior2(6) = E0;
    
    try
        [parF,fval] = fminsearch(@(par) diff_sqr_SEIR(par,data),prior2);
    catch
        parF = [-1 0 0 0 0 0];
    end
    
    if (sum(parF > 0) == 6)
        save prior.mat parF prior2
        amidone = 0;
    end
end