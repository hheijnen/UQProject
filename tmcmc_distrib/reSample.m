function [ IX ] = reSample( w )
%RESAMPLE Residual Resampling (Liu et al)
%   IX = reSample(w), where w is weights and IX is index set of the
%   resulting particles
n = length(w);
w = n*w/sum(w); % normalize to sum up to n
wN = floor(w); % integer parts
wR = w-wN; % residual weigths
wR = wR/sum(wR); % normalize

% filling indexes with integer parts
k = 1;
IX = zeros(1,n);
for i=1:n
    for j = 1:wN(i)
        IX(k) = i;
        k = k+1;
    end
end

% use residuals to fill rest with roulette wheel selection
cs = cumsum(wR);
for j = k:n
    ix = find(cs > rand(),1);
    IX(j) = ix;
end

end