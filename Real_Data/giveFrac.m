function [ frac ] = giveFrac(samples, W, percent)

% samples: as those taken from MCMC
% W: the true value used for generating data
ground_truth = nonzeros(W)';
n=15000;
d=size(samples,2);
frac = zeros(1,d);

for i=1:d
    if ground_truth(1,i) > 0
        ll=ground_truth(1,i)*(1-percent);
        ul=ground_truth(1,i)*(1+percent);
    else
        ul=ground_truth(1,i)*(1-percent);
        ll=ground_truth(1,i)*(1+percent);
    end
    frac(1,i) = sum( (samples(:,i) >= ll).*(samples(:,i) <= ul) )/n;   
end
frac= sort(frac, 'descend');