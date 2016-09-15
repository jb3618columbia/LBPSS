function [ est_dist ] = emp_dist( samples)
% To estimate the probability distribution implied by the samples obtained 

% Will run for:
% 1) Ground truth 
% 2) All other algortihms 

% Input: d x L matrix of samples
% Output: vector d X 1 of prob(s_i=1)


[m,n] = size(samples);
est_dist = zeros(m,1);

for i=1:1:m
    
    est_dist(i,1) = nnz(samples(i,:) == 1)/n;
end


end

