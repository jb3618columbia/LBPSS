function [ sample_prior] = samplePrior(d)

% To generate a sample from the prior. 
% 1) Uniform Prior 
% 2) LBP Prior

% Currently doing from a uniform prior only, to add LBP prior 

% For Ising Model we have states = [-1, +1]

sample_prior = unidrnd(2, d, 1) - 2;
y = sample_prior >= 0;
sample_prior(y) = 1;

