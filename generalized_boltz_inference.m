% Runs experiments for general Boltzman machine

% Load dataset (available from http://artax.karlin.mff.cuni.cz/r-help/library/bnlearn/html/coronary.html)
load('coronary.mat')

% Parameters
n_iterations = 100;
n_samples_per_iteration = 10;
d_s = 5;
d_x = size(coronary,2);
W_sigma_2 = 1000;

%% Run samplers

% Initialize
s = zeros(1 , d_s);
W = zeros(d_s + d_x , d_s + d_x);

% Metropolis-Hastings sampler
% 1. Propose new W' from a symmetric proposal distribution p(W'|W)
% 2. Sample N new S' from p(S'|W') using MCMC methods
% 3. Compute "joint" distribution: J' = \frac{1}{N} * sum_{i=1}^N p(S',W')
% 4. Accept with probability p(W') J' / (p(W) + J). Here p(W) is the prior
% on W. Note that J should not be recomputed each iteration, only J' should
% be.
