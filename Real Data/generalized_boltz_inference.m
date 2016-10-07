% Runs experiments for general Boltzman machine
% Based on "Bayesian Learning in Undirected Graphical Models: Approximate
% MCMC algorithms" (BLUG)

% The goal is to sample W using Metropolis Hastings from its posterior p(W|S).
% Each iteration requires calculating a constant Z(W) which can be
% approximated using MCMC methods if it cannot be calculated exactly by brute force.

% Load dataset (available from http://artax.karlin.mff.cuni.cz/r-help/library/bnlearn/html/coronary.html)
load('coronary.mat');

% Parameters
n = size(coronary,1);      % number of training examples
d = size(coronary,2);      % dimension of S
num_samples = 100;          % number of samples
F = 1000;                   % number of MCMC samples per inner loop

% Prior
mu_prior = zeros(1,d*(d-1)/2);
Sigma_prior = 10*eye(d*(d-1)/2);

% Initial point
W_init = mvnrnd(mu_prior, eye(d*(d-1)/2));  % passed as row vectors 

% Algorithms 
truth = 1;
HMC = 1;
CMH = 1;
AAS = 1;
AAG = 1;

%% Run samplers
% initialW = iwishrnd(A,10)*(10-6-1);  % random PSD matrix, df are kept small 
% Wvec = nonzeros(triu(initialW, 1)');

if true == 1
    tic
    display('True')
    [samples_true, logZ_est_true] = OuterMH( num_samples, F,  W_init, coronary, mu_prior, Sigma_prior, 0);
    samples_true = unique(samples_true, 'rows');
    toc
end

if HMC == 1
    tic
    display('HMC')
    [samples_HMC, logZ_est_HMC] = OuterMH( num_samples, F,  W_init, coronary, mu_prior, Sigma_prior, 1);
    samples_HMC = unique(samples_HMC, 'rows');
    toc
end

if CMH ==1
    tic
    display('CMH')
    [samples_CMH, logZ_est_CMH] = OuterMH( num_samples, F,  W_init, coronary, mu_prior, Sigma_prior, 2);
    samples_CMH = unique(samples_CMH, 'rows');
    toc
end

if AAS ==1
    tic
    display('AAS')
    [samples_AAS, logZ_est_AAS] = OuterMH( num_samples, F,  W_init, coronary, mu_prior, Sigma_prior, 3);
    samples_AAS = unique(samples_AAS, 'rows');
    toc
end

if AAG ==1
    tic
    display('AAG')
    [samples_AAG, logZ_est_AAG] = OuterMH( num_samples, F,  W_init, coronary, mu_prior, Sigma_prior, 4);
    samples_AAG = unique(samples_AAG, 'rows');
    toc
end

% plotRealdata(samples_true, samples_HMC, samples_CMH, samples_AAS, samples_AAG);
plotZerror(logZ_est_true, logZ_est_HMC, logZ_est_CMH, logZ_est_AAS, logZ_est_AAG)

%% Run samplers

% Initialize
W = zeros(d,d);

% Metropolis-Hastings sampler

% 1. Propose new W' from a symmetric proposal distribution p(W'|W)

% 2. Calculate F = Z(W)/Z(W') ~ \frac{1}{K} * \sum_{k=1}^K exp(s_k(W-W')*s_k) 
%    where we have sampled K values of S from p(S|W') using an MCMC method.
%    This is equation (8) of BLUG.
%    Note for the exact case instead of sampling we can calculate
%    F = Z(W)/Z(W') directly using brute force. 

% 3. Accept W' with probability p(W') / p(W) * F^N * exp(\sum_{i=1}^n s_i(W-W')*s_i)
%    where p(W) is the prior on W and s_i are the datapoints with i=1,...,n.
%    This is equation (7) of BLUG.

%% Plots

% Plot empirical distributions of W for the exact sampler vs approximate
% samplers. It is probably only worth comparing LBP-CMH and LBP-AAG with
% the exact sampler.

