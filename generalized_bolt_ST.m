% Runs experiments for general Boltzman machine
% Based on "Bayesian Learning in Undirected Graphical Models: Approximate
% MCMC algorithms" (BLUG)

% The goal is to sample W using Metropolis Hastings from its posterior p(W|S).
% Each iteration requires calculating a constant Z(W) which can be
% approximated using MCMC methods if it cannot be calculated exactly by brute force.
path = '/home/jalaj/Github_LBPSS/Real_Data/data/run5/';

% Load dataset (available from http://artax.karlin.mff.cuni.cz/r-help/library/bnlearn/html/coronary.html)
load('coronary.mat');

% Parameters
n = size(coronary,1);      % number of training examples
d = size(coronary,2);      % dimension of S
num_samples = 100         % number of samples for outre
F = 1000             % number of MCMC samples per inner loop
num_sample_paths = 1000

% Prior
mu_prior = zeros(1,d*(d-1)/2);
Sigma_prior = 3*eye(d*(d-1)/2);

% Algorithms
truth = 1;
AAG_ST = 1;
% Error in estimating the partition function ratio
Z_error_mat = zeros(1, num_sample_paths);
samples_true = zeros((num_samples*num_sample_paths), d*(d-1)/2);
samples_AAG_ST = zeros((num_samples*num_sample_paths), d*(d-1)/2);

for k=1:num_sample_paths
    
    display('Iteration Number =', num2str(k));
    % Initial point
    W_init = mvnrnd(mu_prior, sqrt(2)*eye(d*(d-1)/2));  % passed as row vectors
    
    %% Run samplers
    % initialW = iwishrnd(A,10)*(10-6-1);  % random PSD matrix, df are kept small
    % Wvec = nonzeros(triu(initialW, 1)');
    
    if true == 1
        tic
        display('True')
        [samples_true((k-1)*num_samples + 1:k*num_samples, :), logZ_est_true] = OuterMH( num_samples, F,  W_init, coronary, mu_prior, Sigma_prior, 0);
        toc
    end
    
    if AAG_ST ==1
        tic
        display('AAG_ST')
        [samples_AAG_ST((k-1)*num_samples + 1:k*num_samples, :), logZ_est_AAG_ST] = OuterMH( num_samples, F,  W_init, coronary, mu_prior, Sigma_prior, 4);
        Z_error_mat(1,k) = mean(abs(logZ_est_true - logZ_est_AAG_ST));
        toc
    end
    
end

Z_error_mean = [mean(Z_error_mat,2), std(Z_error_mat,0,2)];
fileName = [path, 'Z_error_mean', '.mat'];
save(fileName, 'Z_error_mean')

fileName = [path, 'tru', '.mat'];
save(fileName, 'samples_true')

samples_true1 = unique(samples_true, 'rows');
fileName = [path, 'tru1', '.mat'];
save(fileName, 'samples_true1')

fileName = [path, 'AAG_ST', '.mat'];
save(fileName, 'samples_AAG_ST')

samples_AAG_ST_1 = unique(samples_AAG_ST, 'rows');
fileName = [path, 'AAG_ST_1', '.mat'];
save(fileName, 'samples_AAG_ST1')

                                

% plotRealdata(samples_true, samples_HMC, samples_CMH, samples_AAS, samples_AAG);
% plotZerror(logZ_est_true, logZ_est_HMC, logZ_est_CMH, logZ_est_AAS, logZ_est_AAG)

%% Run samplers

% Initialize
% W = zeros(d,d);

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

