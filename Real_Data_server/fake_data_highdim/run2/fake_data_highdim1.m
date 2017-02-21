% Runs experiments for general Boltzman machine
% Based on "Bayesian Learning in Undirected Graphical Models: Approximate
% MCMC algorithms" (BLUG)

% The goal is to sample W using Metropolis Hastings from its posterior p(W|S).
% Each iteration requires calculating a constant Z(W) which can be
% approximated using MCMC methods if it cannot be calculated exactly by brute force.
path = '/home/jalaj/Github_LBPSS/Real_Data/fake_data_highdim/run2/';

% Generating fake data using HMC on Boltzmann machine.
% Set dummy parameter values
d=25;
% W = triu(rand(d),1); % making sure the diagonals are zero
conn_strength = 1
W = triu(conn_strength*normrnd(0,1,d,d),1);
scale = 0.1;
bias = scale*(2*rand(d,1)' -1);
bm_original = BM(d, W, bias); % original boltzmann machine object

initial_point = sign(normrnd(0,1,d,1));
t = 1.5;
T=t*pi;
num_data_pts=200;
[samples_hmc, ~, ~] = HMC_binary(bm_original, T, num_data_pts, initial_point);
coronary = samples_hmc'; % fake data generated

% Parameters
n = size(coronary,1);      % number of training examples
d = size(coronary,2);      % dimension of S
num_samples = 100          % number of samples
F = 500                   % number of MCMC samples per inner loop
num_sample_paths = 250

% Prior
% mu_prior = zeros(1,d*(d-1)/2);
mu_prior = nonzeros(W);
Sigma_prior = 2*eye(d*(d-1)/2);

% Algorithms
HMC = 1;
CMH = 1;
AAS = 1;
AAG = 1;
% Error in estimating the partition function ratio
samples_HMC = zeros((num_samples*num_sample_paths), d*(d-1)/2);
samples_CMH = zeros((num_samples*num_sample_paths), d*(d-1)/2);
samples_AAS = zeros((num_samples*num_sample_paths), d*(d-1)/2);
samples_AAG = zeros((num_samples*num_sample_paths), d*(d-1)/2);

for k=1:num_sample_paths
    
    display('Iteration Number =', num2str(k));
    % Initial point
    W_init = mvnrnd(mu_prior, 0.5*eye(d*(d-1)/2));  % passed as row vectors
    
    %% Run samplers
    % initialW = iwishrnd(A,10)*(10-6-1);  % random PSD matrix, df are kept small
    % Wvec = nonzeros(triu(initialW, 1)');
    
    if HMC == 1
        tic
        display('HMC')
        [samples_HMC((k-1)*num_samples + 1:k*num_samples, :), logZ_est_HMC] = OuterMH( num_samples, F,  W_init, coronary, mu_prior, Sigma_prior, 1);
        toc
    end
    
    if CMH ==1
        tic
        display('CMH')
        [samples_CMH((k-1)*num_samples + 1:k*num_samples, :), logZ_est_CMH] = OuterMH( num_samples, F,  W_init, coronary, mu_prior, Sigma_prior, 2);
        toc
    end
    
    if AAS ==1
        tic
        display('AAS')
        [samples_AAS((k-1)*num_samples + 1:k*num_samples, :), logZ_est_AAS] = OuterMH( num_samples, F,  W_init, coronary, mu_prior, Sigma_prior, 3);
        toc
    end
    
    if AAG ==1
        tic
        display('AAG')
        [samples_AAG((k-1)*num_samples + 1:k*num_samples, :), logZ_est_AAG] = OuterMH( num_samples, F,  W_init, coronary, mu_prior, Sigma_prior, 4);
        toc
    end
    
end
samples_HMC = unique(samples_HMC, 'rows');
samples_CMH = unique(samples_CMH, 'rows');
samples_AAS = unique(samples_AAS, 'rows');
samples_AAG = unique(samples_AAG, 'rows');


fileName = [path, 'HMC', '.mat'];
save(fileName, 'samples_HMC')

fileName = [path, 'CMH', '.mat'];
save(fileName, 'samples_CMH')

fileName = [path, 'AAS', '.mat'];
save(fileName, 'samples_AAS')

fileName = [path, 'AAG', '.mat'];
save(fileName, 'samples_AAG')

fileName = [path, 'W', '.mat'];
save(fileName, 'W')

fileName = [path, 'bias', '.mat'];
save(fileName, 'bias')



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

