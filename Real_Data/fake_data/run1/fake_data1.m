% Runs experiments for general Boltzman machine
% Based on "Bayesian Learning in Undirected Graphical Models: Approximate
% MCMC algorithms" (BLUG)

% The goal is to sample W using Metropolis Hastings from its posterior p(W|S).
% Each iteration requires calculating a constant Z(W) which can be
% approximated using MCMC methods if it cannot be calculated exactly by brute force.
path = '/home/jalaj/Github_LBPSS/Real_Data/fake_data/run1/';

d=10;
% Set dummy parameter values
% W = triu(rand(d),1); % making sure the diagonals are zero
conn_strength = 1
W = triu(conn_strength*normrnd(0,1,d,d),1); % mixed connections
% Generating fake data
scale = 0.1
bias = scale*(2*rand(d,1)' -1);
bm_original = BM(d, W, bias); % original boltzmann machine object
coronary = genFakedata(bm_original);

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
truth = 1;
HMC = 1;
CMH = 1;
AAS = 1;
AAG = 1;
% Error in estimating the partition function ratio
Z_error_mat = zeros(4, num_sample_paths);
samples_true = zeros((num_samples*num_sample_paths), d*(d-1)/2);
samples_HMC = zeros((num_samples*num_sample_paths), d*(d-1)/2);
samples_CMH = zeros((num_samples*num_sample_paths), d*(d-1)/2);
samples_AAS = zeros((num_samples*num_sample_paths), d*(d-1)/2);
samples_AAG = zeros((num_samples*num_sample_paths), d*(d-1)/2);

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
    
    if HMC == 1
        tic
        display('HMC')
        [samples_HMC, logZ_est_HMC] = OuterMH( num_samples, F,  W_init, coronary, mu_prior, Sigma_prior, 1);
        Z_error_mat(1,k) = mean(abs(logZ_est_true - logZ_est_HMC));
        toc
    end
    
    if CMH ==1
        tic
        display('CMH')
        [samples_CMH, logZ_est_CMH] = OuterMH( num_samples, F,  W_init, coronary, mu_prior, Sigma_prior, 2);
        Z_error_mat(2,k) = mean(abs(logZ_est_true - logZ_est_CMH));
        toc
    end
    
    if AAS ==1
        tic
        display('AAS')
        [samples_AAS, logZ_est_AAS] = OuterMH( num_samples, F,  W_init, coronary, mu_prior, Sigma_prior, 3);
        Z_error_mat(3,k) = mean(abs(logZ_est_true - logZ_est_AAS));
        toc
    end
    
    if AAG ==1
        tic
        display('AAG')
        [samples_AAG, logZ_est_AAG] = OuterMH( num_samples, F,  W_init, coronary, mu_prior, Sigma_prior, 4);
        Z_error_mat(4,k) = mean(abs(logZ_est_true - logZ_est_AAG));
        toc
    end
    
end
samples_true = unique(samples_true, 'rows');
samples_HMC = unique(samples_HMC, 'rows');
samples_CMH = unique(samples_CMH, 'rows');
samples_AAS = unique(samples_AAS, 'rows');
samples_AAG = unique(samples_AAG, 'rows');

Z_error_mean = [mean(Z_error_mat,2), std(Z_error_mat,0,2)];
fileName = [path, 'Z_error', '.mat'];
save(fileName, 'Z_error_mean')

fileName = [path, 'tru', '.mat'];
save(fileName, 'samples_true')

samples_true1 = unique(samples_true, 'rows');
fileName = [path, 'tru1', '.mat'];
save(fileName, 'samples_true1')

fileName = [path, 'HMC', '.mat'];
save(fileName, 'samples_HMC')

samples_HMC1 = unique(samples_HMC, 'rows');
fileName = [path, 'HMC1', '.mat'];
save(fileName, 'samples_HMC1')


fileName = [path, 'CMH', '.mat'];
save(fileName, 'samples_CMH')

samples_CMH1 = unique(samples_CMH, 'rows');
fileName = [path, 'CMH1', '.mat'];
save(fileName, 'samples_CMH1')

fileName = [path, 'AAS', '.mat'];
save(fileName, 'samples_AAS')

samples_AAS1 = unique(samples_AAS, 'rows');
fileName = [path, 'AAS1', '.mat'];
save(fileName, 'samples_AAS1')


fileName = [path, 'AAG', '.mat'];
save(fileName, 'samples_AAG')

samples_AAG1 = unique(samples_AAG, 'rows');
fileName = [path, 'AAG1', '.mat'];
save(fileName, 'samples_AAG1')


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

