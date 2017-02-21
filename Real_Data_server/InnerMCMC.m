function [ logZratio ] = InnerMCMC( W_curr_mat, W_new_mat, d, F, alg)
% Inputs:
% 1) Current and proposed weight matrices
% 2) K: Number of samples in the inner MCMC loop

initial_point = sign(normrnd(0,1,d,1));
bias = zeros(d,1)';
bm_old = BM(d, W_curr_mat, bias);
bm = BM(d, W_new_mat, bias); % boltzmann machine object
clique_size = d-1;
W = W_new_mat - W_curr_mat;
% W is an upper traingular matrix

% Number of samples for other algorithms
t = 1.5;
T=t*pi;
fn_evlas_hmc = F*( (bm.dim)*t + (bm.dim) + clique_size*(bm.dim)*(t-0.5) );

if alg == 0
    logZratio = groundTruth(bm_old, bm); % passing the current and new obj.
end

if alg == 1
    [samples_hmc, ~, ~] = HMC_binary(bm, T, F, initial_point);
    logZratio = log(sum(diag(exp(samples_hmc'*W*samples_hmc)))/F);
end
if alg == 2
    [samples_CMH, ~, ~, ~] = CMH(bm, fn_evlas_hmc, clique_size, initial_point);
    logZratio = log(sum(diag(exp(samples_CMH'*W*samples_CMH)))/F);
end
if alg == 3
     [samples_AAS, ~, ~, logZratio] = AAS_RB_data(bm, fn_evlas_hmc, clique_size, initial_point, W);
end
if alg == 4
    info_on_off = 1;
    [samples_AAG, ~, ~, logZratio] = AAG_RB_data(bm, fn_evlas_hmc, clique_size, initial_point, W, info_on_off);
end
end

