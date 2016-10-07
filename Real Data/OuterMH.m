function [ samples_W,  logZ_est] = OuterMH( num_samples, F, initialpt, coronary, mu_prior, Sigma_prior, alg )
% Performs the outer MH loop and returns a 3d array of samples
% Will assume the diagonal elements of a Boltzmann machine are 0

% Inputs: initialW vector of upper triangular elemnts of W
n = size(coronary,1);
d = size(coronary,2);

logZ_est = zeros(1, num_samples-1);
samples_W = zeros(num_samples, size(initialpt, 2));  % W is a d x d matrix
samples_W(1,:) = initialpt;
Sigma = .01*eye(size(initialpt, 2));

for i=2:num_samples
    W_curr = samples_W(i-1, :);
    W_new = mvnrnd(W_curr, Sigma);  % symmetric proposal distribution
    
    W_curr_mat = triu(ones(d),1);
    W_curr_mat(~~W_curr_mat)=W_curr';
    
    W_new_mat = triu(ones(d),1);
    W_new_mat(~~W_new_mat) = W_new';
    
    logpriorratio = logmvnpdf(W_new,mu_prior,Sigma_prior) - logmvnpdf(W_curr,mu_prior,Sigma_prior);
    loglik = loglikfn(W_curr_mat, W_new_mat, coronary);
    logZratio = InnerMCMC(W_curr_mat, W_new_mat, d, F, alg); % get this from inner MCMC loop
    logZ_est(1, i-1) = logZratio;
    a = logpriorratio + n*logZratio + loglik;
    a = exp(min(1,a));
    
    if rand < a % accept
        samples_W(i,:) = W_new;
    else % reject
        samples_W(i,:) = W_curr;
    end
    
    
end

end

