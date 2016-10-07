function [ samples, log_likes, i, Zratio_ret] = AAS_RB_data(f, number_fn_evals, clique_size, initial_point, W)

% This sampler analytically finds all the coordinate flips that are above a
% certain threshold level and then uniformly samples one among them

%  There is no angle range here as we will search over the entire slice 

%  Done on a circle with 2d points being considered 
%  The current point and 2d-1 proposed points

%  Inputs:
%    1) Object f: having the log probability function
%    2) Number of fucntion evaluations to run this code for  
%    3) Clique size of the underlying model
%    4) Initial starting point
%    
% Outputs:
%    1) samples: D x L matrix of output samples
%    2) dist: D x L matrix of node marginals (computed by integrating over the slice)
%    3) log_likes: 1 x L vector of log-likelihoods of all samples
%    4) i: total number of samples

%Added the Rao-Blackwellized estimate of pairwise emperical correlations.
%Here a is now a row vector of many (i,j) pairs: nodes for which pairwise error is computed.

d = f.dim;
L = round(number_fn_evals/(clique_size*(2*d-1)));
samples = zeros(d, L);
samples(:,1)=initial_point;   

% Sanity check: the error should start with 0 and remain zero 
% dist(:,1) = 0.5*ones(d,1);  
% This code passes this test

Zratio = zeros(1, L);
Zratio(1,1) = exp(0.5*initial_point'*W*initial_point);

log_likes = zeros(1,L);
log_likes(1,1) = f.logp(initial_point);
cur_log_like = f.logp(initial_point);

for i=2:L
    xx = samples(:, i-1);  % Current point
    acc_samples = zeros(d,2*d);
    acc_samples(:,1) = xx;
    
    p = 2;
    hh = log(rand) + cur_log_like;
    k=randperm(d);
    j=1;
    
    for c=1:(2*d)-1
       xx(k(j)) = -xx(k(j));
       % Efficient way to compute log likeiloohs of the propsoed point 
%        cur_log_like = cur_log_like + sign(xx(k(j)))*f.logp_change(xx,k(j));
       cur_log_like = f.logp(xx);  % Inefficient 
       
       if cur_log_like > hh
           acc_samples(:,p) = xx;
           p = p + 1;
       end
       j = mod(c,d) + 1;
    end
    acc_samples( :, all(~acc_samples,1) ) = [];
    index = unidrnd(p-1);   
    samples(:,i) = acc_samples(:,index);

    Zratio(1,i) = sum(diag(exp(0.5*acc_samples'*W*acc_samples)))/p;  % RB estimate of th epartition function ratio
    cur_log_like = f.logp(samples(:,i));
    log_likes(1,i) = cur_log_like;
    
end
Zratio_ret = log(mean(Zratio));

