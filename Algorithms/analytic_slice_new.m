function [ samples, dist, mag, log_likes, i, emp_count, emperical_counts ] = analytic_slice_new( f, number_fn_evals, clique_size, initial_point, a)

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

dist(:,1) = emp_dist(initial_point);    
pair_size = size(a,2)/2;
emperical_counts = zeros(pair_size, 4, L);
emperical_counts(:,:,1) = empericalCounts(initial_point, a); % this is now a tensor of size K/2 x 4 x L

log_likes = zeros(1,L);
log_likes(1,1) = f.logp(initial_point);
cur_log_like = f.logp(initial_point);
mag = zeros(1,L);
mag(1,1) = sum(initial_point,1)/d;

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
    dist(:,i) = emp_dist(acc_samples(:,1:end));

    emperical_counts(:,:,i) = empericalCounts(acc_samples, a); % this is now a tensor of size K/2 x 4 x L
    mag(1,i) = sum(sample_to_mag(acc_samples));
    cur_log_like = f.logp(samples(:,i));
    log_likes(1,i) = cur_log_like;
    
end
emp_count = mean(emperical_counts, 3);