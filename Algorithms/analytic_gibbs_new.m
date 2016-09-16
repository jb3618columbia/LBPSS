function [ samples, dist, mag, log_likes, i, emp_count, emperical_counts ] = analytic_gibbs_new( f, number_fn_evals, clique_size, info_on_off, initial_point, a )
%Added the Rao-Blackwellized estimate of pairwise emperical correlations.
%Here a is now a row vector of many (i,j) pairs: nodes for which pairwise error is computed.


d = f.dim;
L = round(number_fn_evals/(clique_size*(2*d-1)));
samples = zeros(d, L);
samples(:,1)=initial_point;   

% Sanity check: the error should start with 0 and remain zero 
% dist(:,1) = 0.5*ones(d,1);  
% This code passes that test

dist(:,1) = emp_dist(initial_point);  
pair_size = size(a,2)/2;
emperical_counts = zeros(pair_size, 4, L);
emperical_counts(:,:,1) = empericalCounts(initial_point, a); % this is now a tensor of size K/2 x 4 x L

log_likes = zeros(1,L);
log_likes(1,1) = f.logp(initial_point);
mag = zeros(1,L);
mag(1,1) = sum(initial_point,1)/d;

for i=2:L
    xx = samples(:, i-1);  % Current point
    k=randperm(d);
    prob_vec = zeros(2*d,1);
    prob_vec(1,1) = f.logp(xx);
    count_est = zeros(pair_size,4,2*d);
    count_est(:,:,1) = empericalCounts(xx, a);
    dist_est = zeros(d,2*d);
    dist_est(:,1) = xx > 0;
    mag_est = zeros(1,2*d);
    mag_est(:,1) = sum(xx,1)/d;
    p = 2;
    j=1;
    
    for c=1:(2*d)-1 
        
        xx(k(j)) = -xx(k(j));
        % Efficient way to compute log likeiloohs of the propsoed point
        % cur_log_like = cur_log_like + sign(xx(k(j)))*f.logp_change(xx,k(j));
        prob_vec(p,1) = f.logp(xx);  % Inefficient
        
        if info_on_off ==1
            dist_est(:,p) = xx > 0;
            count_est(:,:,p) = empericalCounts(xx, a);
            mag_est(:,p) = sum(xx,1)/d;
        end
        p = p + 1;
        j = mod(c,d) + 1;
        
    end
    
    max_val = max(prob_vec);
    prob_vec = prob_vec - max_val;
    prob_vec = exp(prob_vec)/(sum(exp(prob_vec)));
    
    if info_on_off == 1
        dist(:,i) = dist_est*prob_vec;   % getting the weighted average
        emperical_counts(:,:,i) = rb_emp_counts(count_est, prob_vec); % getting weighted average of pairwise count estimator
        mag(1,i) = mag_est*prob_vec;  % getting the weighted magnetization
    end
    
%     fn_evals = fn_evals + log(d);     % Extra log d for sampling from a dicrete distribution
    index = discretesample(prob_vec,1);
    
    if index == 1
       point = samples(:, i-1);
    elseif index > d + 1
       point = -samples(:, i-1);
       for j=1:index - (d+1)
           point(k(j)) = -point(k(j));
       end
    else
       point = samples(:, i-1);
       for j=1:index - 1
           point(k(j)) = -point(k(j));
       end
    end
    samples(:, i) = point;
    log_likes(i,1) = f.logp(point);

end
emp_count = mean(emperical_counts, 3);
