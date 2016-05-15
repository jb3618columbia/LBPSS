function [ samples, dist, log_likes, i, tv ] = analytic_gibbs_new( f, number_fn_evals, clique_size, info_on_off, initial_point, true_dist )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here


d = f.dim;
samples(:,1)=initial_point;   

% Sanity check: the error should start with 0 and remain zero 
% dist(:,1) = 0.5*ones(d,1);  
% This code passes that test

dist(:,1) = emp_dist(initial_point);                           
log_likes(1,1) = f.logp(initial_point);
fn_evals = 0;
i=2;

est_dist = zeros(1,2^d);    % for computing the total variation

while fn_evals <= number_fn_evals
    xx = samples(:, i-1);  % Current point
    k=randperm(d);
    prob_vec = zeros(2*d,1);
    prob_vec(1,1) = f.logp(xx);
    dist_est = zeros(d,2*d);
    dist_est(:,1) = xx > 0;
    
    acc_samples = zeros(d,2*d);
    acc_samples(:,1) = xx;
    p = 2;
    j=1;
    
    for c=1:(2*d)-1 
        
        xx(k(j)) = -xx(k(j));
        fn_evals = fn_evals + clique_size;
        % Efficient way to compute log likeiloohs of the propsoed point
        % cur_log_like = cur_log_like + sign(xx(k(j)))*f.logp_change(xx,k(j));
        prob_vec(p,1) = f.logp(xx);  % Inefficient
        
        if info_on_off ==1
            dist_est(:,p) = xx > 0;
            acc_samples(:,p) = xx;
        end
        p = p + 1;
        j = mod(c,d) + 1;
        
    end
    
     
    log_prob = prob_vec;
    max_val = max(prob_vec);
    prob_vec = prob_vec - max_val;
    prob_vec = exp(prob_vec)/(sum(exp(prob_vec)));
    
%     for m=1:p-1
%         yy = (acc_samples(:,m) +1)/2;
%         yy_1 = bi2de(yy');
%         est_dist(1,yy_1+1) = log_prob(m,1)*prob_vec(m,1);
%      end
    
    if info_on_off == 1
        dist(:,i) = dist_est*prob_vec;   % getting the weighted average
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
    log_likes(i,1) = f.logp(xx);
    i = i + 1;

end

est_dist = exp(est_dist)/sum(exp(est_dist));
tv = sum(abs(true_dist - est_dist));
end

