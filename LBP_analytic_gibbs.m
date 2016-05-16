function [ samples, dist, log_likes, i ] = LBP_analytic_gibbs( f, number_fn_evals, clique_size, info_on_off, initial_point, marginals )
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

while fn_evals <= number_fn_evals
    xx = samples(:, i-1);  % Current point
    k=randperm(d);
    prob_vec = zeros(2*d,1);
    prob_vec(1,1) = f.logp(xx);
    dist_est = zeros(d,2*d);
    dist_est(:,1) = xx > 0;
    p = 2;
    j=1;
    z = 2*(rand(d,1) < marginals) - 1;
    
    for c=1:(2*d)-1 
        
        if xx(k(j)) ~= z(k(j))
            xx(k(j)) = -xx(k(j));
        end
        fn_evals = fn_evals + clique_size;
        % Efficient way to compute log likeiloohs of the propsoed point
        % cur_log_like = cur_log_like + sign(xx(k(j)))*f.logp_change(xx,k(j));
        prob_vec(p,1) = f.logp(xx) - (xx==1)'*log(marginals) - (xx==-1)'*log(1-marginals);  % Inefficient
        
        if info_on_off ==1
            dist_est(:,p) = xx > 0;
        end
        p = p + 1;
        j = mod(c,d) + 1;
        
    end
    
    max_val = max(prob_vec);
    prob_vec = prob_vec - max_val;
    prob_vec = exp(prob_vec)/(sum(exp(prob_vec)));
    
    if info_on_off == 1
        dist(:,i) = dist_est*prob_vec;   % getting the weighted average
    end
    
%     fn_evals = fn_evals + log(d);     % Extra log d for sampling from a dicrete distribution
    index = discretesample(prob_vec,1);
    
    % Find next point
    xx = samples(:, i-1);  % Current point
    if index ~= 1
        p = 2;
        j=1;
        for c=1:(2*d)-1 
            if xx(k(j)) ~= z(k(j))
                xx(k(j)) = -xx(k(j));
            end
            if p == index
                break;
            end
            p = p + 1;
            j = mod(c,d) + 1;

        end
    end
    samples(:, i) = xx;
    log_likes(i,1) = f.logp(samples(:, i));
    i = i + 1;

end

