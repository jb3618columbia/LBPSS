function [ samples, log_likes, i, Zratio_ret] = AAG_RB_data(f, number_fn_evals, clique_size, initial_point, W, info_on_off)
%Added the Rao-Blackwellized estimate of pairwise emperical correlations.
%Here a is now a row vector of many (i,j) pairs: nodes for which pairwise error is computed.


d = f.dim;
L = round(number_fn_evals/(clique_size*(2*d-1)));
samples = zeros(d, L);
samples(:,1)=initial_point;   

% Sanity check: the error should start with 0 and remain zero 
% dist(:,1) = 0.5*ones(d,1);  
% This code passes that test

Zratio = zeros(1, L);
Zratio(1,1) = exp(0.5*initial_point'*W*initial_point); 

log_likes = zeros(1,L);
log_likes(1,1) = f.logp(initial_point);

for i=2:L
    xx = samples(:, i-1);  % Current point
    k=randperm(d);
    prob_vec = zeros(2*d,1);
    prob_vec(1,1) = f.logp(xx);
    if info_on_off == 1
        Zratio_est = zeros(1,2*d);
        Zratio_est(1,1) = exp(0.5*xx'*W*xx);
    end
    p = 2;
    j=1;
    
    for c=1:(2*d)-1 
        xx(k(j)) = -xx(k(j));
        % Efficient way to compute log likeiloohs of the propsoed point
        % cur_log_like = cur_log_like + sign(xx(k(j)))*f.logp_change(xx,k(j));
        prob_vec(p,1) = f.logp(xx);  % Inefficient
        
        if info_on_off == 1
            Zratio_est(1,p) = exp(0.5*xx'*W*xx);
        end
        p = p + 1;
        j = mod(c,d) + 1;
    end
    
    max_val = max(prob_vec);
    prob_vec = prob_vec - max_val;
    prob_vec = exp(prob_vec)/(sum(exp(prob_vec)));
    
    if info_on_off == 1
        Zratio(1,i) = Zratio_est*prob_vec; % adding the weighted average of pairwise count estimator
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
    log_likes(1,i) = f.logp(point);

end
Zratio_ret = log(mean(Zratio));
