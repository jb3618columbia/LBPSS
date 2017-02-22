function [ samples, dist, mag, log_likes, i, emp_count, emperical_counts ] = analytic_gibbs_new_ST( f, number_fn_evals, clique_size, info_on_off, initial_point, a )
% Added the Rao-Blackwellized estimate of pairwise emperical correlations.
% Here a is now a row vector of many (i,j) pairs: nodes for which pairwise error is computed.
% f is the Ising object


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
% emperical_counts(:,:,1) = empericalCounts(initial_point, a); % this is now a tensor of size K/2 x 4 x L

% Only on keeping this sum
emp_count = empericalCounts(initial_point, a);  

log_likes = zeros(1,L);
log_likes(1,1) = f.logp(initial_point);
mag = zeros(1,L);
mag(1,1) = sum(initial_point,1)/d;

% Gibbs sampling part changed to incorporate the Suwa-Todo algorithm.

for i=2:L
    xx = samples(:, i-1);  % Current point
    k=randperm(d);
    prob_vec = zeros(2*d,1);
    prob_vec(1,1) = f.logp(xx);
    if info_on_off == 1
        count_est = zeros(pair_size,4,2*d);
        count_est(:,:,1) = empericalCounts(xx, a);
        dist_est = zeros(d,2*d);
        dist_est(:,1) = xx > 0;
        mag_est = zeros(1,2*d);
        mag_est(:,1) = sum(xx,1)/d;
    end
    
    p = 2; j=1;
    
    for c=1:(2*d)-1 
        xx(k(j)) = -xx(k(j));
        % Efficient way to compute log likeiloohs of the propsoed point
        % cur_log_like = cur_log_like + sign(xx(k(j)))*f.logp_change(xx,k(j));
        prob_vec(p,1) = f.logp(xx);  % Inefficient
        if info_on_off == 1
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
    
    
    
    % Obtained the raw probability vector w (as in ST paper)
    % Sort such that w1 is the highest
    
    [m, pos] = max(prob_vec);
    
    % w and S vectors 
    place_holder = zeros(2*d,1);
    for z = 0:2*d -1
       place_holder(z+1,1) = prob_vec(mod(pos+z-1, 2*d)+1, 1);    
    end
        
    S = ones(1,2*d+1); %S_0=1, S_1=w_1   \
    S(2) = m;
    for b=3:2*d+1
        S(b) = S(b-1) + place_holder(b-1,1);      
    end
    
    % Finding the new transition probabilities
    prob_vec_new = zeros(2*d,1);
    if pos == 1 
        u=1;
        for q=0:2*d-1
            v = u+q;
            delta = S(u+1) - S(v) + m;
            prob_vec_new(q+1,1) = max(0, min([delta, place_holder(u,1) + place_holder(v,1) - delta, place_holder(u,1), place_holder(v,1)]) ) / (prob_vec(1,1));
        end
        
    else 
    u = 2*d - (pos-2);
        for q=0:2*d-1
           v = mod(u+q-1, 2*d) + 1;
           delta = S(u+1) - S(v) + m;
           prob_vec_new(q+1,1) = max(0, min([delta, place_holder(u,1) + place_holder(v,1) - delta, place_holder(u,1), place_holder(v,1)]) ) / (prob_vec(1,1));
           % getting the right transition probabilities 
        end 
    end
    
    % prob_vec_new are the right transition probabilities
    % we get the new estimators using prob_vec_new
    

    if info_on_off == 1
        dist(:,i) = dist_est*prob_vec_new;   % getting the weighted average
%         emperical_counts(:,:,i) = rb_emp_counts(count_est, prob_vec); % getting weighted average of pairwise count estimator
        emp_count = emp_count + rb_emp_counts(count_est, prob_vec_new); % adding the weighted average of pairwise count estimator
        mag(1,i) = mag_est*prob_vec_new;  % getting the weighted magnetization
    end
    
%     fn_evals = fn_evals + log(d);     % Extra log d for sampling from a dicrete distribution
    index = discretesample(prob_vec_new,1);
    
    if index == 1
       point = samples(:, i-1);
    elseif index > d + 1
       point = -samples(:, i-1);
       for j=1: (index - (d+1))
           point(k(j)) = -point(k(j));
       end
    else
       point = samples(:, i-1);
       for j=1:(index - 1)
           point(k(j)) = -point(k(j));
       end
    end
    samples(:, i) = point;
    log_likes(1,i) = f.logp(point);

end
% emp_count = mean(emperical_counts, 3);
emp_count = emp_count/L;

