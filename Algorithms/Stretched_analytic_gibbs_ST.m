function [ samples, dist, mag, log_likes, i, emp_count_gibbs, emperical_counts] = Stretched_analytic_gibbs_ST( f, number_fn_evals, clique_size, info_on_off, initial_point, marginals, a )
%Added the Rao-Blackwellized estimate of pairwise emperical correlations.
%Here a is a vector of (i,j): nodes for which pairwise error is computed.


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
emperical_counts(:,:,1) = empericalCounts(initial_point, a); % this is now a tensor of size K/2 x 4 x N

log_likes = zeros(1,L);
log_likes(1,1) = f.logp(initial_point);
mag = zeros(1,L);
mag(1,1) = sum(initial_point,1)/d;

for i=2:L
    xx = samples(:, i-1);  % Current point
    t_bracket_width = ((xx==1).*marginals+(xx==-1).*(1-marginals))*pi;
    t = (2*rand(d,1)-1).*t_bracket_width;
    t_brackets = mod([t-t_bracket_width,t+t_bracket_width],2*pi);
    [times, flip_indices_sorted] = sort([t_brackets(1:d), t_brackets(d+1:2*d)]);
    k = mod(flip_indices_sorted-1,d)+1;
    
    prob_vec = zeros(2*d,1);
    prob_vec(1,1) = log(times(1)+2*pi-times(2*d))+f.logp(xx) - (xx==1)'*log(marginals) - (xx==-1)'*log(1-marginals);
    count_est = zeros(pair_size,4,2*d);
    count_est(:,:,1) = empericalCounts(xx, a);
    dist_est = zeros(d,2*d);
    dist_est(:,1) = (xx > 0);
    mag_est = zeros(1,2*d);
    mag_est(:,1) = sum(xx,1)/d;
    
    for c=1:(2*d)-1 
        
        xx(k(c)) = -xx(k(c));
        % Efficient way to compute log likeiloohs of the propsoed point
        % cur_log_like = cur_log_like + sign(xx(k(j)))*f.logp_change(xx,k(j));
        prob_vec(c+1,1) = log(times(c+1)-times(c)) + f.logp(xx) - (xx==1)'*log(marginals) - (xx==-1)'*log(1-marginals);  % Inefficient
        
        if info_on_off ==1
            dist_est(:,c+1) = xx > 0;
            count_est(:,:,c+1) = empericalCounts(xx, a);
            mag_est(:,c+1) = sum(xx,1)/d;
        end
    end
    
    max_val = max(prob_vec);
    prob_vec = prob_vec - max_val;
    prob_vec = exp(prob_vec)/(sum(exp(prob_vec)));
    
    % Obtained the raw probability vector w (as in ST paper)
    % Sort such that w1 is the highest
    % Tip to check if coordinate system rotations are correct: notice if
    % pos=1 then the two systems are the same. Then notice that rotations
    % are proportional to either +pos or -pos where appropriate.
    [~, pos] = max(prob_vec);
    w = prob_vec( mod( (1:2*d) + pos - 2, 2*d) + 1 ); % Rotate entries to new coordinate system with max being the first entry, equivalent to [prob_vec(pos:length(prob_vec)) ; prob_vec(1:pos-1)]
    S = cumsum(w);
    ii = mod(1-pos, 2*d)+1; % This is the index we are transitioning from in the new coordinate system
    v = zeros(2*d,1); % Transition probabilities from ii in the new coordinate system
    for j = 1:2*d
        delta = S(ii) - S(mod(j-2,2*d)+1) + w(1);  % Note mod(j-1,2*d)+1 = j unless j=0 in which case = 2*d. This is for S_0 = S_n.
        v(j) = max(0, min([delta, w(ii) + w(j) - delta, w(ii), w(j)]));
    end
    v = v( mod( (1:2*d) - pos, 2*d) + 1 ) /sum(v) ; % Rotate entries back to original coordinate system and normalize

    
    if prob_vec(1,1) <= 10^(-5)
        index = discretesample(prob_vec,1);
        rb_prob_vec = prob_vec;
    else
        index = discretesample(v,1);
        rb_prob_vec = v;
    end
    
    if info_on_off == 1
        dist(:,i) = dist_est*rb_prob_vec;   % getting the weighted average
        emperical_counts(:,:,i) = rb_emp_counts(count_est, prob_vec); % getting weighted average of pairwise count estimator
        mag(1,i) = mag_est*rb_prob_vec;  % getting the weighted magnetization
    end
    
    % Find next point
    xx = samples(:, i-1);  % Current point
    if index ~= 1
        for c=1:(2*d)-1 
            xx(k(c)) = -xx(k(c));
            if (c+1) == index
                break;
            end
        end
    end
    samples(:, i) = xx;
    log_likes(1,i) = f.logp(samples(:, i));

end
emp_count_gibbs = mean(emperical_counts, 3);


