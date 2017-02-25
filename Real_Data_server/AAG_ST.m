function [ samples, log_likes, i, Zratio_ret] = AAG_ST(f, number_fn_evals, clique_size, initial_point, W, info_on_off)

d = f.dim;
L = round(number_fn_evals/(clique_size*(2*d-1)));
samples = zeros(d, L);
samples(:,1)=initial_point;   

% Sanity check: the error should start with 0 and remain zero 
% dist(:,1) = 0.5*ones(d,1);  
% This code passes that test

Zratio = zeros(1, L);
Zratio(1,1) = exp(initial_point'*W*initial_point); 

log_likes = zeros(1,L);
log_likes(1,1) = f.logp(initial_point);

for i=2:L
    xx = samples(:, i-1);  % Current point
    k=randperm(d);
    prob_vec = zeros(2*d,1);
    prob_vec(1,1) = f.logp(xx);
    if info_on_off == 1
        Zratio_est = zeros(1,2*d);
        Zratio_est(1,1) = exp(xx'*W*xx);
    end
    p = 2;
    j=1;
    
    for c=1:(2*d)-1 
        xx(k(j)) = -xx(k(j));
        % Efficient way to compute log likeiloohs of the propsoed point
        % cur_log_like = cur_log_like + sign(xx(k(j)))*f.logp_change(xx,k(j));
        prob_vec(p,1) = f.logp(xx);  % Inefficient  
        if info_on_off == 1
            Zratio_est(1,p) = exp(xx'*W*xx);
        end
        p = p + 1;
        j = mod(c,d) + 1;
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
    index = discretesample(v,1);
    
    
%     fn_evals = fn_evals + log(d);     % Extra log d for sampling from a dicrete distribution
    
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
    if info_on_off == 1
        Zratio(1,i) = Zratio_est*prob_vec; % adding the weighted average of pairwise count estimator
    else
        Zratio(1,i) = exp(point'*W*point);
    end

end
Zratio_ret = log(mean(Zratio));
