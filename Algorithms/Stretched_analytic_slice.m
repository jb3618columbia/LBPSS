function [ samples, dist, log_likes, i, emp_count, emperical_counts ] = Stretched_analytic_slice( f, number_fn_evals, clique_size, info_on_off, initial_point, marginals, a )
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

for i=2:L
    xx = samples(:, i-1);  % Current point
    t_bracket_width = ((xx==1).*marginals+(xx==-1).*(1-marginals))*pi;
    t = (2*rand(d,1)-1).*t_bracket_width;
    t_brackets = mod([t-t_bracket_width,t+t_bracket_width],2*pi);
    [times, flip_indices_sorted] = sort([t_brackets(1:d), t_brackets(d+1:2*d)]);
    k = mod(flip_indices_sorted-1,d)+1;

    log_y = log(rand)+ log(times(1)+2*pi-times(2*d))+f.logp(xx) - (xx==1)'*log(marginals) - (xx==-1)'*log(1-marginals);% Slice height
    prob_vec = zeros(2*d,1);
    prob_vec(1,1) = 1;
    acc_samples = zeros(d,2*d);
    acc_samples(:,1) = xx;
    p = 2;
    
    for c=1:(2*d)-1 
        xx(k(c)) = -xx(k(c));
        % Efficient way to compute log likeiloohs of the propsoed point
        % cur_log_like = cur_log_like + sign(xx(k(j)))*f.logp_change(xx,k(j));
        % prob_vec(c+1,1) = ((log(times(c+1)-times(c)) + f.logp(xx) - (xx==1)'*log(marginals) - (xx==-1)'*log(1-marginals)) > log_y);  % Inefficient
        prob_vec(c+1,1) = (times(c+1)-times(c)) * ((f.logp(xx) - (xx==1)'*log(marginals) - (xx==-1)'*log(1-marginals)) > log_y);  % Inefficient
        if prob_vec(c+1,1) == 1
            acc_samples(:,p) = xx;
            p = p + 1;
        end
    end
    acc_samples( :, all(~acc_samples,1) ) = [];
    index = unidrnd(p-1);   
    samples(:,i) = acc_samples(:,index); 
    log_likes(1,i) = f.logp(samples(:,i));
    
    if info_on_off == 1
        dist(:,i) = emp_dist(acc_samples);
        emperical_counts(:,:,i) = empericalCounts(acc_samples, a); % this is now a tensor of size K/2 x 4 x L
    end
    
%     fn_evals = fn_evals + log(d);     % Extra log d for sampling from a dicrete distribution

end
emp_count = mean(emperical_counts, 3);
