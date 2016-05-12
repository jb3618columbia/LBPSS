function [ samples, dist, loglik, number_fn_evals, sample_index ] = ussSampler(f, line_on_off, analytic_on_off, number_fn_evals, k, initial_point)

% Slice Sampling with a uniform prior 

% Input 
% f is Ising 1D object

% Output
% samples - d x number_samples matrix; each column is a sample

% Initialize 
fn_evals = 0;
samples(:,1) = initial_point;
dist(:,1) = emp_dist(initial_point);
cur_log_like = logp(f, initial_point);
loglik(1, 1) = cur_log_like;
sample_index = 2;

while fn_evals < number_fn_evals
   
    if mod(sample_index, 1000) == 0
        sample_index
    end
    
%     prior = samplePrior(d);
    prior = -samples(:, sample_index-1); % Valid Sampler
    if analytic_on_off ==1
        [samples(:, sample_index), dist(:, sample_index), cur_log_like, curr_fn_evals] = analytic_slice( samples(:, sample_index-1), prior, f, cur_log_like, k);
    else
        
        if line_on_off == 1
            [samples(:, sample_index), cur_log_like, curr_fn_evals] = discrete_slice( samples(:, sample_index-1), prior, f, cur_log_like, k);
        else
            [samples(:, sample_index), cur_log_like, curr_fn_evals] = discrete_ess( samples(:, sample_index-1), prior, f, cur_log_like, k);
        end
    end
    fn_evals = fn_evals + curr_fn_evals;
    
    ll = f.logp(samples(:, sample_index));
    loglik(1, sample_index) = ll;
    sample_index = sample_index + 1;
    
    
end

end


 