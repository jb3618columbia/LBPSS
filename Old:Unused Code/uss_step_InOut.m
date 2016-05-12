function [ samples, loglik, number_fn_evals, sample_index ] = uss_step_InOut(f, number_fn_evals, k)

% Slice Sampling with a uniform prior
% Stepping In, Stepping Out is done

% Input 
% f is Ising 1D object
% Output
% samples - d x number_samples matrix; each column is a sample

% Only on a line is shown 

% Initialize 
d = f.dim;
fn_evals = 0;
initial_point = ones(d,1);
samples(:,1) = initial_point;
cur_log_like = logp(f, initial_point);
loglik(1, 1) = cur_log_like;
sample_index = 2;


while fn_evals < number_fn_evals
    
    if mod(sample_index, 1000) == 0
        sample_index
    end
    
    prior = -samples(:, sample_index-1, chain_index);
    [samples(:, sample_index), cur_log_like, curr_fn_evals] = dss_step_InOut( samples(:, sample_index-1), prior, f, cur_log_like, k);
    fn_evals = fn_evals + curr_fn_evals;
    
    ll = logp(f, samples(:, sample_index));
    loglik(1, sample_index) = ll;
    sample_index = sample_index + 1;
    
end

end
