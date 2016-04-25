function [ samples, loglik, number_fn_evals ] = uss_step_InOut(f, number_samples, d, number_chains, line_on_off)

% Slice Sampling with a uniform prior
% Stepping In, Stepping Out is done

% Input 
% f is Ising 1D object
% Output
% samples - d x number_samples matrix; each column is a sample

% On a line as well as on a circle are shown

% Initialize 
samples = zeros(d, number_samples, number_chains);
% energy = zeros(1, number_samples, number_chains);
loglik = zeros(1, number_samples, number_chains);
number_fn_evals = 0;

for chain_index = 1: number_chains
    
    initial_point = ones(d,1);
    samples(:,1, chain_index) = initial_point;
    cur_log_like = logp(f, initial_point);
    
    for sample_index = 2 : number_samples
        
        sample_index
        
        prior = samplePrior(d);
        if line_on_off == 1
            [samples(:, sample_index, chain_index), cur_log_like, curr_fn_evals] = dss_step_InOut( samples(:, sample_index-1, chain_index), prior, f, cur_log_like);
        else
            [samples(:, sample_index, chain_index), cur_log_like, curr_fn_evals] = discrete_ess( samples(:, sample_index-1, chain_index), prior, f, cur_log_like);
        end
        number_fn_evals = number_fn_evals + curr_fn_evals;
        
        ll = logp(f, samples(:, sample_index, chain_index));
        loglik(1, sample_index, chain_index) = ll;                
        
        
    end
    
end


end