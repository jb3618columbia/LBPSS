function [samples, log_likes, total_samples] = CMH_LBP_DES_DIRTY(f,  number_fn_evals, clique_size, initial_point, marginals, varargin)

d = f.dim;

% K is the number of binary flip proposals between recorded d-dimensional samples 
% We can take a sample after evry d flips or at every flip
% To make it comparable to BPS, we will take after evry d samples

if length(varargin) > 0
    K = varargin{1};
else
    K = d;
end


S = initial_point;  % Initial point
samples(:,1)=S;    
log_likes(1,1) = f.logp(S) - (S==1)'*log(marginals) - (S==-1)'*log(1-marginals);
L = number_fn_evals/(clique_size*K);           % Since every fn_eval is O(2), this denotes the total number of samples
indices = randsample(d,L*K,true);
it_saved = 0;

i = 1;
new_points = 0;
while new_points < L*K
    i=i+1;
    
    j = unidrnd(d,1);
    qi = exp(-S(j)*f.logp_change(S,j)) * ( marginals(j) / (1-marginals(j)) )^S(j) ; % Change in log likelihood

    if (rand < (marginals(j)^(S(j)==-1) * (1-marginals(j))^(S(j)==1)))
        new_points = new_points + 1; % Count new point
        if (rand < min(1,qi)) % Accept and propose point
            S(j) = -S(j);    % flip it
        end
    end
    
    % Add every K-th sample
    if mod(i,K)==0
        log_likes(1,i/K) = f.logp(S);
        samples(:,i/K)=S;
    end
end
total_samples = length(log_likes);