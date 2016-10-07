function [samples, log_likes, i] = independence_LBP(f,  number_fn_evals, clique_size, initial_point, marginals, varargin)

% Parameters
d = f.dim;
% Here K is 1, we are taking every consecutive sample
K = 1;
L = 2*number_fn_evals/(clique_size*d); 
% The justification here is simple. For 1D Ising models, the cost is O(d),
% and for 2D Ising models, the cost is ~ 2D
% Clique size for 1D is 2 and for 2D is 4.

% Initialize point and log-likelihood
S = initial_point; 
log_p = f.logp(S) - (S==1)'*log(marginals) - (S==-1)'*log(1-marginals);

% Create data structures
samples = zeros(d,L); 
samples(:,1)=S;    
log_likes = zeros(1,L);   
log_likes(1,1) = log_p;
% indices = randsample(d,L*K,true); 

for i=2:L
   
    if mod(i,1000) == 0
    end
    for k = 1:K
        S_prop = 2*(rand(d,1)< marginals) -1;
        logp_prop = f.logp(S_prop) - (S_prop==1)'*log(marginals) - (S_prop==-1)'*log(1-marginals);
        
        if (rand < exp(logp_prop - log_p) ) % Accept and propose point
            S = S_prop;
            log_p = logp_prop;
        end
    end
    
    log_likes(1,i) = f.logp(S);
    samples(:,i)=S;        
    
end
