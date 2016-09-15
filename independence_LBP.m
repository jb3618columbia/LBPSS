function [samples, log_likes, i] = independence_LBP(f,  number_fn_evals, clique_size, initial_point, marginals, varargin)



% Parameters
d = f.dim;
if length(varargin) > 0
    K = varargin{1};
    % K is the number of binary flip proposals between recorded d-dimensional samples 
    % We can take a sample after every d flips or at every flip
    % To make it comparable to BPS, we will take after every d samples
else
    K = d;
end
L = number_fn_evals/(clique_size*K);           % Since every fn_eval is O(2), this denotes the total number of samples

% Initialize point and log-likelihood
S = initial_point; 
log_p = f.logp(S) - (S==1)'*log(marginals) - (S==-1)'*log(1-marginals);

% Create data structures
samples = zeros(d,L); 
samples(:,1)=S;    
log_likes = zeros(1,L);   
log_likes(1,1) = log_p;
indices = randsample(d,L*K,true);

for i=2:L
   
    if mod(i,1000) == 0
    end
    for k = 1:K
        S_prop = 2*(rand(d,1)< marginals) -1;
        logp_prop = f.logp(S_prop) - (S_prop==1)'*log(marginals) - (S_prop==-1)'*log(1-marginals);;
        
        
        if (rand < exp(logp_prop - log_p) ) % Accept and propose point
            S = S_prop;
            log_p = logp_prop;
        end
    end
    
    log_likes(1,i) = f.logp(S);
    samples(:,i)=S;        
    
end
