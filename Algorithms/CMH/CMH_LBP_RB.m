function [samples, log_likes, i, dist] = CMH_LBP_RB(f,  number_fn_evals, clique_size, initial_point, marginals, varargin)

d = f.dim;

% K is the number of binary flip proposals between recorded d-dimensional samples 
% We can take a sample after evry d flips or at every flip
% To make it comparable to BPS, we will take after evry d samples

if length(varargin) > 0
    K = varargin{1};
else
    K = d;
end
basis = eye(d); % Used for basis vectors in dist computation

S = initial_point;  % Initial point
samples(:,1)=S;    
log_likes(1,1) = f.logp(S) - (S==1)'*log(marginals) - (S==-1)'*log(1-marginals);
L = number_fn_evals/(clique_size*K);           % Since every fn_eval is O(2), this denotes the total number of samples
indices = randsample(d,L*K,true);
dist(:,1) = emp_dist(initial_point); 

for i=2:L
   
    for k = 1:K
        
        j=indices((i-2)*K+k);
        qi = exp(-S(j)*f.logp_change(S,j)) * ( marginals(j) / (1-marginals(j)) )^S(j) ; % Change in log likelihood
        alpha = min(1,qi) * (marginals(j)^(S(j)==-1) * (1-marginals(j))^(S(j)==1));
        
        if k==K
            dist(:,i) = (S==1) - alpha*S(j)*basis(:,j);
        end
        
        if (rand < alpha)  % Accept and propose point
            S(j) = -S(j);    % flip it
        end
    end
    
    log_likes(1,i) = f.logp(S);
    samples(:,i)=S;
    
end
