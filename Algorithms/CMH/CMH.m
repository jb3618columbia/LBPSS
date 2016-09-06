function [samples, log_likes, i] = CMH(f,  number_fn_evals, clique_size, initial_point, varargin)

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
log_likes(1,1) = f.logp(S);
L = number_fn_evals/(clique_size*K);           % Since every fn_eval is O(2), this denotes the total number of samples
indices = randsample(d,L*K,true);

for i=2:L
   
    if mod(i,1000) == 0
    end
    for k = 1:K
        
        j=indices((i-2)*K+k);
%         j = unidrnd(d,1);
        qi = exp(-S(j)*f.logp_change(S,j));
        
        if rand() < min(1,qi)
            S(j) = -S(j);    % flip it
        end
    end
    
    log_likes(1,i) = f.logp(S);
    samples(:,i)=S;        
    
end
