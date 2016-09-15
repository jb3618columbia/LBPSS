function [samples, log_likes, i, max_val] = CMH(f,  number_fn_evals, clique_size, initial_point, varargin)

d = f.dim;
% K is the number of binary flip proposals between recorded d-dimensional samples 
% We can take a sample after evry d flips or at every flip
% To make it comparable to BPS, we will take after evry d samples

% samples: d X L matrix

if length(varargin) > 0
    K = varargin{1};
else
    K = d;
end

S = initial_point;  % Initial point
L = number_fn_evals/(clique_size*K);           % Since every fn_eval is O(2), this denotes the total number of samples
samples = zeros(d, L);
samples(:,1)=S;    

log_likes = zeros(1,L);
log_likes(1,1) = f.logp(S);
L = number_fn_evals/(clique_size*K);            % Since every fn_eval is O(2), this denotes the total number of samples. 
                                                % Clique size is the cost
                                                % per proposed flip and
                                                % clique_size*K is the cost
                                                % per sample.
indices = randsample(d,L*K,true);

max_val = 0;
for i=2:L
   
    if mod(i,1000) == 0
    end
    for k = 1:K
        
        j=indices((i-2)*K+k);
%         Using Ari's code for probability change
        qi_1 = exp(-S(j)*f.logp_change(S,j));
        
        new_point = S;
        new_point(j) = -new_point(j);
        qi = exp(f.logp(new_point) - f.logp(S));  % Slower approach for computing probability change
        
        if max_val < abs(qi_1 - qi)
           max_val = abs(qi_1 - qi);
        end
        
        if rand() < min(1,qi)
            S(j) = -S(j);    % flip it
        end
    end
    
    log_likes(1,i) = f.logp(S);
    samples(:,i)=S;        
    
end
