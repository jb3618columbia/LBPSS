function [samples, log_likes, total_samples, max_val] = CMH_LBP_DES(f,  number_fn_evals, clique_size, initial_point, marginals, varargin)

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
L = number_fn_evals/(clique_size*K);           % Since every fn_eval is O(2), this denotes the total number of samples, since we are taking a sample after K proposals 
% samples = zeros(d,L*K);
samples(:,1)=S;    
log_likes(1,1) = f.logp(S) - (S==1)'*log(marginals) - (S==-1)'*log(1-marginals);

i=0;
max_val = 0;

for new_samples = 2:L*K
   
    % Sample number of rejected flips and add those counts to the previous
    % sample
    for j = 0:geornd( sum( marginals.^(S==-1) .* (1-marginals).^(S==1) ) / d )
        i=i+1;
        
        % Add every K-th sample
        if mod(i,K)==0
            log_likes(1,i/K) = f.logp(S);
            samples(:,i/K)=S;
        end 
    end
    
    % Sample new coordinate to be flipped
    j = randsample(d,1,true, marginals.^(S==-1) .* (1-marginals).^(S==1));
    
    new_point = S;
    new_point(j) = -new_point(j);
    r = exp(f.logp(new_point) - f.logp(S));
    % Change in log likelihood
    q_j = r*(marginals(j) / (1-marginals(j)) )^S(j) ; 

    % Ari's code    
    q_j_1 = exp(-S(j)*f.logp_change(S,j)) * ( marginals(j) / (1-marginals(j)) )^S(j) ; 

    if max_val < abs(q_j - q_j_1)
           max_val = abs(q_j - q_j_1);
    end
        
    % Accept the point with standard metropolis probabilities
    if (rand < min(1,q_j))
        S(j) = -S(j);    % flip it
    end       
    
end

total_samples = length(log_likes);
   



