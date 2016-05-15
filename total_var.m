function [ error ] = total_var( f, samples, logliks, true_dist )

% Input: vector 1 X L of log likelihoods, vector of size 2^d of true probabilites 
% Compute a vector of size 2^d of implied probabilites
% Output: error: vector of difference 

d = f.d;
est_dist = zeros(1,2^d);
[~,n] = size(samples);

for i=1:n
    yy = (samples(:,i) +1)/2;
    yy_1 = bi2de(yy');
    est_dist(1,yy_1+1) = logliks(1,i);
    
end
est_dist = exp(est_dist)/sum(exp(est_dist));
error = sum(abs(true_dist - est_dist));
end

