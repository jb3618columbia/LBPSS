function [ gen_data ] = genFakedata(obj)
% Comptes the partition function by brute force

d = obj.dim;
weight = zeros(1,2^d);
for j=0:2^d-1
    yy = (2*de2bi(j,d) - 1)';
    weight(1,j+1) = obj.logp(yy);
end
prob_vec = exp(weight)/sum(exp(weight));
num_samples = 200;
X = discretesample(prob_vec, num_samples);

gen_data = zeros(num_samples, d);
for i=1:size(X,2)
    gen_data(i,:) = (2*de2bi(X(1,i)-1,d)-1);
end

end

