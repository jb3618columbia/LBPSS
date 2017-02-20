function [ Z ] = truePartition(obj)
% Comptes the partition function by brute force

d = obj.dim;
weight = zeros(1,2^d);
for j=0:2^d-1
    yy = (2*de2bi(j,d) - 1)';
    weight(1,j+1) = obj.logp(yy);
end
Z = sum(exp(weight));

end

