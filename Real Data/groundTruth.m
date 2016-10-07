function [ logZratio ] = groundTruth(bm_old, bm)
% Computes the partition function ratio by brute force

Z_old = truePartition(bm_old);
Z_new = truePartition(bm);
logZratio = log(Z_old) - log(Z_new);

end

