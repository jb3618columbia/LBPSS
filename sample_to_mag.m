function [ mag ] = sample_to_mag( samples )
% Given the sample matrix, get the implied magnetization at each sample
% point

% Input: d x N, output: 1 x N
d = size(samples, 1);
mag = sum(samples,1)/d;

end

