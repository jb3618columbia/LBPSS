function [ count ] = empericalCounts(samples, a)

% samples: matrix of MCMC samples: d X N 
% a: vector of size two, specify the 2 coordinates
% count: 4 x 1 vector giving emperical counts for [(1,1), (-1,1), (1,-1), (-1,-1)]

i = a(1); j=a(2);
m = samples(i, :) + samples(j,:);
count(1) = size(m(m==2),2);
count(4) = size(m(m==-2),2);

n = samples(i, :) - samples(j,:);
count(2) = size(n(n==-2),2);
count(3) = size(n(n==2),2);