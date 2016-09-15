function [ count ] = empericalCountsOld(samples, a)

% samples: matrix of MCMC samples: d X N 
% a: vector of size two, specify the 2 coordinates
% count: 4 x 1 vector giving emperical counts for [(1,1), (-1,1), (1,-1), (-1,-1)]

i = a(1); j=a(2);
count = [sum((samples(i, :) > 0) .* (samples(j, :) > 0)) , ... 
         sum((samples(i, :) < 0) .* (samples(j, :) < 0)) , ...
         sum((samples(i, :) < 0) .* (samples(j, :) > 0)) , ...
         sum((samples(i, :) > 0) .* (samples(j, :) < 0)) ] ...
         / size(samples,2) ;



m = samples(i, :) + samples(j,:);
count(1,1) = size(m(m==2),2);
count(4,1) = size(m(m==-2),2);

n = samples(i, :) - samples(j,:);
count(2,1) = size(n(n==-2),2);
count(3,1) = size(n(n==2),2);

% Now a is a vector of more than 2 elements, want to compare across
% multiple pairs
% K = size(a,2);
% for k=1:K-1
%    i = a(k); j=a(k+1);
%    m = samples(i, :) + samples(j,:);
%    count(k,1) = size(m(m==2),2);
%    count(k,4) = size(m(m==-2),2);
%    
%    n = samples(i, :) - samples(j,:);
%    count(k,2) = size(n(n==-2),2);
%    count(k,3) = size(n(n==2),2);
%    
% end

count = count/size(samples,2);
