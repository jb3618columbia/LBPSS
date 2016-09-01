function [ marginal ] = BruteMarginals( a , b, repeat )
% Computes marginals of 2d Ising models using junction trees

% Likelihoods are of the form \prod_{i=1}^n \prod_{j=1}^m \exp{a_{ij} s_{ij} } \prod_{kl \neq ij} \exp{b_{ijkl} s_{ij} s_{kl}}
% where s_{ij} \in \{+1,-1\}
% Returns marginal vector of coordinate of index 1,1 (i.e. top left in Ising model)

% The algorithm finds the marginal using a junction tree. Their is a
% factor in the junction tree between every pair of consecutive columns.
% The factors are eliminated one by one until only one column remains. Then
% the elements in the column are marginalized over, leaving the marginal
% behind.

n = size(a,1);              % Number of rows in Ising model
m = size(a,2);              % Number of columns in Ising model

% Find marginal of first element in column
f_minus_1 = 0;                           % p(s_{1,1} = -1)
f_plus_1 = 0;                            % p(s_{1,1} = 1)

for k = 0:(2^(m*n)-1)
    s = generateS(k,n,m);
    if s(1,1) == -1;
        f_minus_1 = f_minus_1 + evaluate(s,a,b,m,n);
    else
        f_plus_1 = f_plus_1 + evaluate(s,a,b,m,n);
    end
end

marginal = f_plus_1 / (f_plus_1 + f_minus_1) ;

end


function result = evaluate(s,a,b,m,n)
    result = exp(sum(sum(s.*a)));
    for j = 1:m
        for i = 1:n
            for i2 = 1:n
                result = result * exp( s(i,j) * s(i2,j) * b(i,j,i2,j) +  s(i,j) * s(i2,mod(j,m)+1) * b(i,j,i2,mod(j,m)+1) );
            end
        end
    end
end

function s = generateS(k,n,m)
    s = zeros(n,m);
    s_dec = 2*de2bi(k,n*m)-1;
    for i = 1:m
        s(:,i) = s_dec( (1+(i-1)*n):(i*n) );
    end
end