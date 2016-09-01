function [ marginal ] = JunctionTreeMarginals( a , b )
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

% Initialize junction tree factor node
F = zeros(2^n,2^n);
for s1_dec = 1:2^n             % These are the elements of the first column in the Ising model in decimal format.
    for s2_dec = 1:2^n         % These are the elements of the non-first column in the junction tree factor node. Initially this will be the second column.
        F(s1_dec, s2_dec) = decimals2weights( s1_dec, s2_dec, 1, b, n ) * decimal2bias( s1_dec, 1, a, n );
    end
end

% Iterate eliminating column from Ising model
for j = 2:(m-1)
    G = zeros(2^n,2^n);
    for s1_dec = 1:2^n             % These are the elements of the first column in the Ising model
        for s_next_dec = 1:2^n-1         % These are the elements of the next column in Ising model to be added to the junction tree factor node
            for s_cur_dec = 1:2^n         % These are the elements of the curent (non-first) column in the junction tree factor node
                G(s1_dec, s_next_dec) = G(s1_dec, s_next_dec) ...
                                + F(s1_dec, s_cur_dec) ...
                                * decimals2weights( s_cur_dec, s_next_dec, j, b, n ) ...
                                * decimal2bias( s_next_dec, 1, a, n );
            end
        end
    end
    F = G;
end

% When only two columns left, sum out last column
f = zeros(2^n);
for s1_dec = 1:2^n             % These are the elements of the first column in the Ising model
    for s_cur_dec = 1:2^n         % These are the elements of the curent (non-first) column in the junction tree factor node
        f(s1_dec) = f(s1_dec) + F(s1_dec, s_cur_dec);
    end
end

% Find marginal of first element in column
f_minus_1 = 0;                           % p(s_{1,1} = -1)
f_plus_1 = 0;                           % p(s_{1,1} = 1)
for s1_dec_not_1 = 1:2^(n-1)             % These are the elements of the first column in the Ising model
    f_minus_1 = f_minus_1 + f(s1_dec_not_1);
    f_plus_1 = f_plus_1 + f(2^(n-1) + s1_dec_not_1);
end

marginal = f_plus_1 / (f_plus_1 + f_minus_1) ; 

end

function [ result ] = decimals2weights( s1_dec, s2_dec, j, b, n ) 
% \prod_{i1=1}^n \prod_{i2=1}^n \exp{b_{i1,j,i2,j+1} s_{i1,j} s_{i2,j+1}}
%  = \exp{  \sum_{i1=1}^n \sum_{i2=1}^n   b_{i1,j,i2,j+1} s_{i1,j} s_{i2,j+1}}
    s1_bin = 2 * de2bi(s1_dec-1,n) - 1;
    s2_bin = 2 * de2bi(s2_dec-1,n) - 1;
    result = 1;
    for s1_sign = [1,-1]
        for s2_sign = [1,-1]
            vec = b(s1_bin == s1_sign , j , s2_bin == s2_sign, j+1);
            if ~isempty(vec)
                result = result * exp(s1_sign * s2_sign * sum(sum(vec)));
            end
        end
    end
end

function [ result ] = decimal2bias( s_dec, j, a, n ) 
% \prod_{i=1}^n \exp{a_{ij} s_{ij} } = \exp{ \sum_{i=1}^n a_{ij} s_{ij} } 
    s_bin = 2 * de2bi(s_dec-1,n) - 1;
    result = exp( sum(a( s_bin == 1 , j)) - sum(a( s_bin == -1 , j)));
end