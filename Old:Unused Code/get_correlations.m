function [ marginals_JT, marginal_f ] = get_correlations( a , b, k )

nStates = 2;

% Make dimensions for vector a correct
[a_rows, a_col] = size(a);
if a_col > a_rows
    a = a';
end

% Make b symmetric if not already
if ~isequal(b,b') % If not symmetric
    disp('b matrix not symmetric, making it symmetric')
    b = b+b'; % Make symmetric
end

% Find the normalizing constant
[ ~, ~, marginalsJT, log_Z_JT ] = get_LBP_marginals( a , b );
marginals_JT = marginalsJT(k);

% a_f = a;
% a_f(k) = 10^2;
% [ ~, ~, ~, Z_JT_f ] = get_LBP_marginals( a_f , b );
% marginal_f = Z_JT_f * exp(a(k)) / Z_JT;


% Find the marginal by fixing a coordinate
b_f = b( 1:length(a) ~= k , 1:length(a) ~= k );
a_f = a( 1:length(a) ~= k ).*exp(b( 1:length(a) ~= k , k));

% Find the new normalizing constant
[ ~, ~, ~, log_Z_JT_f ] = get_LBP_marginals( a_f , b_f );

% The new marginal
marginal_f = exp(log_Z_JT_f + a(k) - log_Z_JT);