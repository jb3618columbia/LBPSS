function [sigma nr] = RcToIsing(omega, U, sh)
% [SIGMA NR] = RCTOISING(OMEGA,U,SH) produce an Ising configuration SIGMA 
%   from a Random-cluster state OMEGA with the random numbers U. 
%   Every connected component of OMEGA get the same spin with probability
%   1/2.
%   The second output NR is the number of connected components in OMEGA.
%   Absence of U will cause the generation by the time.
%   (U can be absent, but SH present anyway.)
%   If SH==0 or SH is missing, there will be no plot of SIGMA.


%% Control the input arguments
if nargin<2
%     disp('Generation of U in time!');
    sh = 0;
elseif nargin<3 && length(U) == 1 && round(U)==U
%     disp('Generation of U in time!');
    sh = U;
elseif nargin<3
    sh = 0;
else
    sh = 1;
end


%% Function
N = length(omega); % number of vertices

sigma = zeros(sqrt(N));

[nr Comp] = graphconncomp(omega,'Directed',false);

if nargin<2 || (nargin<3 && length(U) == 1 && round(U)==U)
    U = 2*randi(2,nr,1)-3;
elseif length(U)<N
    error('Length of U must be the number of vertices!')
end

for i = 1:N
    sigma(i) = U(Comp(i));
end

if sh == 1
    title = sprintf('%dx%d Ising model generated from OMEGA',sqrt(N),sqrt(N));
    IsingPlot(sigma,title)
end
