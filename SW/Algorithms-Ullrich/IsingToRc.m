function omega = IsingToRc(sigma, beta, U, sh)
% [OMEGA M E] = ISINGTORC(SIGMA,U,SH) produce a random cluster state OMEGA 
%   from a Ising configuration SIGMA with the random numbers U. 
%   (U must be a 0-1-matrix with length 2*(number of vertices))
%   Each two neighbors with equal spin get a edge between them with 
%   probability p = 1 - exp(-BETA).
%   The second output NR is the number of connected components in OMEGA.
%   Absence of U will cause the generation by the time.
%   (U can be absent, but SH present anyway.)
%   If SH==0 or SH is missing, there will be no plot of SIGMA.


%% Control the input arguments
if nargin<2
    error('We need SIGMA and BETA!');
elseif nargin<3
%     disp('Generation of U in time!');
    U = binornd(1, 1-exp(-beta), 1, 2*size(sigma,1)*size(sigma,2));
    sh = 0;
elseif nargin<4 && length(U) == 1 && round(U)==U
%     disp('Generation of U in time!');
    sh = U;
    U = binornd(1, 1-exp(-beta), 1, 2*size(sigma,1)*size(sigma,2));
else
    sh = 1;
end


%% Function

N = size(sigma,1);    % length of one site of the square lattice

omega = sparse([],[],[],N^2,N^2,2*N^2); % empty set

ind = 0;    % index for the random numbers

for i = 1:N^2
    if mod(i,N) == 0 && ceil(i/N) == N      % lower right vertex
        ind = ind + 2;
        omega(N^2,N) = 1 * (sigma(N^2)==sigma(N)) * U(ind-1); % right
        omega(N^2,(N-1)*N+1) = 1 * (sigma(N^2)==sigma((N-1)*N+1)) * U(ind); % down
    elseif mod(i,N) == 0                    % lower vertex
        ind = ind + 2;
        omega(i+N,i) = 1 * (sigma(i+N)==sigma(i)) * U(ind-1); % right
        omega(i,i-N+1) = 1 * (sigma(i)==sigma(i-N+1)) * U(ind); % down
    elseif ceil(i/N) == N                   % right vertex
        ind = ind + 2;
        omega(i,i-((N-1)*N)) = 1 * (sigma(i)==sigma(i-((N-1)*N))) * U(ind-1); % right
        omega(i+1,i) = 1 * (sigma(i+1)==sigma(i)) * U(ind); % down
    else                                    % other vertex
        ind = ind + 2;
        if ((sigma(i+N)==sigma(i)) && U(ind-1)) 
            omega(i+N,i) = 1; % right
        end
        if ((sigma(i+1)==sigma(i)) && U(ind))
            omega(i+1,i) = 1; % down
        end
    end
end
    
%% Plot
if sh == 1
    title = sprintf('BETA=%f',beta);
    RcPlot(omega,title)
end
    
end