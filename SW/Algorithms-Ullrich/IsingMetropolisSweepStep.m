function [sigma M E] = IsingMetropolisSweepStep(sigma,beta,randTol,U,E)

% [SIGMA2 M E] = ISINGMETROPOLISSWEEPSTEP(SIGMA,COORD,BETA,B,U) generates a 
%   new Ising state with the Metropolis algorithm from state SIGMA on all 
%   components by a sweep.
%   SIGMA - generated Ising state (matrix)
%   M - magnetization of the new state
%   E - energy of the new state 
%   SIGMA - initial Ising state (matrix)
%   BETA - inverse temperatur time interaction strength (BETA >= 0)
%           (critical value = 0.8813736)
%	RANDTOL - Determines how much steps would be accepted. This is only 
%       necessary for the sweep algorithm and is only for the look of the 
%       output. Should be 0.5, but 1 for full speed.
%   U - necessary random numbers in a NxN matrix (0 <= U <= 1 / 
%               U = -1: random numbers will be generated in time)
%   E - initial energy (absence of E will cause no calculation)


[row col] = size(sigma);

%% Control the input arguments
CalcE = 1;

if nargin<4 
    error('We need SIGMA, BETA, RANDTOL, U!');
elseif (~isa(beta,'numeric') || (~isequal(size(sigma)>1,[1 1])) ...
            || beta < 0 || ...
            ( (~isequal(size(U)==size(sigma),[1 1]))&&...
            (~isequal(size(U),[1 1])) )  )
    error('Format error!');
elseif U == -1
    U = rand(row,col);
end
if nargin<5
    CalcE = 0;
end

%% Metropolis Step

% Change first a chessboard and the the complement
for choice = 0:1

    % Calculate number of neighbored nodes unequal to the spin of all nodes
%     neighbors = ( 1-sigma.*circshift(sigma, [ 0 1]) + ... 
%         1-sigma.*circshift(sigma, [ 0 -1]) + ... 
%         1-sigma.*circshift(sigma, [ 1 0]) + ... 
%         1-sigma.*circshift(sigma, [-1 0]) )./2; 
    neighbors = circshift(sigma, [ 0 1]) + ... 
        circshift(sigma, [ 0 -1]) + ... 
        circshift(sigma, [ 1 0]) + ... 
        circshift(sigma, [-1 0]); 
    
    
    % Calculate the change in energy of flipping a spin 
%     DeltaE = 2 .* neighbors + 4 + 2 .* B .* sigma; 
%     DeltaE = 2*sigma .* neighbors + 2 .* B .* sigma;
%    DeltaE = 2 * sigma .* neighbors; %(critical value for beta = 0.4406868)
    DeltaE = sigma .* neighbors;  %(critical value for beta = 0.8813736)
    
    % Calculate the transition probabilities 
    p = exp(-beta*DeltaE); 

    % Decide which transitions will occur 
    %transitions = (U < p ).*(rand(N) < randTol) * -2 + 1; 
    transitions = (U < p ).* (rand(row,col) < randTol) .*... 
      (circshift(repmat(eye(2),ceil(row/2),ceil(col/2)),[0 choice]))* -2+1;

    % Calculate the change of energy
    if CalcE == 1
        Echange = sum(sum((U < p ).* (rand(row,col) < randTol) .*... 
        (circshift(repmat(eye(2),ceil(row/2),ceil(col/2)),[0 choice])) .* ...
        DeltaE));
        E = E + Echange;
    end
    % Perform the transitions 
        sigma = sigma .* transitions; 
end


%% Calculate magnetization M and energy E of the new state
M = sum(sum(sigma)); 
%E = IsingEnergy(sigma); 

