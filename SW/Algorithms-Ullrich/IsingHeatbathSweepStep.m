function [sigma M E] = IsingHeatbathSweepStep(sigma,beta,U,E)

% [SIGMA M E] = ISINGHEATBATHSWEEPSTEP(SIGMA,BETA,RANDTOL,U,E) generates a 
%   new Ising state with the heat bath algorithm from state SIGMA on all 
%   components by a sweep.
%   SIGMA - generated Ising state (matrix)
%   M - magnetization of the new state
%   E - energy of the new state 
%   SIGMA - initial Ising state (matrix)
%   BETA - inverse temperatur time interaction strength (BETA >= 0)
%           (critical value = 0.8813736)
%   U - necessary random numbers in a NxN matrix (0 <= U <= 1)
%           (if U = -1 it will be generated)
%   E - initial energy (absence of E will cause no calculation)


[row col] = size(sigma);

%% Control the input arguments
CalcE = 1;

if nargin<3
    error('We need SIGMA, BETA, U!');
elseif (~isa(beta,'numeric') || (~isequal(size(sigma)>1,[1 1])) ...
            || beta < 0 )
    error('Format error!');
end
if U == -1
    U = rand(row,col);
end
if nargin<4
    CalcE = 0;
end

if size(U,1)==1 || size(U,2)==1
    U = reshape(U, row, col);   % produce matrix instead of vector
end

%% Heat bath Step

transitions=zeros(row,col);

% Change first a chessboard and the the complement
for choice = 0:1

    % Calculate number of neighbored nodes unequal to the spin of all nodes
%     neighbors = ( 1-sigma.*circshift(sigma, [ 0 1]) + ... 
%         1-sigma.*circshift(sigma, [ 0 -1]) + ... 
%         1-sigma.*circshift(sigma, [ 1 0]) + ... 
%         1-sigma.*circshift(sigma, [-1 0]) )./2; 

    % Sum up all spin values of the neighbors for all nodes
    neighbors = circshift(sigma, [ 0 1]) + ... 
        circshift(sigma, [ 0 -1]) + ... 
        circshift(sigma, [ 1 0]) + ... 
        circshift(sigma, [-1 0]); 
    
    
    % Calculate the change in energy of flipping a spin 
%     DeltaE = 2 .* neighbors + 4 + 2 .* B .* sigma; 
%     DeltaE = 2*sigma .* neighbors + 2 .* B .* sigma;
%    DeltaE = 2 * sigma .* neighbors; %(for critical value for beta = 0.4406868)
%     DeltaE = sigma .* neighbors;  %(for critical value for beta = 0.8813736)
    
    % Calculate the transition probabilities 
    %p = 1./(1+exp(-beta*neighbors)); 
    p = (1+tanh((beta/2).*neighbors))./2;

    % Decide which transitions will occur 
    %transitions = (U < p ).*(rand(N) < randTol) * -2 + 1; 
    transitions = (2*(U < p)-1) .*... 
      (circshift(repmat(eye(2),ceil(row/2),ceil(col/2)),[0 choice]));

    % Calculate the change of energy
    if CalcE == 1
        Echange = sum(sum((U < p ) .*... 
        (circshift(repmat(eye(2),ceil(row/2),ceil(col/2)),[0 choice])) .* ...
        DeltaE));
        E = E + Echange;
    end
    % Perform the transitions 
        sigma = sigma.*...
            (circshift(repmat(eye(2),ceil(row/2),ceil(col/2)),[0 1-choice]))...
            + transitions;  
        
end


%% Calculate magnetization M and energy E of the new state
M = sum(sum(sigma)); 
%E = IsingEnergy(sigma); 

