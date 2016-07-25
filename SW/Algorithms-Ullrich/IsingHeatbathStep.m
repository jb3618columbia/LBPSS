function [sigma M E] = IsingHeatbathStep(sigma,coord,beta,B,U,E)

% [SIGMA M E] = ISINGHEATBATHSTEP(SIGMA,COORD,BETA,B,U,E) generates a new 
%   Ising state with the heat bath algorithm from state SIGMA on the 
%   coordinate COORD.
%   SIGMA - new Ising state
%   M - magnetization of the new state
%   E - energy of the new state
%   SIGMA - initial Ising state (matrix)
%   COORD - coordinate to change
%           (2d vector or integer 
%       ( (x,y) = x+(y-1)*rows(SIGMA) ) 
%       ( x = mod(i-1,rows(SIGMA))+1, y = ceil(i/rows(SIGMA)) ) )
%   BETA - inverse temperatur time interaction strength (BETA >= 0)
%           (critical value = 0.8813736)
%   B - external field
%   U - necessary random number (0 <= U <= 1)
%   E - initial energy (absence of E will cause no calculation)


[row col] = size(sigma);

%% Control the input arguments
CalcE = 1;

if nargin<5 
    error('We need SIGMA, COORD ,BETA, B, U!');
elseif (~(isequal(size(coord),[1 2])) && ~(isequal(size(coord),[2 1])) ...
            && ~(isequal(size(coord),[1 1])))
    error('COORD must be a 2-dimensional vector or integer!');
elseif isequal(size(coord),[1 1])
    x = coord;
    coord(1) = mod(x-1,row)+1; 
    coord(2) = ceil(x/row);
end
if (~isa([beta B U],'numeric') || (~isequal(size(sigma)>1,[1 1])) ...
            || beta < 0 || U < 0 || U > 1 )
    error('Format error!');
end
if nargin<6
    CalcE = 0;
end

%% Heat bath Step

% Calculate number of neighbored nodes unequal to the spin of the node
spin = sigma(coord(1),coord(2));
% neighbors = (spin~=sigma(mod(coord(1)-2,row)+1, coord(2))) + ...
%         (spin~=sigma(coord(1), mod(coord(2),col)+1)) + ...
%         (spin~=sigma(mod(coord(1),row)+1, coord(2))) + ...
%         (spin~=sigma(coord(1), mod(coord(2)-2,col)+1));
neighbors = sigma(mod(coord(1)-2,row)+1, coord(2)) + ...
        sigma(coord(1), mod(coord(2),col)+1) + ...
        sigma(mod(coord(1),row)+1, coord(2)) + ...
        sigma(coord(1), mod(coord(2)-2,col)+1);    

% Calculate the change in energy of flipping a spin 
DeltaE = neighbors + B; 

% Calculate the transition probabilities 
p = 1/(1+exp(-beta*DeltaE)); 
%p = (1+tanh((beta/2).*neighbors))./2;

% Make the step
if U < p
    sigma(coord(1),coord(2)) = 1;
else
    sigma(coord(1),coord(2)) = -1;
end
if CalcE ==1 && spin ~= sigma(coord(1),coord(2))
    E = E + DeltaE;
end

%% Calculate magnetization M and energy E of the new state
M = sum(sum(sigma)); 
%E = IsingEnergy(sigma); 


