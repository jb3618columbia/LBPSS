function [omega M E] = RcSingleBondStep(omega,coord,beta,U)

% [OMEGA2 M E] = RCSINGLEBONDSTEP(OMEGA,COORD,BETA,B,U) generates a new 
%   Random-cluster state with the single bond algorithm from state OMEGA 
%   on the coordinate COORD.
%   OMEGA2 - new Random cluster state
%   M - magnetization of the new state
%   E - energy of the new state
%   OMEGA - initial Random cluster state (lower triangle matrix (sparse))
%   COORD - valid coordinate to change (must be a 2D-vector)
%       Valid means: COORD(1) any vertex, COORD(2) right or lower neighbor!
%   BETA - inverse temperatur times interaction strength (BETA >= 0)
%           (critical value = 0.8813736)
%   U - necessary random number (0 <= U <= 1)


N = length(omega);    % size of connection matrix (number of vertices)


%% Control the input arguments
% CalcE = 1;

% % sort the coordinates
% x = coord;
% coord(1) = min(x);
% coord(2) = max(x);

if nargin<4
    error('We need OMEGA, COORD ,BETA, U!');
% elseif (~(isequal(size(coord),[1 2])) && ~(isequal(size(coord),[2 1])))
%     error('COORD must be a 2-dimensional vector!');
end
% if (~isa([beta U],'numeric') || (~isequal(size(omega)>1,[1 1])) ...
%             || beta < 0 || U < 0 || U > 1 )
%     error('Format error!');
% end
if coord(1)~=coord(2)-1 && coord(2)~=coord(1)-1 && ...
   coord(1)~=coord(2)-sqrt(N) && coord(2)~=coord(1)-sqrt(N) && ...
        coord(1)~=coord(2)+(sqrt(N)-1)*sqrt(N) && ...
        coord(2)~=coord(1)+(sqrt(N)-1)*sqrt(N) && ...
        coord(1)~=coord(2)+sqrt(N)-1 && coord(2)~=coord(1)+sqrt(N)-1
    error('%d and %d are no direct neighbors!',coord(1),coord(2));
end
% if nargin<5
%     CalcE = 0;
% end

%% Single bond Step

p = 1 - exp(-beta); % standard definition

% Make the step 

omega(max(coord),min(coord)) = 0;

if U < p/(2-p) % && ~Connected(omega,coord(1),coord(2))
    omega(max(coord),min(coord)) = 1;
elseif U < p 
    
%     % MatLab routine:
%     [nr_clusters Comp] = graphconncomp(omega,'Directed',false);
%     
%     if Comp(coord(1)) == Comp(coord(2))     % connected in omega\e
%         omega(max(coord),min(coord)) = 1;
%     end

%   Self implemented:
    if Connected(omega,coord(1),coord(2))     % connected in omega\e
        omega(max(coord),min(coord)) = 1;
    end
    
else
    omega(max(coord),min(coord)) = 0;
end


%% Calculate magnetization M and energy E of the new state
% if CalcE ==1 || spin ~= sigma(coord(1),coord(2))
%     E = E + DeltaE;
% end
%M = sum(sum(sigma)); 
%E = IsingEnergy(sigma); 
M=[];
E=[];

