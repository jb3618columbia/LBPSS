function [sigma M E] = IsingWolffStep(sigma,beta,U1,U2,E)

% [SIGMA M E] = ISINGWOLFFSTEP(SIGMA,BETA,U1,U2,E) generates a new 
%   Ising state with Wolff's algorithm from state SIGMA.
%   SIGMA - new Ising state
%   M - magnetization of the new state
%   E - energy of the new state
%   SIGMA - initial Ising state (matrix)
%   BETA - inverse temperatur times interaction strength (BETA >= 0)
%           (critical value = 0.8813736)
%   U1 - necessary random number for the coordinate (U1 \in \{1,...,N^2\})
%   U2 - necessary random number for the update (U2 \in [0,1]^(N^2))
%           Smaller size for U2 is possible! If more random numbers are 
%           necessary, the programm stops.
%   E - energy of initial state (absence of E will cause no calculation)

% Program from: W.Krauth. Algorithms and Computations (p. 257)

[row col] = size(sigma);

%% Control the input arguments
CalcE = 1;

if nargin<4 
    error('We need SIGMA, BETA, B, U!');
elseif (~isa(beta,'numeric') || (~isequal(size(sigma)>1,[1 1])) ...
            || beta < 0 || round(U1)~=U1)
    error('Format error!');
end
if nargin<5
    CalcE = 0;
end

%% Wolff Step

% C = [U1];
%P = [U1];
%P = sparse(2,length(U2));
% P(:,1) = [mod(U1-1,row)+1; ceil(U1/row)];
P = [mod(U1-1,row)+1; ceil(U1/row)];
totalsteps = 1;
C = sparse(row,col); C(U1) = 1;
step = 0;

while numel(P)
    step = step+1;
%     k = P(1);
%     k(1) = P(1,step); 
%     k(2) = P(2,step);    
    k(1) = P(1,1); 
    k(2) = P(2,1);
    
    
%     Erklärung des Programms
% %     for l=1:numel(sigma)
%         if (l is not in C) && sigma(k)==sigma(l) && (U2(l) < 1-exp(-beta))
%             P = [P l];
%             C = [C l];
%         end
%         P = P(2:length(P));
% %     end


% Go over all neighbors
    
    l(1) = mod(k(1)-2,row)+1; 
    l(2) = k(2);
    if ~(C(l(1),l(2))) && sigma(k(1),k(2))==sigma(l(1),l(2)) && ...
                (U2(4*step-1) < 1-exp(-beta))
        totalsteps = totalsteps + 1;
%         P(:,totalsteps) = [l(1); l(2)];
        P = [P [l(1); l(2)]];
        C(l(1),l(2)) = 1;
    end
    
    l(1) = k(1); 
    l(2) = mod(k(2),col)+1;
    if ~(C(l(1),l(2))) && sigma(k(1),k(2))==sigma(l(1),l(2)) && ...
                (U2(4*step-1) < 1-exp(-beta))
        totalsteps = totalsteps + 1;
%         P(:,totalsteps) = [l(1); l(2)];
        P = [P [l(1); l(2)]];
        C(l(1),l(2)) = 1;
    end
    
    l(1) = mod(k(1),row)+1; 
    l(2) = k(2);
    if ~(C(l(1),l(2))) && sigma(k(1),k(2))==sigma(l(1),l(2)) && ...
                (U2(4*step-1) < 1-exp(-beta))
        totalsteps = totalsteps + 1;
%         P(:,totalsteps) = [l(1); l(2)];
        P = [P [l(1); l(2)]];
        C(l(1),l(2)) = 1;
    end
    
    l(1) = k(1); 
    l(2) = mod(k(2)-2,col)+1;
    if ~(C(l(1),l(2))) && sigma(k(1),k(2))==sigma(l(1),l(2)) && ...
                (U2(4*step-1) < 1-exp(-beta))
        totalsteps = totalsteps + 1;
%         P(:,totalsteps) = [l(1); l(2)];
        P = [P [l(1); l(2)]];
        C(l(1),l(2)) = 1;
    end
    
    P = P(:,2:numel(P)/2);
%     P(:,step) = [0;0];

end
% for i=1:numel(find(C))
%     sigma(C(i)) = -sigma(C(i));
% end
sigma(find(C)) = (-1) .* sigma(find(C));

%% Calculate magnetization M and energy E of the new state
M = sum(sum(sigma)); 
if CalcE ==1
    E = IsingEnergy(sigma);
end

