function [M, E] = IsingWolff(N,beta,B,stps,sh,CalcE) 
% ISINGWOLFF(N,BETA,B,STPS,SH,CALCE) initializes and 
%   iterates an Ising array for the given values.
%   e.g. IsingWolff(250,log(1+sqrt(2)),0,1000,1,0)
%       ( log(1+sqrt(2)) \approx 0.8813736 )
%   N - number of rows
%   BETA - inverse temperatur time interaction strength
%   B - external field (must be 0 for the sweep algorithm)
%   STPS - number of iterations
%   SH - 1 for showing all spin states
%        0 for showing only the last spin state
%       -1 for no output
%   CALCE - 1 for calculation of the energy (expensive)
%           0 for no calculation

M = 0;  %mag field
E = 0;  %energy


%% Control the input arguments




%% Initial spin configuration
%sigma = ones(N); E = 0;
%sigma = -ones(N); E = 0;
sigma=(-1).^(round(rand(N))); E = IsingEnergy(sigma);



%% Evolve the system for a fixed number of steps 
%stps=5000;
for i=1:stps, 
    % Make one step
    if CalcE == 1
    [sigma M E] = IsingWolffStep(sigma,beta,randi(N^2),rand(4*N^2,1),E);
    else
    [sigma M] = IsingWolffStep(sigma,beta,randi(N^2),rand(4*N^2,1));
    end

    % Display the current state of the system (optional) 
    if sh==1||(sh==0&&i==stps)

        title = sprintf('beta = %0.2f, M = %0.2f, E = %0.2f, B = %0.2f, i = %d', beta, M/N^2, E/N^2,B,i); 

        IsingPlot(sigma,title);
        if sh&&i==stps
            title=strcat('IsingW-',num2str(beta),'+',num2str(M/N^2),'+',num2str(B));
            IsingSave(sigma,title);
        end
    end
    if sh==0 && ~mod(i,100)
        fprintf('%d\n',i);
    end
    
   % pause
end 


