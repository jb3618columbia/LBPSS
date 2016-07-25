function [sigma, M, E] = IsingHeatbath(N,beta,B,start,stps,sh,CalcE,method,U) 
% ISINGHEATBATH(N,BETA,B,START,STPS,SH,CALCE,METHOD,U) initializes and 
%   iterates an Ising array for the given values.
%   e.g. IsingHeatbath(128,log(1+sqrt(2)),0,1,1000,1,0,0)
%       ( log(1+sqrt(2)) \approx 0.8813736 )
%   N - number of rows
%   BETA - inverse temperatur times interaction strength
%   B - external field (must be 0 for the sweep algorithm)
%   START - 0 for random choice 
%           1 for all spins up
%          -1 for all spins down
%   STPS - number of iterations
%   SH - 1 for showing all spin states
%        0 for showing only the last spin state
%       -1 for no output
%   CALCE - 1 for calculation of the energy (expensive)
%           0 for no calculation
%   METHOD - 0 for sweep update / 1 for random spin flip
%   U - all random numbers for the sweep algorithm (for Propp-Wilson). 
%       Size of U must equal STPS*N^2. (Absence causes generating new 
%       random numbers.)

M = 0;  %mag field
E = 0;  %energy


seed_start = 2010;
RandStream.setDefaultStream(RandStream('mt19937ar','seed',seed_start));


%% Control the input arguments
% if nargin<9
%     U = rand(stps*N^2, 1);  % definition of random numbers
% end
if N/2 ~= ceil(N/2)
    warning('MATLAB:paramAmbiguous','N must be even! Set N = N+1.');
    N=N+1;
end



%% Initial spin configuration
if start == -1
    sigma = -ones(N); E = 0;
elseif start == 1
    sigma = ones(N); E = 0;
else
    sigma=(-1).^(round(rand(N))); E = IsingEnergy(sigma);
end



%% Evolve the system for a fixed number of steps 
%stps=5000;
for i=1:stps, 
    % Make one step
    if method == 0
        if B~=0
            error('Sweep algorithm only works for B=0!')
        end
        if CalcE == 1
%         [sigma M E] = IsingHeatbathSweepStep(sigma,beta,...
%             U((i-1)*N^2+1:i*N^2),E);
        [sigma M E] = IsingHeatbathSweepStep(sigma,beta,...
            rand(N),E);
        else
%         [sigma M] = IsingHeatbathSweepStep(sigma,beta,...
%             U((i-1)*N^2+1:i*N^2));
        [sigma M] = IsingHeatbathSweepStep(sigma,beta,...
            rand(N));
        end
    else
        for temp = 1:N^2
            if CalcE == 1
%             [sigma M E] = IsingHeatbathStep(sigma,randi(N^2),beta,B,...
%                 U((i-1)*N^2+temp),E);
            [sigma M E] = IsingHeatbathStep(sigma,randi(N^2),beta,B,...
                rand(),E);
            else
%             [sigma M] = IsingHeatbathStep(sigma,randi(N^2),beta,B,...
%                 U((i-1)*N^2+temp));
            [sigma M] = IsingHeatbathStep(sigma,randi(N^2),beta,B,...
                rand());
            end
        end
    end
    % Display the current state of the system (optional) 
    if sh==1||(sh==0&&i==stps)

        title = sprintf('beta = %0.2f, M = %0.2f, E = %0.2f, B = %0.2f, i = %d', beta, M/N^2, E/N^2,B,i); 

        IsingPlot(sigma,title);
        
%         if sh&&i==stps
%             title=strcat('IsingH-',num2str(method),'_',num2str(beta),...
%                 '_',num2str(M/N^2),'_',num2str(B));
%             IsingSave(sigma,title);
%         end
    end
    if sh==0 && ~mod(i,100)
        fprintf('%d\n',i);
    end
end 


