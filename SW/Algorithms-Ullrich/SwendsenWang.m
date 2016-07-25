function [sigma omega] = SwendsenWang(N, beta, stps, start, sh)
% [SIGMA OMEGA] = SWENDSENWANG(N,BETA,STPS,START,SH) initializes and 
%   iterates an Ising and random cluster array for the given values.
%   e.g. SwendsenWang(100,log(1+sqrt(2)),1000,1,1)
%       ( log(1+sqrt(2)) \approx 0.8813736 )
%   SIGMA - last Ising configuration
%   OMEGA - last random cluster state
%   N - number of rows
%   BETA - inverse temperatur time interaction strength
%   STPS - number of iterations
%   START - For the initial Ising state:
%           1 for the "all spins up" state
%          -1 for the "all spins down" state
%           0 for random state
%   SH - 1 for showing all spin states
%        0 for showing only the last spin state
%       -1 for no output


%% Define waitbar
if sh ~= 1
h = waitbar(0,'step = 0','Name','Swendsen-Wang algorithm',...
            'CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)');
end

%% Initial spin configuration
if start == -1
    sigma = -ones(N);
elseif start == 1
    sigma = ones(N);
else
    sigma=(-1).^(round(rand(N)));
end


%% Evolve the system for a fixed number of steps 

omega = sparse([],[],[],N^2,N^2,2*N^2); % empty set

for i=1:stps 
    
    % Check for Cancel button press
    if sh ~= 1
    if getappdata(h,'canceling')
        break
    end
    end
    
    % Make one step
    omega = IsingToRc(sigma,beta);
    sigma = RcToIsing(omega);
    

    % Display the current state of the system (optional) 
    if sh==1 || (sh==0&&i==stps)

        title = sprintf('beta = %f, i = %d', beta,i); 

        IsingPlot(sigma,title);

    end
%     if sh==0 && ~mod(i,100)
%         fprintf('%d\n',i);
%     end
    
    % Refresh the waitbar
    if sh ~= 1
        waitbar(i/stps,h,sprintf('step = %d',i))
    end
end

if sh ~= 1
    delete(h)
end
