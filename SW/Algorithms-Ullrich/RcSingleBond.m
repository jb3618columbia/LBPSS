function omega = RcSingleBond(N,beta,start,stps,sh,method,U)
% OMEGA = RCSINGLEBOND(N,BETA,START,STPS,SH, METHOD,U) initializes and 
%   iterates an Random-cluster OMEGA array for the given values.
%   e.g. RcSingleBond(64,log(1+sqrt(2)),0,1000,1,1)
%       ( log(1+sqrt(2)) \approx 0.8813736 )
%   N - number of rows (number of vertices is N^2)
%   BETA - inverse temperatur time interaction strength
%   START - 0 for empty set 
%           1 for all edges
%          -1 for random cluster (with prob. 1/2) (don't work at the moment)
%   STPS - number of iterations
%   SH - 1 for showing all Ising states at the time and the last RC state
%        0 for showing only the last RC state
%       -1 for no output
%   METHOD - 0 for random choice of edge
%            1 for sweep (deterministic from upper left)
%   U - all random numbers for the sweep algorithm (for Propp-Wilson). 
%       Size of U must equal STPS*(number of edges). 
%       (Absence causes generating new random numbers.)


%% Define waitbar
if sh ~= 1
h = waitbar(0,'step = 0','Name','Single-bond algorithm',...
            'CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)');
end

%% Initial cluster 
%(omega is a lower triangle sparse matrix representing the graph)

omega = sparse([],[],[],N^2,N^2,2*N^2); % empty set

if start == 1   % full graph
    for i = 1:N^2
        if mod(i,N) == 0 && ceil(i/N) == N      % lower right vertex
            omega(N^2,N) = 1; % right
            omega(N^2,(N-1)*N+1) = 1; % down
        elseif mod(i,N) == 0                    % lower vertex
            omega(i+N,i) = 1; % right
            omega(i,i-N+1) = 1; % down
        elseif ceil(i/N) == N                   % right vertex
            omega(i,i-((N-1)*N)) = 1; % right
            omega(i+1,i) = 1; % down
        else                                    % other vertex
            omega(i+N,i) = 1; % right
            omega(i+1,i) = 1; % down
        end
    end
end



%% Preparation of the coordinates

iter = 2*N^2;   % number of iterations in one step (2*N^2 for a sweep)

if method == 0
    right = zeros(iter*stps,1);
    lower = zeros(iter*stps,1);
    direction = randi(2,iter*stps,1)-1;  % 0 for right / 1 for down
    
    coord1 = randi(N^2,iter*stps,1);
    
    right(ceil(coord1./N)==N) = 1;
    lower(mod(coord1,N)==0) = 1;
    
    coord2 = coord1 + direction.*(1 -N*lower) + ...
        (1-direction).*(N -N^2*right);
    
    coord = [coord1 coord2];
else
    iter = N^2;
    right = zeros(iter,1);
    lower = zeros(iter,1);
    
    coord1 = (1:N^2)';
    
    right(ceil(coord1./N)==N) = 1;
    lower(mod(coord1,N)==0) = 1;
    
    coord2r = coord1 + (N -N^2*right);
    coord2d = coord1 + (1 -N*lower);
    
    coord = [coord1 coord2r; coord1 coord2d];
    iter = 2*N^2;
end



%% Control the input arguments
if nargin<7
    U = rand(iter*stps, 1);  % definition of random numbers
end



%% Evolve the system for a fixed number of steps 
for i=1:stps 
    
    % Check for Cancel button press
    if sh ~= 1
    if getappdata(h,'canceling')
        break
    end
    end
    
    % Make one step
    if method==0
        for j = 1:iter
            [omega M E] = RcSingleBondStep(omega,coord((i-1)*iter+j,:),...
                beta,U((i-1)*iter+j));
        end
    else
        for j = 1:iter
            [omega M E] = RcSingleBondStep(omega,coord(j,:),... 
                beta,U((i-1)*iter+j));
        end
    end
    
    % Display the current state of the system (optional) 
    if (sh==0&&i==stps || sh==1&&i==stps)

        title = sprintf('beta = %0.2f, i = %d', beta,i); 
       
        RcPlot(omega,title);

%         if sh&&i==stps
%             title=strcat('Rc-',num2str(N),'_',num2str(beta));
%             IsingSave(sigma,title);
%         end
    elseif sh==1
        
        title = sprintf('beta = %0.2f, i = %d', beta,i); 
        
        sigma = RcToIsing(omega);
        
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

delete(h)
