function [perfect time] = ProppWilson(N,beta)
% [PERFECT TIME] = PROPPWILSON(N,BETA,MODEL) produces a state of the 
%   random cluster model according exactly the right distribution using 
%	the single-bond algorithm.
%   e.g. perfect = ProppWilson(50,log(1+sqrt(2)))
%       ( log(1+sqrt(2)) \approx 0.8813736 )
%   N - number of rows in the square lattice
%   BETA - inverse temperatur times interaction strength
%   TIME - number of required time steps


seed_start = 2009;
RandStream.setDefaultStream(RandStream('mt19937ar','seed',seed_start));

NrE = 2*N^2;   % number of possible edges

finish = 0;
m = 1;
time = 2^(m-1);

%% Initial cluster 
%(omega is a lower triangle sparse matrix representing the graph)

omega = sparse([],[],[],N^2,N^2,2*N^2); % empty set

start_bottom = omega;

% full graph
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

start_top = omega;

clear omega;


%% Preparation of the coordinates
  
    iter = N^2;    % 0.5 * number of iterations in one step (N^2 for a sweep)
    right = zeros(iter,1);
    lower = zeros(iter,1);
    
    coord1 = (1:N^2)';
    
    right(ceil(coord1./N)==N) = 1;
    lower(mod(coord1,N)==0) = 1;
    
    coord2r = coord1 + (N -N^2*right);
    coord2d = coord1 + (1 -N*lower);
    
    coord = [coord1 coord2r; coord1 coord2d];

    clear iter right lower coord1 coord2r coord2d
    

%% Propp-Wilson algorithm
while finish == 0
    
    seed = seed_start - 1;
    
    fprintf('Start in -%d\n',2^(m-1));
    
    bottom = start_bottom;
    top = start_top;

        for i=2^(m-1):-1:1   
            disp(strcat('         step: ',num2str(-i)));
            if log2(i) == ceil(log2(i))
                seed = seed + 1;
                RandStream.setDefaultStream(RandStream('mt19937ar','seed',seed));
                U = rand(1,ceil(i/2)*NrE);
                ind = 0;    % index to read a random number
            end
 
            
            % Make one sweep
            for j = 1:NrE
                ind = ind+1;
                bottom = RcSingleBondStep(bottom,coord(j,:),... 
                    beta,U(ind));
                top = RcSingleBondStep(top,coord(j,:),... 
                    beta,U(ind));
            end
        end

        
        if isequal(bottom,top)
            finish = 1;
            perfect = top;
            break;
        end
        
        m = m + 1;
        time = 2^(m-1);
        
        seed_start = seed_start - 1;
end

disp('FINISH')

    
    title = sprintf('beta = %0.3f, start = -%d', beta,time); 
    
    RcPlot(perfect,title)
    

