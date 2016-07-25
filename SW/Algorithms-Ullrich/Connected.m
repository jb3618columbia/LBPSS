function conn = Connected(G,v,w)
% CONN = CONNECTED(G,V,W) checks, if the vertices V and W are connected 
%   in the Graph G. The graph is represented as a lower triangle sparse 
%   matrix of size N^2xN^2 (N^2 is the number of vertices), which have 
%   nonzero entries for edges in the graph. For details see "Graph Theory" 
%   in the MatLab Help.
% CONN - Is V and W connected or not (1 / 0)
% G - Graph (sparse matrix)
% V,W - Vertices, represented as number in {1,...,N}

% I take a method similar to Wolff's algorithm.

N = sqrt(length(G));

P = [v];
totalsteps = 0;
C = sparse(N,N); C(v) = 1;
step = 0;

conn = 0;

while numel(P) && conn == 0
    step = step+1;
    
    s = P(1);       % current coordinate in {1,N^2}
    k(1) = mod(s-1,N)+1;    % row of current coordinate
    k(2) = ceil(s/N);       % column of current coordinate

% Go over all neighbors
    
    % up
    l(1) = mod(k(1)-2,N)+1; 
    l(2) = k(2);
    u = (l(2)-1)*N+l(1);
    if ~(C(l(1),l(2))) && G(max(s,u),min(s,u))
        totalsteps = totalsteps + 1;
        if u == w
            conn = 1;
        end
        P = [P u];
        C(l(1),l(2)) = 1;
    end
    
    % right
    l(1) = k(1); 
    l(2) = mod(k(2),N)+1;
    u = (l(2)-1)*N+l(1);
    if ~(C(l(1),l(2))) && G(max(s,u),min(s,u))
        totalsteps = totalsteps + 1;
        if u == w
            conn = 1;
        end
        P = [P u];
        C(l(1),l(2)) = 1;
    end
    
    % down
    l(1) = mod(k(1),N)+1; 
    l(2) = k(2);
    u = (l(2)-1)*N+l(1);
    if ~(C(l(1),l(2))) && G(max(s,u),min(s,u))
        totalsteps = totalsteps + 1;
        if u == w
            conn = 1;
        end
        P = [P u];
        C(l(1),l(2)) = 1;
    end
    
    % left
    l(1) = k(1); 
    l(2) = mod(k(2)-2,N)+1;
    u = (l(2)-1)*N+l(1);
    if ~(C(l(1),l(2))) && G(max(s,u),min(s,u))
        totalsteps = totalsteps + 1;
        if u == w
            conn = 1;
        end
        P = [P u];
        C(l(1),l(2)) = 1;
    end
    
    P = P(2:numel(P));

end
