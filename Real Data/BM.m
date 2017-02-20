classdef BM < handle  % inherit from handle so that we can pass by reference
    
    % A fully connected Boltzman machine object
    % Inputs: a symmetric weights matrix and optional bias terms.
    
    % This class is created for every outer MH iteration.
    
    properties(SetObservable = true)
        % These properties are public by default
        dim;   % will have a d x d grid
        M;
        mlp=-Inf;
        Neis;
        bias;
    end
    
    
    methods
        % These methods are public by default.
        
        function obj = BM(d,W,b)
            % class constructor
            
            obj.dim = d;    %linear dimension of the 1D grid
            obj.M = W;
            obj.bias = b;
            
        end
        % Conventions:
        % S is a signs vector  (+1,-1)
        % H = -log P(S) = S'*obj.M*S
        
        function lp = logp(obj,S)
            lp = - obj.bias*S - S'*obj.M*S; % half is added to compensate for double counting
            if lp > obj.mlp
                obj.mlp = lp;
            end
            
        end
        
        function lpc= logp_change(obj,S,j)
            % returns the difference in the log probability when
            % S(j) == +1 and S(j) == -1
            if S(j) == 1
                a = obj.logp(S);
                S(j,:) = -S(j,:);
                b = obj.logp(S);
                lpc = a-b;
            else
                b = obj.logp(S);
                S(j,:) = -S(j,:);
                a = obj.logp(S);
                lpc = a-b;
            end
            
        end
        
    end
    
end