classdef Ising2D < handle  % inherit from handle so that we can pass by reference

        % This class implements a 1d Ising model, bias and correlation strength 
        % terms have been added to encourage multimodality 
     
        properties(SetObservable = true)
        % These properties are public by default
           d;   % will have a d x d grid
           beta;
           dim;
           mlp=-Inf;   
           M;
           Neis;
           bias;
           a;
        end
        
             
        methods
        % These methods are public by default. 
        
        function obj = Ising2D(d,temp, scale)
            % class constructor
            
            obj.d = d;    %linear dimension of the 1D grid
            obj.beta=1/temp;
            obj.dim = d*d;
            obj.M = sparse(obj.dim,obj.dim);
            obj.Neis = zeros(obj.dim,4);
            obj.bias = zeros(1,obj.dim);            
            
            for j=1:d*d
                nei = neighbors(j);
                obj.Neis(j,:) = nei;
                obj.M(j, nei) =-1;
                obj.bias(1,j) = -scale*rand;
            end
            obj.M = obj.beta*obj.M/2;
            obj.a = full(obj.M);
            
            obj.bias = obj.beta*obj.bias;
            function nei = neighbors(j)
                if j == 1
                    nei = [2, d+1, d, d*(d-1)+1];
                elseif j== d
                    nei = [1, d-1, 2*d, d*d];
                elseif j== d*(d-1) + 1
                    nei = [j-d, j+1, 1, d*d];
                elseif j== d*d
                    nei = [j-1, j-d, d, d*(d-1)+1];
                elseif j == d
                    nei = [d-1,1  ];
                elseif  j>1 && j<d
                    nei = [j-1, j+1, d*(d-1)+j, j+d];
                elseif  j>(d*(d-1)) && j<d*d
                    nei = [j-1, j+1, j- (d*(d-1)), j-d];
                elseif mod(j,d) == 0
                    nei = [j-1, j-d, j-(d-1), j+d];
                elseif mod(j,d) == 1
                    nei = [j+1, j-d, j+(d-1), j+d];
                else
                    nei = [j-1, j+1, j-d, j+d];
                end
            end
        end
    
            % Conventions:
            % S is a signs vector  (+1,-1)
            % H = -log P(S) = S'*obj.M*S 
            
             function lp = logp(obj,S)                             
                lp = - obj.bias*S -S'*obj.M*S;
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
                 
%                  nei = obj.Neis(j,:);
% %                  lpc = -2*obj.bias(j) - 4*sum(S(nei).*obj.M(j,nei)');  
%                  lpc = -2*obj.bias(j) + 1*sum(S(nei))*obj.beta;  
                     
             end

        end
    
end