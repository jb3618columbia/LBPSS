classdef Ising2D_rand_weight < handle  % inherit from handle so that we can pass by reference

        % this class implements a 1d Ising model with zero magnetic field. 
     
        properties(SetObservable = true)
        % These properties are public by default
           d;   % will have a d x d grid
           beta;
           dim;
           mlp=-Inf;   
           M;
           Neis;
        end
        
             
        methods
        % These methods are public by default. 
        
        function obj = Ising2D_rand_weight(d,temp)
            % class constructor
            
            obj.d = d;    %linear dimension of the 1D grid
            obj.beta=1/temp;
            obj.dim = d*d;
            obj.M = zeros(obj.d*obj.d,obj.d*obj.d);
            obj.Neis = zeros(obj.d*obj.d,4);
            
            
            for j=1:d*d
                nei = neighbors(j);
                obj.Neis(j,:) = nei;
                for k=1:length(nei)
                    if nei(k) > j
                        w = -5*rand;
                        obj.M(j,nei(k)) = w;
                        obj.M(nei(k),j) = w;
                    end
                end
            end
            obj.M = obj.beta*obj.M/2;
            
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
                lp = -S'*obj.M*S;
                if lp > obj.mlp
                    obj.mlp = lp;
                end                
                    
             end
         
             function lpc= logp_change(obj,S,j)
                 % returns the difference in the log probability when 
                 % S(j) == +1 and S(j) == -1     
                 
                 nei = obj.Neis(j,:);
                 weight = -2*obj.M(j,nei);
                 lpc = 2*sum(weight'.*S(nei));                 
                     
             end
             
             
 
             
        end
        
        
        
end