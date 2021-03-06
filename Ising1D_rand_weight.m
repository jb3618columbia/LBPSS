classdef Ising1D_rand_weight < handle  % inherit from handle so that we can pass by reference

        % this class implements a 1d Ising model with zero magnetic field. 
     
        properties(SetObservable = true)
        % These properties are public by default
           d;
           beta;
           dim;
           mlp=-Inf;   
           M;
           Neis;
        end
        
             
        methods
        % These methods are public by default. 
        
            function obj = Ising1D_rand_weight(d,temp)
            % class constructor            
            
                obj.d = d;    %linear dimension of the 1D grid
                obj.beta=1/temp;     
                obj.dim = d*1;
                obj.M = zeros(obj.d,obj.d);
                obj.Neis = zeros(obj.d,2); 
                
                    for j=1:d                            
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
                    obj.M = obj.beta*obj.M/2;    % the factor of 1/2 is because pairs are counted twice                

                function nei = neighbors(j)
                    if j == 1
                        nei = [d, 2];
                    elseif j == d
                        nei = [d-1, 1]; 
                    else
                        nei = [j-1, j+1];
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
                 % For fixed weights
%                  nei = obj.Neis(j,:);
%                  lpc = 2*obj.beta*sum(S(nei));
                 
                 nei = obj.Neis(j,:);
                 weight = -2*obj.M(j,nei);
                 lpc = 2*(sum(weight'.*S(nei)));  
                                                   
             end
             
             
 
             
        end
        
        
        
end