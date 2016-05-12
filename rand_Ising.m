classdef rand_Ising < handle  % inherit from handle so that we can pass by reference

        % this class implements a 1d Ising model with zero magnetic field. 
     
        properties(SetObservable = true)
            % These properties are public by default
            dim;
            beta;
            mlp=-Inf;
            M;
            Neis;
            
        end
        
        
        methods
            % These methods are public by default.
            
            function obj = rand_Ising(d, temp)
                % class constructor
                
                obj.dim = d;    %linear dimension of the 1D grid
                obj.beta=1/temp;
                obj.M = randi(2, obj.dim,obj.dim);
                obj.Neis = zeros(obj.dim,obj.dim);
                obj.M = obj.beta*obj.M;    % the factor of 1/2 is because pairs are counted twice             
            end
            
            function lp = logp(obj,S)
                lp = -S'*obj.M*S;
                if lp > obj.mlp
                    obj.mlp = lp;
                end
                
            end
            
            function lpc= logp_change(obj,S,j)
                 % returns the difference in the log probability when 
                 % S(j) == +1 and S(j) == -1     
                 nei = obj.M(j,:);
                 lpc = obj.dim*obj.beta*sum(S(nei));                 
                     
             end
            
        end
            
end



                
                
         
             
             