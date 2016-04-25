function [xx, cur_log_like, curr_fn_evals] = discrete_slice( xx, prior, f, cur_log_like, angle_range )
%  [xx, cur_log)lik] = discrete_slice(xx, prior, log_lik_fn, cur_log_like)
%  Optionally angle range can be specified too

%  Inputs:
%  Assume that the current point xx and the prior are all of D X 1 size
%    1) Current point xx
%    2) Sample from the prior 
%    3) Log likelihood function - for computing llik for proposed points:
%    here we pass Ising 1D object which specifies this function
%    4) Log ikelihood of the current point 
%    5) Optionally: specfy an angle range
%    
% Outputs:
%    1) New point xx
%    2) ikelihood of the new point


% We will use this for 2 purposes:
% 
% 1) Uniform prior
% 2) LBP prior
% ? Stepping-in, stepping-out in theta space


   
dim = numel(xx);

if numel(prior) ==  dim
    z=reshape(prior, dim,1);
    
else 
    error('check that prior is of correct dimensions')
    
end



if (nargin < 5) || isempty(angle_range)
    angle_range = 0;
end
   
if angle_range==0
   theta = rand;
   theta_min = 0;
   theta_max = theta;

else 
    theta_min = angle_range*rand;
    theta_max = theta_min + angle_range;
    theta = rand*(theta_max - theta_min) + theta_min;
end


hh = log(rand) + cur_log_like;
curr_fn_evals = 0;

% Generate Coordinate threshholds

t = rand(dim,1);


% Slice Sampling Loop

while true 
   
    phi = t-theta;
    positive = phi >=0;
    negative = phi < 0;
    
    xx_prop = xx.*positive + z.*negative;
    coordinate_change = sum(xx_prop == xx); 
    cur_log_like = logp(f,xx_prop);
    curr_fn_evals = curr_fn_evals + 2*coordinate_change;
    
    if cur_log_like > hh
%         theta
%         positive 
%         negative
        % New point is on slice, ** EXIT LOOP **
        break;
    end
    
    % Slice shrinkage step
    if theta < 0
        theta_min = theta;
    elseif theta >= 0
        theta_max = theta;
    else
        error('BUG DETECTED: Shrunk to current position and still not acceptable.');
    end
    
    % Propose a new point
    theta = rand*(theta_max-theta_min) + theta_min;
end

xx=xx_prop;
end

