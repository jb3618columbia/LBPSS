function [ xx, cur_log_like, curr_fn_evals ] = dss_step_InOut( xx, prior, f, cur_log_like, clique_size)
% Implements stepping-in, stepping-out for slicesampling on a line

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


dim = numel(xx);

if numel(prior) ==  dim
    z=reshape(prior, dim,1);
    
else 
    error('check that prior is of correct dimensions')
    
end

hh = log(rand) + cur_log_like;
curr_fn_evals = 0;


% Generate Coordinate threshholds 
% This is equivalent to choosing a random search direction in the
% underlying space.


t = (1/(dim+1)):(1/(dim+1)):1-(1/(dim+1));
k=randperm(dim);
t=t(k)';
w=1/((dim+1));

theta=rand*w; % Draw an initial right hand side boundary

for i=1:1:dim-1
   
    theta=theta+w;
    phi = t-theta;
    positive = phi >=0;
    negative = phi < 0;
        
    xx_prop = xx.*positive + z.*negative;
    cur_log_like = logp(f,xx_prop);
    
    if xx(i) ~= z(i)
        curr_fn_evals = curr_fn_evals + clique_size;
    end
    
    if cur_log_like < hh
        theta=theta-w;
        % We found the first point not on the slice, ** EXIT LOOP **
        break;
    end
    
    
    
end

% We get [0, theta] range in which all points should be accpeted with
% probabaility 1

% New point being uniformly chosen from the accpetable region
theta_prop = rand*(theta);
phi_prop = t-theta_prop;
positive = phi_prop >=0;
negative = phi_prop < 0;

xx_prop = xx.*positive + z.*negative;
cur_log_like = logp(f,xx_prop);

if cur_log_like < hh
   
    error('Bug in this code: this point should have been accepted')
    
end

xx=xx_prop;


end

