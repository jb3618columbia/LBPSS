function [ xx, dist_est, cur_log_like, curr_fn_evals ] = analytic_slice( xx, prior, f, cur_log_like, clique_size)
% This sampler analytically finds all the coordinate flips that are above a
% certain threshold level and then uniformly samples one among them

%  [xx, cur_log)lik] = analytic_slice(xx, prior, log_lik_fn, cur_log_like)
%  There is no angle range here as we will search over the entire slice 

%  Inputs:
%  Assume that the current point xx and the prior are all of D X 1 size
%    1) Current point xx
%    2) Sample from the prior 
%    3) Log likelihood function - for computing llik for proposed points:
%    here we pass Ising 1D object which specifies this function
%    4) Log ikelihood of the current point 
%    
% Outputs:
%    1) New point xx
%    2) Likelihood of the new point
%    3) Estimate of the node marginals (from the curent set of slices)
%    4) Total number of function evaluations
    


% We will use this for 2 purposes:
% 
% 1) Uniform prior
% 2) LBP prior


d = numel(xx);
if numel(prior) ==  d
    z=reshape(prior, d,1);
else 
    error('check that prior is of correct densions')
end

hh = log(rand) + cur_log_like; % likelihood threshold
curr_fn_evals = 0;

t = (1/(d+1)):(1/(d+1)):1-(1/(d+1));
k=randperm(d);
t=t(k)';                               
w=1/((d+1));

% To be used later
% Random permutations using correlations
% t = (1/(d+1)):(1/(d+1)):1-(1/(d+1));
% sigma = abs(f.M/f.beta) + diag(ones(d,1));
% mu = zeros(d,1);
% vec = mvnrnd(mu, sigma);
% [~, k] = sort(vec);
% t=t(k)';
% w=1/((d+1));


theta=rand*w; % Draw an initial right hand side boundary
theta_accept = zeros(d+1,1);
theta_accept(1,1) = theta;
acc_samples(:,1) = xx;
j=2;
for i=2:d+1
   
    phi = t - (theta + ((i-1)*w));
    positive = phi >=0;
    negative = phi < 0;
    
    xx_prop = xx.*positive + z.*negative;
    cur_log_like = logp(f,xx_prop);
    
    if xx(i-1) ~= z(i-1)
        curr_fn_evals = curr_fn_evals + clique_size;
    end
    
    if cur_log_like > hh
        theta_accept(i,1) = theta + ((i-1)*w);
        acc_samples(:,j) = xx_prop;
        j = j + 1;
    end
    
       
end

% propose uniformly from the acceptable region
% acc_samples
theta_accept = theta_accept(theta_accept~=0);
k=numel(theta_accept);
theta_prop = theta_accept(unidrnd(k),1);    
phi_prop = t-theta_prop;
positive = phi_prop >=0;
negative = phi_prop < 0;

xx_prop = xx.*positive + z.*negative;
cur_log_like = logp(f,xx_prop);

if cur_log_like < hh
    error('Bug in this code: this point should have been accepted')
end

xx=xx_prop;

% More information from the slices
% Value of any fucntion for the slices 

% Getting the estimate of distribution 
dist_est = emp_dist(acc_samples(:,1:end));

end

