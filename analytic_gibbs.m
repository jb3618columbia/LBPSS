function [ samples, dist, log_likes, i ] = analytic_gibbs( f, number_fn_evals, clique_size, info_on_off, initial_point )
% This implements analytic gibbs sampling instead of analytic slice
% sampling. We will search over all the coordinate flips proposed 

%  [xx, cur_log)lik] = analytic_gibbs(xx, prior, log_lik_fn, cur_log_like)


%  Inputs:
%  Log likelihood function - for computing llik for proposed points:
%  here we pass Ising 1D object which specifies this function

%    
% Outputs:
% A matrix of samples: d X L

% We will use it with LBP and Uniform Prior 

d = f.dim;
S = initial_point;  % Initial point
samples(:,1)=S;    
% dist_est is a matrix which gives the estimated distribution after each
% sample
dist(:,1) = emp_dist(S);
log_likes(1,1) = f.logp(S);
fn_evals = 0;
i=2;

while fn_evals <= number_fn_evals
   
    xx = samples(:, i-1);  % Current point
    z = -xx;
%     z = samplePrior(d);  % Use LBP or uniform prior    
    
    t = (1/(d+1)):(1/(d+1)):1-(1/(d+1));
    k=randperm(d);
    t=t(k)';   % Order of coordinates 
    w=1/((d+1));
    
    theta=0.5*w; 
    prob_vec = zeros(d+1,1);
    prob_vec(1,1) = f.logp(xx);
    dist_est(:,1) = emp_dist(xx); 

    
    for j=2:d+1
        
        phi = t - (theta + ((j-1)*w));
        positive = phi > 0;
        negative = phi < 0;
        
        xx_prop = xx.*positive + z.*negative;
        if info_on_off ==1
            dist_est(:,j) = emp_dist(xx_prop);
        end
        prob_vec(j,1) = f.logp(xx_prop);
        
        if xx(j-1) ~= z(j-1)
            fn_evals = fn_evals + clique_size;
        end
      
    end
%     prob_vec;
    max_val = max(prob_vec);
    prob_vec = prob_vec - max_val;
    prob_vec = exp(prob_vec)/(sum(exp(prob_vec)));
    
    if info_on_off == 1
        dist(:,i) = dist_est*prob_vec;
    end
    
    fn_evals = fn_evals + log(d);     % Extra log d
    point = discretesample(prob_vec,1);
    phi = t - (theta + ((point-1)*w) );
    positive = phi >=0;
    negative = phi < 0;
    samples(:,i) = xx.*positive + z.*negative;
    log_likes(1,i) = f.logp(samples(:,i));
    i = i + 1;
    
    
end