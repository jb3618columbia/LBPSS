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
log_likes(1,1) = f.logp(S);
fn_evals = 0;
i=2;

% dist_est is a matrix of estimated distribution after each sample
dist(:,1) = emp_dist(S);
% This code pases the sanity check
% dist(:,1) = 0.5*ones(d,1);


while fn_evals <= number_fn_evals
   
    xx = samples(:, i-1);  % Current point
    z = -xx;
    
    t = (1/(d+1)):(1/(d+1)):1-(1/(d+1));
    k=randperm(d);
    t=t(k)';   % Order of coordinates 
    w=1/((d+1));
    
    theta=0.5*w; 
    prob_vec = zeros(2*d,1);
    prob_vec(1,1) = f.logp(xx);
    dist_est(:,1) = xx > 0; 
    
    for j=2:d+1
        
        phi = t - (theta + ((j-1)*w));
        positive = phi > 0;
        negative = phi < 0;
        
        xx_prop = xx.*positive + z.*negative;
        fn_evals = fn_evals + clique_size;
        if info_on_off ==1
            dist_est(:,j) = xx_prop > 0;
        end
        prob_vec(j,1) = f.logp(xx_prop);
        
    end
    
    for j=d+2:(2*d)
        
        phi = t - (theta + (j-1-d)*w);
        positive = phi >=0;
        negative = phi < 0;
        
        xx_prop = z.*positive + xx.*negative;
        fn_evals = fn_evals + clique_size;
        if info_on_off ==1
            dist_est(:,j) = xx_prop > 0;
        end
        prob_vec(j,1) = f.logp(xx_prop);
        
    end
    
    max_val = max(prob_vec);
    prob_vec = prob_vec - max_val;
    prob_vec = exp(prob_vec)/(sum(exp(prob_vec)));

    if info_on_off == 1
        dist(:,i) = dist_est*prob_vec;   % getting the weighted average
    end
    
%     fn_evals = fn_evals + log(d);     % Extra log d
    index = discretesample(prob_vec,1);
    
    if index > d+1
        phi = t - (theta + ((index-1-d)*w) );
        positive = phi >=0;
        negative = phi < 0;
        samples(:,i) = z.*positive + xx.*negative;
        log_likes(1,i) = f.logp(samples(:,i));
    else
        phi = t - (theta + ((index-1)*w) );
        positive = phi >=0;
        negative = phi < 0;
        samples(:,i) = xx.*positive + z.*negative;
        log_likes(1,i) = f.logp(samples(:,i));
    end
    
    i = i + 1;
    
end


