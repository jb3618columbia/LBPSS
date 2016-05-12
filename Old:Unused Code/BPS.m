function [ output_samples, samples, loglik] = BPS( f, N, num_fn_evals )

% Implements Bouncy particle sampler on a hypecube centred at the origin

% Input:
% 1) f: Ising 1D object which has the log likelihood function
% 2) N: number of samples desired
% 3) d: dimension of the problem
% 
% Output
% sample: d x number_samples matrix, each column is a sample



d = f.d;
y = abs(normrnd(0,1,d,1)); % This is the initial point
s = sign(y); % This denotes the last sample, we will update it every time
y = y/norm(y,2);


v = normrnd(0,1,d,1);  % velocity vector
% v = ones(d,1);
v = abs(v/norm(v,2));

ind = y.*v > 0;
dist = ind  - y.*sign(v);
hit_time = dist./v;    % Initial hit times

i=1;
curr_fn_evals=0;
t=0;

while curr_fn_evals < num_fn_evals
    
    
[value, min_index] = min(hit_time);

% We will try to flip the coordinate corresponding to this index
% First propose the new point
xx_prop = s;
curr = s(min_index,1);
xx_prop(min_index,1) = -curr;

% We have the proposed and the current point
% Get acceptance ratio

ratio = min(1, logp(f,xx_prop)/logp(f,s));
u = rand;
t = t+value;

if (u<= ratio) % Accept flip and update the current point
    samples(:,i) = xx_prop; 
    s = xx_prop;
    t_bounce(1,i) = t;
    i=i+1;
    % Else reject flip
end

hit_time = hit_time - value;
hit_time(min_index,1) = 1/v(min_index,1);
% curr_fn_evals = curr_fn_evals + (2+log(d));
curr_fn_evals = curr_fn_evals + (2);


end

% We have the accepted flips and coresponding bounce times
% Take samples by dividing the trajectory into regular time intervals

% Divide the total running time into N equal parts
output_samples = zeros(d,N);
loglik = zeros(1,N);
for j=1:N
    
    time = j*(t/N);
    [m, n] = max(t_bounce(t_bounce < time));
    output_samples(:,j) = samples(:,n);    
    loglik(:,j) = f.logp(output_samples(:,j));
end



end

