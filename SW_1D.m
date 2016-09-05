function [ samples, dist, log_likes, l ] = SW_1D( is_obj, number_fn_evals,  initial_point)
% This implements SW algorithm for 1D Ising model systems:
% general connection matrix and bias vector 

% Input: 1D Ising object, starting point, 

d = is_obj.dim;
samples(:,1)=initial_point; 
log_likes(1,1) = is_obj.logp(initial_point);
fn_evals = 0;
l=2;

while fn_evals <= number_fn_evals
    curr_point = samples(:,l-1);
    bond = zeros(d,d);
    prob_matrix = 1 - exp(-1*abs(is_obj.M));
    bond_prob = (rand(d) < prob_matrix);
    
    % function to form the bonds between nodes
    % done by wrapping around
    counter=1;
    clust = zeros(1,d);
    for i=1:d 
        clust(i) = counter;
        j = is_obj.Neis(i,2);
        bond(i,j) = (curr_point(i) ==  curr_point(j))*bond_prob(i,j);
        if bond(i,j) == 1
            clust(j) = counter;
        else
            counter=counter+1;
        end
    end
    flip = unidrnd(2,1,max(clust));
    flip(flip==2)=-1;
    for k=1:max(clust)
       clust(clust==k) = flip(k); 
    end
    new_point = curr_point'.*clust;
    samples(:,l) = new_point';
    l=l+1;
    log_likes(l,1) = is_obj.logp(new_point');
    fn_evals = fn_evals + d;  % The computational cost here is O(d) 
end
dist = emp_dist(samples);                           


