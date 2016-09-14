function [ rb_count ] = rb_emp_counts(count_est, prob_vec)
 % To compute the emperical counts in the inner loop of the
 % rao-blackwelization step.
 
 % count_est is a tensor: K/2 x 4 x 2*d
 % prob_vec is a column vecor: 2*d x 1
 % count_est is returned as a K/2*4 matrix: rb estimate

  
 C = zeros(1,1,size(prob_vec,1));
 C(:) = prob_vec;
 rb_count = sum(bsxfun(@times,count_est,C),3);
 
 