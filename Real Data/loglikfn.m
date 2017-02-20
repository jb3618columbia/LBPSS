function [ loglik ] = loglikfn(W_curr_mat, W_new_mat, coronary)
% W_curr_mat, W_new_mat are upper triangular
n = size(coronary,1);
loglik = 0;

for j=1:n
    s = coronary(j,:)';
    loglik = loglik + s'*(W_new_mat - W_curr_mat)*s;   
end

end

