% Ising model with LBP diverging

d=6;
temp=-10*pi;
scale = 20;
is2_rand = Ising2D_rand_weight(d,temp);
is2_nonrand = Ising2D(d,temp, scale);
[lbp ,jt] = get_LBP_marginals( -is2_nonrand.bias , -is2_rand.M * 100);
RMSE = sqrt(mean((lbp-jt).^2))