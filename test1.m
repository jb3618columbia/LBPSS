ising_1d = 0;
ising_2d = 1;

d=25;
temp=5*pi;
scale_bias = 2;
scale_corr = 5;

is1 = Ising2D(sqrt(d),temp, scale_corr, scale_bias);  % Create 2D Ising Object
[row,col] = find(abs(tril(is1.afull)) > 0); % (i,j) pairs
aa = [col,row];
a = reshape(aa',1,size(aa,1)*size(aa,2));
clique_size=4;

[edgeStruct, dist_LBP, edgeBelLBP, Z_LBP, marginalsJT, edgeBelJT, Z_JT] = get_LBP_marginals( -is1.bias , -is1.M);

dist_LBP - marginalsJT