% Test Junction Tree Marginals
 n = 2;
 m = 2;
 
 marginal = 0.75;
 a = -ones(n,m)*log(1/marginal-1)/2;
 b = rand(n,m,n,m);
 
 JunctionTreeMarginals( a , b )
 BruteMarginals( a , b , 10000)