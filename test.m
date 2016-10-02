d=36;
temp_vec=[1];
scale_vec = linspace(0,0.15,3);
scale_conn = linspace(0.05,0.25,5);

for v=1:1:length(scale_vec);
    scale_bias = scale_vec(v);
    
    for w = 1:1:length(scale_conn)
        scale_corr = scale_conn(w);
        
        is1 = Ising2D(sqrt(d),temp, scale_corr, scale_bias);  % Create 2D Ising Object
        [edgeStruct, dist_LBP, edgeBelLBP, Z_LBP, marginalsJT, edgeBelJT, Z_JT] = get_LBP_marginals( -is1.bias , -is1.M);
        error(v,w) = sum(abs(dist_LBP - marginalsJT));
    end
end