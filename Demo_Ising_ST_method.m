% Implementation of different algorithms for 1D and 2D Ising Models
% Path 
% Path for server '/home/jalaj/Github_LBPSS/Outputs_after_NIPS/data/'

%addpath(genpath('/Users/Jalaj/Documents/Github - LBPSS'));
% path = '/Users/Jalaj/Documents/Github - LBPSS/Outputs_AISTATS_final/Ising_ST_data/';
% path = '/Users/Francois/Documents/LBPSS_results/';

% Model:
ising_1d = 0;
ising_2d = 1;

%Parameters
d=81;
temp_vec=[1];
scale_vec = 0;
scale_conn = [0.25, 0.5, 0.75, 1, 2, 3, 4, 5];
number_samples = 1000;
num_examples = 15;
% rng(50)

% Algorithms:
ana_gibbs_st = 1;
ana_gibbs_rb_lbp_st = 1;
info_on_off = true;
% legend_lab = {'MH-LBP','HMC','CMH','CMH + LBP','AAS','AAS+RB','AAS+RB+LBP','AAG','AAG+RB','AAG+RB+LBP'};
% legend_lab = {'HMC','CMH','CMH + LBP','AAS','AAS+RB','AAS+RB+LBP','AAG','AAG+RB','AAG+RB+LBP'};


% This computes the errors in node and pairwise marginals
node_marginals = 1;
num_alg = 2;
marg_array = zeros(num_alg+1,2,size(scale_conn,2), size(scale_vec,2)); % for each bias and connection, store mean and std. deviation.
% +1 as we are also showing the LBP approximation
plot_marginals = 1;
pair_marg = 1;
pw_marg_array = zeros(num_alg,2,size(scale_conn,2), size(scale_vec,2));
plot_marg_pairwise = 1;


for u=1:1:length(temp_vec)
    temp = temp_vec(u);
    
    for v=1:1:length(scale_vec);
        scale_bias = scale_vec(v);
        
        for w = 1:1:length(scale_conn)
            scale_corr = scale_conn(w);
            
            if ising_1d == 1
                disp('1D Ising models');
                is1 = Ising1D_new(d,temp,scale_corr, scale_bias);  % Create 1D Ising Object
                a = linspace(1,d,d); % (i,j) for getting the pairwise marginals
                clique_size=2;
            end
            
            if ising_2d == 1
                disp('2D Ising models');
                is1 = Ising2D(sqrt(d),temp, scale_corr, scale_bias);  % Create 2D Ising Object
                [row,col] = find(abs(tril(is1.afull)) > 0); % (i,j) pairs
                aa = [col,row];
                a = reshape(aa',1,size(aa,1)*size(aa,2));
                clique_size=4;
            end
            
            % a has the pairs for which we compute the RMSE of pairwise
            % node marginals
            
            [edgeStruct, dist_LBP, edgeBelLBP, Z_LBP, marginalsJT, edgeBelJT, Z_JT] = get_LBP_marginals( -is1.bias , -is1.M);
            % initial_point = sign(normrnd(0,1,d,1));
            % Choosing the initial point from the LBP prior
            initial_point = binornd(ones(d,1), dist_LBP);
            initial_point( initial_point==0 )=-1;
            
            P = size(a,2)/2;
            true_marg = zeros(P, 4); % pairwise marginals
            for p=1:P
                true_marg(p,:) = reshape(edgeBelJT(:,:,find(ismember(edgeStruct.edgeEnds, [a(2*p-1), a(2*p)] ,'rows'))),1,4);
                % Reshape operation here gives the answer for [(1,1), (-1,1), (1,-1), (-1,-1)
            end
            
            
            dist_truth = marginalsJT; % LBP marginals
            % Error for various samplers
            
            % Node marginals
            error_ana_gibbs_st = zeros(num_examples,1);
            error_ana_gibbs_rb_lbp_st = zeros(num_examples,1);
            
            % Pairwise node marginals
            err_pw_ana_gibbs_st = zeros(num_examples,1);
            err_pw_ana_gibbs_rb_lbp_st = zeros(num_examples,1);
            
            for q=1:num_examples
                q
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Exact-HMC to decide the computational cost
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                t = 1.5;
                if ising_2d ==1
                    t = 2.5;
                end
                fn_evlas_hmc = number_samples*((is1.dim)*t + (is1.dim) + clique_size*(is1.dim)*(t-0.5));
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % AAST + RB
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if ana_gibbs_st == 1
                    tic
                    disp('AAST')
                    [samples_ana_gibbs, dist_ana_gibbs_rb, mag_ana_gibbs, loglik_ana_gibbs, nu_samples_ana_gibbs, emp_count_ana_gibbs, emp_counts_gibbs] = analytic_gibbs_new_ST_2( is1, fn_evlas_hmc, clique_size, info_on_off, initial_point, a);
                    toc
                    % doing the calculations through samples to avoid RB
                    if node_marginals==1
                        error_ana_gibbs_st = sqrt(mean( (dist_truth - emp_dist(samples_ana_gibbs(:,1:end))) .^2));
                    end
                    if pair_marg==1
                        err_pw_ana_gibbs_st(q,1) = sqrt(mean((sum(abs(empericalCounts(samples_ana_gibbs, a) - true_marg), 2)).^2));
                    end
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % AAST + RB + LBP 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if ana_gibbs_rb_lbp_st == 1
                    tic
                    disp('AAST + RB + LBP')
                    %[samples_ana_gibbs, dist_ana_gibbs_rb, mag_ana_gibbs, loglik_ana_gibbs, nu_samples_ana_gibbs, emp_count_ana_gibbs, emp_counts_gibbs] = %
                    [samples_ana_gibbs_rb_lbp_st, dist_ana_gibbs_rb_lbp_st, mag_ana_gibbs_rb_lbp_st, loglik_ana_gibbs_rb_lbp_st, nu_samples_ana_gibbs_rb_lbp_st, emp_count_ana_gibbs_rb_lbp_st, emp_counts_gibbs_rb_lbp_st] = Stretched_analytic_gibbs_ST( is1, fn_evlas_hmc, clique_size, info_on_off, initial_point, dist_LBP, a);
                    toc
                    if node_marginals==1
                        error_ana_gibbs_rb_lbp_st = sqrt(mean(  (dist_truth - mean(dist_ana_gibbs_rb_lbp_st(:,1:end),2))  .^2));
                    end
                    if pair_marg==1
                        err_pw_ana_gibbs_rb_lbp_st(q,1) = sqrt(mean((sum(abs(emp_count_ana_gibbs_rb_lbp_st - true_marg),2)).^2));
                    end
                end
                
            end   % End for the number of examples
            
            % Node marginals
            if  plot_marginals==1
                mean_err = [mean(error_ana_gibbs_st), mean(error_ana_gibbs_rb_lbp_st)];
                std_err =  [std(error_ana_gibbs_st), std(error_ana_gibbs_rb_lbp_st)];
                fileName = [path, 'data_err_st', num2str(v), num2str(w), '.mat'];
                mean_std = [mean_err', std_err'];
                save(fileName, 'mean_std')
            end
            
            % Pairwise node marginals
            if  plot_marg_pairwise ==1
                mean_err_pw = [mean(err_pw_ana_gibbs_st), mean(err_pw_ana_gibbs_rb_lbp_st)];
                std_err_pw = [std(err_pw_ana_gibbs_st), std(err_pw_ana_gibbs_rb_lbp_st)];
                fileName = [path, 'data_err_pw_st', num2str(v), num2str(w), '.mat'];
                mean_std = [mean_err_pw', std_err_pw'];
                save(fileName, 'mean_std')
            end
                        
        end
        
    end
    
end








