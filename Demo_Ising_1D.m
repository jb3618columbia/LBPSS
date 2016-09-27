% Implementation of different algorithms for 1D and 2D Ising Models
% Parameters
ising_1d = 0;
ising_2d = 1;

d=4;
temp_vec=[5*pi];
scale_vec = [50];
scale_conn = [1];
% rng(50)

% This computes the errors in node and pairwise marginals
node_marginals = 1;
pair_marg = 1;

plot_marg_pairwise = 1;
plot_marginals = 1;

% This plots the log-likelihood of samples
plot_log_liks=0;

% This computes and plots the auto-correlation time for magnetization
act_mag=1;
plot_act_mag=1;


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
            
            % Using Ari's object for zero bias
            %         is1 = Ising1D(d,temp);
            %         bias = zeros(1,d);
            
            if ising_2d == 1
                disp('2D Ising models');
                is1 = Ising2D(sqrt(d),temp, scale_bias);  % Create 2D Ising Object
                [row,col] = find(abs(tril(is1.a)) > 0); % (i,j) pairs
                aa = [col,row];
                a = reshape(aa',1,size(aa,1)*size(aa,2));
                clique_size=4;
            end
            
            number_samples = 400;
            num_examples = 1;
            
            % Algorithms:
            truth = 0;     % brute force gorund truth; no longer needed
            ind_sampler = 1;
            exact_hmc = 1;
            cmh = 1;
            cmh_lbp = 1;
            ana_slice = 1;
            ana_slice_rb_lbp = 1;
            ana_gibbs = 1;
            ana_gibbs_rb_lbp = 1;
            sw = 0;
            info_on_off = true;
            
            %%%%%%%%%%%%%%%%%%%%%%%
            % True Node Marginals
            %%%%%%%%%%%%%%%%%%%%%%%
            if truth == 1
                weight = zeros(1,2^d);
                for j=0:2^d-1
                    if mod(j,10000)==0
                    end
                    yy = (2*de2bi(j,d) - 1)';
                    weight(1,j+1) = is1.logp(yy);
                end
                weight = exp(weight)/sum(exp(weight));
                dist_truth = zeros(d,1);
                for j=0:2^d-1
                    if mod(j,10000)==0
                    end
                    yy = (2*de2bi(j,d) - 1)';
                    dist_truth = dist_truth + weight(1,j+1)*(yy >0);
                end
            end
            
            %         If using zero bias as in Ari's original code
            %         [edgeStruct, dist_LBP, edgeBelLBP, Z_LBP, marginalsJT, edgeBelJT, Z_JT] = get_LBP_marginals( bias , -is1.M);
            
            [edgeStruct, dist_LBP, edgeBelLBP, Z_LBP, marginalsJT, edgeBelJT, Z_JT] = get_LBP_marginals( -is1.bias , -is1.M);
            % initial_point = sign(normrnd(0,1,d,1));
            initial_point = binornd(ones(d,1), dist_LBP);
            initial_point( initial_point==0 )=-1;
            
            P = size(a,2)/2;
            true_marg = zeros(P, 4);
            for p=1:P
                true_marg(p,:) = reshape(edgeBelJT(:,:,find(ismember(edgeStruct.edgeEnds, [a(2*p-1), a(2*p)] ,'rows'))),1,4);
                % Reshape operation here gives the answer for [(1,1), (-1,1), (1,-1), (-1,-1)
            end
            
            
            dist_truth = marginalsJT;
            % Error for various samplers
            % Node marginals
            error_ind = zeros(num_examples,1);
            error_hmc = zeros(num_examples,1);
            error_ana = zeros(num_examples,1);
            error_ana_rb = zeros(num_examples,1);
            error_ana_rb_lbp = zeros(num_examples,1);
            error_cmh = zeros(num_examples,1);
            error_cmh_lbp = zeros(num_examples,1);
            error_ana_gibbs = zeros(num_examples,1);
            error_ana_gibbs_rb = zeros(num_examples,1);
            error_ana_gibbs_rb_lbp = zeros(num_examples,1);
            error_sw = zeros(num_examples,1);
            
            % Pairwise node marginals
            err_pw_ind = zeros(num_examples,1);
            err_pw_hmc = zeros(num_examples,1);
            err_pw_ana = zeros(num_examples,1);
            err_pw_ana_rb = zeros(num_examples,1);
            err_pw_ana_rb_lbp = zeros(num_examples,1);
            err_pw_cmh = zeros(num_examples,1);
            err_pw_cmh_lbp = zeros(num_examples,1);
            err_pw_ana_gibbs = zeros(num_examples,1);
            err_pw_ana_gibbs_rb = zeros(num_examples,1);
            err_pw_ana_gibbs_rb_lbp = zeros(num_examples,1);
            err_pw_sw = zeros(num_examples,1);
            
            % Act of magnetization
            act_mag_mh_lbp = zeros(num_examples,1);
            act_mag_hmc = zeros(num_examples,1);
            act_mag_cmh = zeros(num_examples,1);
            act_mag_cmh_lbp = zeros(num_examples,1);
            act_mag_ana = zeros(num_examples,1);
            act_mag_ana_rb = zeros(num_examples,1);
            act_mag_ana_rb_lbp = zeros(num_examples,1);
            act_mag_ana_gibbs = zeros(num_examples,1);
            act_mag_ana_gibbs_rb = zeros(num_examples,1);
            act_mag_ana_gibbs_rb_lbp = zeros(num_examples,1);
            
            for q=1:num_examples
                q
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Exact-HMC
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if exact_hmc == 1
                    tic
                    disp('Exact_HMC')
                    t = 1.5;
                    if ising_2d ==1
                        t = 2.5;
                    end
                    T=t*pi;
                    [samples_hmc, loglik_hmc, energy_hmc] = HMC_binary(is1,T,number_samples, initial_point);
                    fn_evlas_hmc = number_samples*((is1.dim)*t + (is1.dim) + clique_size*(is1.dim)*(t-0.5));
                    toc
                    if node_marginals==1
                        error_hmc(q,1)= sqrt(mean( (dist_truth - emp_dist(samples_hmc(:,end))) .^2));
                    end
                    if pair_marg==1
                        err_pw_hmc(q,1) = sqrt(mean((sum(abs(empericalCounts(samples_hmc, a) - true_marg), 2)).^2));
                    end
                    if act_mag==1
                        act_mag_hmc(q,1) = acorrtime(sample_to_mag(samples_hmc)'); % pass as a column vector to acorrtime
                    end
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Independence sampler LBP
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if ind_sampler == 1
                    tic
                    disp('Independent sampler')
                    [samples_ind, log_lik_ind, nu_samples_ind] = independence_LBP(is1,  fn_evlas_hmc, clique_size, initial_point, dist_LBP);
                    r_ind = nu_samples_ind/number_samples;
                    toc
                    if node_marginals==1
                        error_ind(q,1)= sqrt(mean( (dist_truth - emp_dist(samples_ind(:,end))) .^2));
                    end
                    if pair_marg==1
                        err_pw_ind(q,1) = sqrt(mean((sum(abs(empericalCounts(samples_ind, a) - true_marg), 2)).^2));
                    end
                    if act_mag==1
                        act_mag_mh_lbp(q,1) = acorrtime(sample_to_mag(thin(samples_ind', 0, r_ind, nu_samples_ind-1)')'); % pass as a column vector to acorrtime
                    end
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Coordiante MH
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if cmh == 1
                    tic
                    disp('Simple Coordinate MH')
                    [samples_CMH, log_lik_CMH, nu_samples_cmh, max_CMH] = CMH(is1,  fn_evlas_hmc, clique_size, initial_point);
                    r_cmh = nu_samples_cmh/number_samples;
                    toc
                    if node_marginals==1
                        error_cmh(q,1) = sqrt(mean( (dist_truth - emp_dist(samples_CMH(:,end))) .^2));
                    end
                    if pair_marg==1
                        err_pw_cmh(q,1) = sqrt(mean((sum(abs(empericalCounts(samples_CMH, a) - true_marg), 2)).^2));
                    end
                    if act_mag==1
                        act_mag_cmh(q,1) = acorrtime(sample_to_mag(thin(samples_CMH', 0, r_cmh, nu_samples_cmh-1)')'); % pass as a column vector to acorrtime
                    end
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Coordiante MH with LBP implemented using DES
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if cmh_lbp == 1
                    tic
                    disp('Coordinate MH with LBP')
                    [samples_CMH_lbp, log_lik_CMH_lbp, nu_samples_cmh_lbp, max_CMH_LBP] = CMH_LBP_DES(is1,  fn_evlas_hmc, clique_size, initial_point, dist_LBP);
                    r_cmh_lbp = nu_samples_cmh_lbp/number_samples;
                    toc
                    if node_marginals==1
                        error_cmh_lbp(q,1) = sqrt(mean( (dist_truth - emp_dist(samples_CMH_lbp(:,end))) .^2));
                    end
                    if pair_marg==1
                        err_pw_cmh_lbp(q,1) = sqrt(mean((sum(abs(empericalCounts(samples_CMH_lbp, a) - true_marg), 2)).^2));
                    end
                    if act_mag==1
                        act_mag_cmh_lbp(q,1) = acorrtime(sample_to_mag(thin(samples_CMH_lbp', 0, r_cmh_lbp, nu_samples_cmh_lbp-1)')');
                    end
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Analytic Slice Sampling (with RB)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if ana_slice == 1
                    tic
                    disp('Analytic Slice Sampling')
                    [samples_ana, dist_ana, mag_ana, loglik_ana, nu_samples_ana, emp_count_ana, emp_counts_ana]= analytic_slice_new( is1, fn_evlas_hmc, clique_size, initial_point, a);%2*(dist_LBP>0.5)-1%analytic_slice_new( is1, fn_evlas_hmc, clique_size, initial_point);%
                    r_a = nu_samples_ana/number_samples;
                    toc
                    if node_marginals==1
                        error_ana(q,1)= sqrt(mean( (dist_truth - emp_dist(samples_ana(:,end))) .^2));
                        error_ana_rb(q,1) = sqrt(mean(  (dist_truth - (dist_ana(:,end)))  .^2));
                    end
                    if pair_marg==1
                        err_pw_ana(q,1) = sqrt(mean((sum(abs(empericalCounts(samples_ana, a) - true_marg), 2)).^2));
                        err_pw_ana_rb(q,1) = sqrt(mean((sum(abs(emp_count_ana - true_marg), 2)).^2));
                    end
                    if act_mag==1
                        act_mag_ana(q,1) = acorrtime(sample_to_mag(thin(samples_ana', 0, r_a, nu_samples_ana-1)')');
                        act_mag_ana_rb(q,1) = acorrtime(thin(mag_ana', 0, r_a, nu_samples_ana-1));
                    end
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Analytic Slice Sampling with RB + LBP
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if ana_slice_rb_lbp == 1
                    tic
                    disp('Analytic Slice Sampling with LBP')
                    [samples_ana_lbp, dist_ana_lbp, mag_ana_lbp, loglik_ana_lbp, nu_samples_ana_lbp, emp_count_ana_lbp, emp_counts_ana_lbp]= Stretched_analytic_slice_new( is1, fn_evlas_hmc, clique_size, info_on_off, initial_point, dist_LBP, a);%2*(dist_LBP>0.5)-1%analytic_slice_new( is1, fn_evlas_hmc, clique_size, initial_point);%
                    r_a_lbp = nu_samples_ana_lbp/number_samples;
                    toc
                    if node_marginals==1
                        error_ana_rb_lbp(q,1) = sqrt(mean(  (dist_truth - (dist_ana_lbp(:,end)))  .^2));
                    end
                    if pair_marg==1
                        err_pw_ana_rb_lbp(q,1) = sqrt(mean((sum(abs(emp_count_ana_lbp - true_marg),2)).^2));
                    end
                    if act_mag==1
                        act_mag_ana_rb_lbp(q,1) = acorrtime(thin(mag_ana_lbp', 0, r_a_lbp, nu_samples_ana_lbp-1));
                    end
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Analytic Gibbs Sampling (with RB)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if ana_gibbs == 1
                    tic
                    disp('Analytic Gibbs Sampling')
                    [samples_ana_gibbs, dist_ana_gibbs_rb, mag_ana_gibbs, loglik_ana_gibbs, nu_samples_ana_gibbs, emp_count_ana_gibbs, emp_counts_gibbs] = analytic_gibbs_new( is1, fn_evlas_hmc, clique_size, info_on_off, initial_point, a);
                    r_a_g = nu_samples_ana_gibbs/number_samples;
                    toc
                    if node_marginals==1
                        error_ana_gibbs = sqrt(mean( (dist_truth - emp_dist(samples_ana_gibbs(:,end))) .^2));
                        error_ana_gibbs_rb = sqrt(mean(  (dist_truth - (dist_ana_gibbs_rb(:,end)))  .^2));
                    end
                    if pair_marg==1
                        err_pw_ana_gibbs(q,1) = sqrt(mean((sum(abs(empericalCounts(samples_ana_gibbs, a) - true_marg), 2)).^2));
                        err_pw_ana_gibbs_rb(q,1) = sqrt(mean((sum(abs(emp_count_ana_gibbs - true_marg),2)).^2));
                    end
                    if act_mag==1
                        act_mag_ana_gibbs(q,1) = acorrtime(sample_to_mag(thin(samples_ana_gibbs', 0, r_a_g, nu_samples_ana_gibbs-1)')');
                        act_mag_ana_gibbs_rb(q,1) = acorrtime(thin(mag_ana_gibbs', 0, r_a_g, nu_samples_ana_gibbs-1));
                    end
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Analytic Gibbs Sampling with LBP
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if ana_gibbs_rb_lbp == 1
                    tic
                    disp('Analytic Gibbs Sampling with LBP')
                    [samples_ana_gibbs_lbp, dist_ana_gibbs_lbp, mag_ana_gibbs_lbp, loglik_ana_gibbs_lbp, nu_samples_ana_gibbs_lbp, emp_count_ana_gibbs_lbp, emp_counts_gibbs_lbp] = Stretched_analytic_gibbs_new( is1, fn_evlas_hmc, clique_size, info_on_off, initial_point, dist_LBP, a);
                    r_a_g_lbp = nu_samples_ana_gibbs_lbp/number_samples;
                    toc
                    if node_marginals==1
                        error_ana_gibbs_rb_lbp(q,1) = sqrt(mean(  (dist_truth - (dist_ana_gibbs_lbp(:,end)))  .^2));
                    end
                    if pair_marg==1
                        err_pw_ana_gibbs_rb_lbp(q,1) = sqrt(mean((sum(abs(emp_count_ana_gibbs_lbp - true_marg),2)).^2));
                    end
                    if act_mag==1
                        act_mag_ana_gibbs_rb_lbp(q,1) = acorrtime(thin(mag_ana_gibbs_lbp', 0, r_a_g_lbp, nu_samples_ana_gibbs_lbp-1));
                    end
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Swendsen Wang
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if sw == 1
                    tic
                    disp('Swendsen Wang')
                    [samples_sw, dist_sw, loglik_sw, nu_samples_sw] = SW_1D( is1, fn_evlas_hmc, initial_point);
                    r_sw = nu_samples_sw/number_samples;
                    toc
                    if node_marginals==1
                        error_sw(q,1) = sqrt(mean( (dist_truth - emp_dist(samples_sw(:,end))) .^2));
                    end
                    if pair_marg==1
                        err_pw_sw(q,1) = sqrt(mean((sum(abs(empericalCounts(samples_sw, a) - true_marg), 2)).^2));
                    end
                end
                
            end
            
            % Node marginals
            if  plot_marginals==1
                mean_err = [sum(abs(dist_LBP-dist_truth)), mean(error_ind), mean(error_hmc), mean(error_cmh), mean(error_cmh_lbp), mean(error_ana), mean(error_ana_rb), mean(error_ana_rb_lbp),  ...
                    mean(error_ana_gibbs), mean(error_ana_gibbs_rb), mean(error_ana_gibbs_rb_lbp)];
                std_err =  [0, std(error_ind), std(error_hmc), std(error_cmh), std(error_cmh_lbp), std(error_ana), std(error_ana_rb), std(error_ana_rb_lbp),  ...
                    std(error_ana_gibbs), std(error_ana_gibbs_rb), std(error_ana_gibbs_rb_lbp)];
                plot_marg(mean_err, std_err, 'Node marginals', 0);
                name = strcat('Temp', num2str(temp), 'Bias', num2str(scale_bias), '.fig');
                path = '/Users/Jalaj/Documents/Github - LBPSS/Outputs_after_NIPS/Node_marginals';
                savefig(gcf, fullfile(path, name))
            end
            
            % Pairwise node marginals
            if  plot_marg_pairwise ==1
                mean_err_pw = [mean(error_ind), mean(err_pw_hmc), mean(err_pw_cmh), mean(err_pw_cmh_lbp), mean(err_pw_ana), mean(err_pw_ana_rb), ...
                    mean(err_pw_ana_rb_lbp),  mean(err_pw_ana_gibbs), mean(err_pw_ana_gibbs_rb), mean(err_pw_ana_gibbs_rb_lbp)];
                std_err_pw = [std(error_ind), std(err_pw_hmc), std(err_pw_cmh), std(err_pw_cmh_lbp), std(err_pw_ana), std(err_pw_ana_rb), ...
                    std(err_pw_ana_rb_lbp), std(err_pw_ana_gibbs), std(err_pw_ana_gibbs_rb), std(err_pw_ana_gibbs_rb_lbp)];
                plot_fn(mean_err_pw, std_err_pw,'Pairwise marginals', 0);
                name = strcat('Temp', num2str(temp), 'Bias', num2str(scale_bias), '.fig');
                path = '/Users/Jalaj/Documents/Github - LBPSS/Outputs_after_NIPS/Pairwise_marginals';
                savefig(gcf, fullfile(path, name))
                
            end
            
            % Plotting act for magnetization
            if  plot_act_mag ==1
                mean_act= [mean(act_mag_mh_lbp), mean(act_mag_hmc), mean(act_mag_cmh), mean(act_mag_cmh_lbp), mean(act_mag_ana), mean(act_mag_ana_rb), ...
                    mean(act_mag_ana_rb_lbp), mean(act_mag_ana_gibbs), mean(act_mag_ana_gibbs_rb), mean(act_mag_ana_gibbs_rb_lbp)];
                std_act= [std(act_mag_mh_lbp), std(act_mag_hmc), std(act_mag_cmh), std(act_mag_cmh_lbp), std(act_mag_ana), std(act_mag_ana_rb), ...
                    std(act_mag_ana_rb_lbp), std(act_mag_ana_gibbs), std(act_mag_ana_gibbs_rb), std(act_mag_ana_gibbs_rb_lbp)];
                plot_fn(mean_act, std_act,'Auto-correlation time: magnetization', 1);
                name = strcat('Temp', num2str(temp), 'Bias', num2str(scale_bias), '.fig');
                path = '/Users/Jalaj/Documents/Github - LBPSS/Outputs_after_NIPS/Act_mag';
                savefig(gcf, fullfile(path, name))
            end
            
            % Plotting log-likes of the samples
            if plot_log_liks==1
                figure
                Z=200;
                plot(loglik_hmc(1:Z),'r')
                hold on
                plot(log_lik_CMH(1:Z),'k')
                hold on
                plot(log_lik_CMH_lbp(1:Z),'g')
                hold on
                plot(loglik_ana(1:Z),'b')
                hold on
                plot(loglik_ana_lbp(1:Z),'y')
                hold on
                plot(loglik_ana_gibbs(1:Z),'c')
                hold on
                plot(loglik_ana_gibbs_lbp(1:Z),'m')
                hold on
                %             plot(loglik_sw(1:Z), 'g');
                %             hold on
                h=legend('HMC', 'CMH', 'CMH LBP', 'AAS', 'RBS LBP', 'AAG','RBG LBP');
                set(h);
                xlabel('Iterations');
                ylabel('Log-likelihood');
            end
                   
        end
        
    end
    
end


%% Old PLotting fucntions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RMSE for node marginals with iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% N=25; % Check error after every N equivalent iterations
% error_hmc = zeros(num_examples, number_samples/N);
% error_ana = zeros(num_examples, number_samples/N);
% error_ana_lbp = zeros(num_examples, number_samples/N);
% error_cmh = zeros(num_examples, number_samples/N);
% error_cmh_lbp = zeros(num_examples, number_samples/N);
% error_ana_gibbs = zeros(num_examples, number_samples/N);
% error_ana_gibbs_rb = zeros(num_examples, number_samples/N);
% error_ana_gibbs_rb_lbp = zeros(num_examples, number_samples/N);
% error_sw = zeros(num_examples, number_samples/N);

% if node_marginals==1
%     init = sqrt(mean( (dist_truth - emp_dist(initial_point))  .^2));
%     error_hmc(q,1)=init;
%     error_ana(q,1)=init;
%     error_ana_lbp(q,1)=init;
%     error_cmh(q,1)=init;
%     error_cmh_lbp(q,1)=init;
%     error_ana_gibbs=init;
%     error_ana_gibbs_rb=init;
%     error_ana_gibbs_rb_lbp(q,1)=init;
%     error_sw(q,1)=init;
%
%     for j=2:(number_samples/N)
%         error_hmc(q,j) = sqrt(mean(   (dist_truth - emp_dist(samples_hmc(:,2:j*N-1)))  .^2));
%         %                 error_ana(q,j) = sqrt(mean(   (dist_truth - mean(cat(1,dist_ana(:,2:round(j*N*r_a-1))),2))   .^2));
%         error_ana_lbp(q,j) = sqrt(mean(   (dist_truth - mean(cat(1,dist_ana_lbp(:,2:round(j*N*r_a_lbp-1))),2))   .^2));
%         %                 error_cmh(q,j) = sqrt(mean(  (dist_truth - emp_dist(samples_CMH(:,2:round(j*N*r_cmh-1))))  .^2));
%         error_cmh_lbp(q,j) = sqrt(mean(  (dist_truth - emp_dist(samples_CMH_lbp(:,2:round(j*N*r_cmh_lbp-1))))  .^2));
%         %                 error_ana_gibbs(q,j) = sqrt(mean(  (dist_truth - mean(cat(1,dist_ana_gibbs(:,2:round(j*N*r_a_g-1))),2))  .^2));
%         %                 error_ana_gibbs_rb(q,j) = sqrt(mean(  (dist_truth - mean(cat(1,dist_ana_gibbs_rb(:,2:round(j*N*r_a_g-1))),2))  .^2));
%         error_ana_gibbs_rb_lbp(q,j) = sqrt(mean(  (dist_truth - mean(cat(1,dist_ana_gibbs_lbp(:,2:round(j*N*r_a_g_lbp-1))),2))  .^2));
%         error_sw(q,j) = sqrt(mean(  (dist_truth - emp_dist(samples_sw(:,2:round(j*N*r_sw-1))))  .^2));
%     end
%
% end

% if  plot_marginals==1
%     figure
%     semilogy(0:length(error_ana_lbp)-1, mean(error_ana_lbp,1), '--', 'Color', [152,78,163]/255);
%     hold on
%     semilogy(0:length(error_ana_gibbs_rb_lbp)-1, mean(error_ana_gibbs_rb_lbp,1),'--', 'Color', [77,175,74]/255);
%     hold on
%     semilogy(0:length(error_hmc)-1, mean(error_hmc,1), '--', 'Color', [228,26,28]/255);
%     hold on
%     semilogy(0:length(error_cmh_lbp)-1, mean(error_cmh_lbp,1), '--', 'Color', [55,126,184]/255);
%     hold on
%     semilogy(0:length(error_sw)-1, mean(error_sw,1), '--', 'Color', [55,126,184]/255);
%     hold on
%     h=legend('RBS LBP', 'RBG LBP','HMC', 'CMH LBP', 'SW');
%     set(h);
%     xlabel('Function Evaluations');
%     ylabel('RMSE (Marginals)');
%     str=sprintf('Bias scale = %d', scale);
%     title(str);
%     set(gcf,'units','points','position',[10,10,800,800]);
%     name = strcat('RMSE Temp', num2str(temp), 'Bias', num2str(scale), '.fig');
%     path = '/Users/Jalaj/Documents/Github - LBPSS/New Outputs';
%     savefig(gcf, fullfile(path, name))
%
% end



% Node Marginal Errors
%         if  plot_marginals==1
%             figure
%             mean_err = [mean(error_ind), mean(error_hmc), mean(error_cmh), mean(error_cmh_lbp), mean(error_ana), mean(error_ana_rb), mean(error_ana_rb_lbp),  ....
%                 mean(error_ana_gibbs), mean(error_ana_gibbs_rb), mean(error_ana_gibbs_rb_lbp)];
%             std_err =  [std(error_ind), std(error_hmc), std(error_cmh), std(error_cmh_lbp), std(error_ana), std(error_ana_rb), std(error_ana_rb_lbp),  ....
%                 std(error_ana_gibbs), std(error_ana_gibbs_rb), std(error_ana_gibbs_rb_lbp)];
%             hold on
%             for k = 1:length(mean_err)
%                 e1 = errorbar(k,mean_err(k),std_err(k),'o');
%                 set(e1,'Color', [228,26,28]/255)
%                 set(e1,'MarkerEdgeColor',[228,26,28]/255)
%             end
%             hold off
%             %             str=sprintf('Node marginals:Bias scale = %d', scale);
%             %             title(str)
%             title('Node marginals');
%             ylabel('RMSE')
%             box on
%             % Change the labels for the tick marks on the x-axis
%             xtl1 = '\begin{tabular}{c} MH-LBP \end{tabular}';
%             xtl2 = '\begin{tabular}{c} HMC \end{tabular}';
%             xtl3 = '\begin{tabular}{c} CMH \end{tabular}';
%             xtl4 = '\begin{tabular}{c} CMH \\+\\LBP\end{tabular}';
%             xtl5 = '\begin{tabular}{c} AAS \end{tabular}';
%             xtl6 = '\begin{tabular}{c} AAS \\+\\RB\end{tabular}';
%             xtl7 = '\begin{tabular}{c} AAS \\ + \\ RB\\ + \\ LBP\end{tabular}';
%             xtl8 = '\begin{tabular}{c} AAG \end{tabular}';
%             xtl9 = '\begin{tabular}{c} AAG\\+\\ RB \end{tabular}';
%             xtl10 = '\begin{tabular}{c} AAG \\ + \\ RB\\ + \\ LBP\end{tabular}';
%             Algorithms = {xtl1, xtl2, xtl3, xtl4, xtl5, xtl6, xtl7, xtl8, xtl9, xtl10};
%             set(gca, 'XTick', 1:10, 'XTickLabel', Algorithms, 'TickLabelInterpreter', 'latex');
%             name = strcat('Temp', num2str(temp), 'Bias', num2str(scale), '.fig');
%             path = '/Users/Jalaj/Documents/Github - LBPSS/Outputs_after_NIPS/Pairwise_marginals';
%             savefig(gcf, fullfile(path, name))
%             %             set(gcf,'units','points','position',[10,10,800,800]);
%
%         end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ACtime for node-marginal estimates
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% if act_marginals==1
%     %Getting the baseline: Ana_gibbs_LBP
%     thin_dist_ana_gibbs_lbp = thin(dist_ana_gibbs_lbp', 0, r_a_g_lbp, nu_samples_ana_gibbs_lbp-1);
%     corr_time_dist_ana_gibbs_lbp = acorrtime(thin_dist_ana_gibbs_lbp);
%     
%     %HMC
%     corr_time_dist_hmc = acorrtime(samples_hmc');
%     err_hmc(q,1) = sum(corr_time_dist_hmc-corr_time_dist_ana_gibbs_lbp);
%     
%     %AAS
%     thin_dist_ana = thin(dist_ana', 0, r_a, nu_samples_ana-1);
%     corr_time_dist_ana = acorrtime(thin_dist_ana);
%     err_ana(q,1) = sum(corr_time_dist_ana-corr_time_dist_ana_gibbs_lbp);
%     
%     %RBS+LBP
%     thin_dist_ana_lbp = thin(dist_ana_lbp', 0, r_a_lbp, nu_samples_ana_lbp-1);
%     corr_time_dist_ana_lbp = acorrtime(thin_dist_ana_lbp);
%     err_ana_lbp(q,1) = sum(corr_time_dist_ana_lbp-corr_time_dist_ana_gibbs_lbp);
%     
%     %CMH
%     thin_dist_cmh = thin(samples_CMH', 0, r_cmh, nu_samples_cmh-1);
%     corr_time_dist_cmh = acorrtime(thin_dist_cmh);
%     err_cmh(q,1) = sum(corr_time_dist_cmh-corr_time_dist_ana_gibbs_lbp);
%     
%     %CMH+LBP
%     thin_dist_cmh_lbp = thin(samples_CMH_lbp', 0, r_cmh_lbp, nu_samples_cmh_lbp-1);
%     corr_time_dist_cmh_lbp = acorrtime(thin_dist_cmh_lbp);
%     err_cmh_lbp(q,1) = sum(corr_time_dist_cmh_lbp-corr_time_dist_ana_gibbs_lbp);
%     
%     %AAG
%     thin_dist_ana_gibbs = thin(samples_ana_gibbs', 0, r_a_g, nu_samples_ana_gibbs-1);
%     corr_time_dist_ana_gibbs = acorrtime(thin_dist_ana_gibbs);
%     err_ana_gibbs(q,1) = sum(corr_time_dist_ana_gibbs-corr_time_dist_ana_gibbs_lbp);
%     
%     %RBG
%     thin_dist_ana_gibbs_rb = thin(dist_ana_gibbs_rb', 0, r_a_g, nu_samples_ana_gibbs-1);
%     corr_time_dist_ana_gibbs_rb = acorrtime(thin_dist_ana_gibbs_rb);
%     err_ana_gibbs_rb(q,1) = sum(corr_time_dist_ana_gibbs_rb-corr_time_dist_ana_gibbs_lbp);
% end