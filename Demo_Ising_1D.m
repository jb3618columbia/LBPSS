% Implementation of Different Algorithms for 1D Ising Models
% Parameters
d=50;
temp_vec=[-10*pi];
scale_vec = [5];

% This computes and plots the errors in node marginals
node_marginals = 0;
plot_marginals = 0;  

% This computes and plots the errors in pariwise marginals
pair_marg=0;
plot_pair_marg=0;  

% This plots the log-likelihood of samples 
plot_log_liks=0;  

% This computes and plots the auto-correlation time for node marginals
act_marginals=1;
plot_act_marginals=1;


for u=1:1:length(temp_vec)
    temp = temp_vec(u);
    
    for v=1:1:length(scale_vec);
        scale = scale_vec(v);
        is1 = Ising1D_new(d,temp, scale);  % Create 1D Ising Object
        clique_size=2; %Clique size
        number_samples = 2000;
        num_examples = 1;
        initial_point = sign(normrnd(0,1,d,1));
        
        % Algorithms:
        truth = 0;
        ana_slice = 1;  % plan to only show slice with rb and lbp
        ana_gibbs = 1;
        ana_gibbs_rb = 1;
        ana_slice_rb_lbp = 1;
        ana_gibbs_rb_lbp = 1; % in gibbs sampling we show rb, rb+lbp
        exact_hmc = 1; 
        cmh = 1;
        cmh_lbp = 1;
        sw = 1;
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
           
        [edgeStruct, dist_LBP, edgeBelLBP, Z_LBP, marginalsJT, edgeBelJT, Z_JT] = get_LBP_marginals( -is1.bias , -is1.M);
        a = [6,7]; % (i,j) for getting the pairwise marginals
        true_marg = reshape(edgeBelJT(:,:,find(ismember(edgeStruct.edgeEnds, a ,'rows'))),1,4);
        % Reshape operation here gives the answer for [(1,1), (-1,1), (1,-1), (-1,-1)]
        
        % Error for various samplers
        N=25; % Check error after every N equivalent iterations
        error_hmc = zeros(num_examples, number_samples/N);
        error_ana = zeros(num_examples, number_samples/N);
        error_ana_lbp = zeros(num_examples, number_samples/N);
        error_cmh = zeros(num_examples, number_samples/N);
        error_cmh_lbp = zeros(num_examples, number_samples/N);
        error_ana_gibbs = zeros(num_examples, number_samples/N);
        error_ana_gibbs_rb = zeros(num_examples, number_samples/N);
        error_ana_gibbs_rb_lbp = zeros(num_examples, number_samples/N);
        error_sw = zeros(num_examples, number_samples/N);   
        
        err_pw_hmc = zeros(num_examples,1);
        err_pw_ana = zeros(num_examples,1);
        err_pw_ana_lbp = zeros(num_examples,1);
        err_pw_cmh = zeros(num_examples,1);
        err_pw_cmh_lbp = zeros(num_examples,1);
        err_pw_ana_gibbs = zeros(num_examples,1);
        err_pw_ana_gibbs_rb = zeros(num_examples,1);
        err_pw_ana_gibbs_rb_lbp = zeros(num_examples,1);   
        err_pw_sw = zeros(num_examples,1);   
        
        for q=1:num_examples
            q
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Exact-HMC
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if exact_hmc == 1
                tic
                disp('Exact_HMC')
                t = 1.5; T=t*pi;
                [samples_hmc, loglik_hmc, energy_hmc] = HMC_binary(is1,T,number_samples, initial_point);
                fn_evlas_hmc = number_samples*((is1.dim)*t + (is1.dim) + clique_size*(is1.dim));
                % Doing the thinning step 
                corr_time_dist_hmc = acorrtime(samples_hmc');
                toc
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Analytic Slice Sampling
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ana_slice == 1
                tic
                disp('Analytic Slice Sampling')
                [samples_ana, dist_ana, loglik_ana, nu_samples_ana]= analytic_slice_new( is1, fn_evlas_hmc, clique_size, initial_point);%2*(dist_LBP>0.5)-1%analytic_slice_new( is1, fn_evlas_hmc, clique_size, initial_point);%
                r_a = nu_samples_ana/number_samples;
                % Doing the thinning step 
                thin_dist_ana = thin(dist_ana', 0, r_a, nu_samples_ana-1);
                corr_time_dist_ana = acorrtime(thin_dist_ana);
                toc
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Analytic Slice Sampling with LBP
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ana_slice_rb_lbp == 1
                tic
                disp('Analytic Slice Sampling with LBP')
                [samples_ana_lbp, dist_ana_lbp, loglik_ana_lbp, nu_samples_ana_lbp, emp_count_ana_lbp]= Stretched_analytic_slice_new( is1, fn_evlas_hmc, clique_size, info_on_off, initial_point, dist_LBP, a);%2*(dist_LBP>0.5)-1%analytic_slice_new( is1, fn_evlas_hmc, clique_size, initial_point);%
                r_a_lbp = nu_samples_ana_lbp/number_samples;
                % Doing the thinning step 
                thin_dist_ana_lbp = thin(dist_ana_lbp', 0, r_a_lbp, nu_samples_ana_lbp-1);
                corr_time_dist_ana_lbp = acorrtime(thin_dist_ana_lbp);
                toc
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Coordiante MH 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if cmh == 1
                tic
                disp('Simple Coordinate MH')
                [samples_CMH, log_lik_CMH, nu_samples_cmh] = CMH(is1,  fn_evlas_hmc, clique_size, initial_point);
                r_cmh = nu_samples_cmh/number_samples;
                % Doing the thinning step, now on the samples
                thin_dist_cmh = thin(samples_CMH', 0, r_cmh, nu_samples_cmh-1);
                corr_time_dist_cmh = acorrtime(thin_dist_cmh);
                toc
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Coordiante MH with LBP implemented using DES 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if cmh_lbp == 1
                tic
                disp('Coordinate MH with LBP')
                [samples_CMH_lbp, log_lik_CMH_lbp, nu_samples_cmh_lbp] = CMH_LBP_DES(is1,  fn_evlas_hmc, clique_size, initial_point, dist_LBP);
                r_cmh_lbp = nu_samples_cmh_lbp/number_samples;
                toc
                % Doing the thinning step, now on the samples
                thin_dist_cmh_lbp = thin(samples_CMH_lbp', 0, r_cmh_lbp, nu_samples_cmh_lbp-1);
                corr_time_dist_cmh_lbp = acorrtime(thin_dist_cmh_lbp);
            end
           
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Analytic Gibbs Sampling
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ana_gibbs == 1
               tic
               disp('Analytic Gibbs Sampling')
               [samples_ana_gibbs, dist_ana_gibbs_rb, loglik_ana_gibbs, nu_samples_ana_gibbs, emp_count_ana_gibbs] = analytic_gibbs_new( is1, fn_evlas_hmc, clique_size, info_on_off, initial_point, a);
               r_a_g = nu_samples_ana_gibbs/number_samples;
               % To get the estimates without rao-blackwellization
               dist_ana_gibbs = emp_dist(samples_ana_gibbs);
               
               % Doing the thinning step for rao-blackwellized version
                thin_dist_ana_gibbs = thin(dist_ana_gibbs_rb', 0, r_a_g, nu_samples_ana_gibbs-1);
                corr_time_dist_ana_gibbs = acorrtime(thin_dist_ana_gibbs);
               toc
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Analytic Gibbs Sampling with LBP
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ana_gibbs_rb_lbp == 1
                tic
                disp('Analytic Gibbs Sampling with LBP')
                [samples_ana_gibbs_lbp, dist_ana_gibbs_lbp, loglik_ana_gibbs_lbp, nu_samples_ana_gibbs_lbp, emp_count_ana_gibbs_lbp] = Stretched_analytic_gibbs_new( is1, fn_evlas_hmc, clique_size, info_on_off, initial_point, dist_LBP, a);
                r_a_g_lbp = nu_samples_ana_gibbs_lbp/number_samples;
                % Doing the thinning step 
                thin_dist_ana_gibbs_lbp = thin(dist_ana_gibbs_lbp', 0, r_a_g_lbp, nu_samples_ana_gibbs_lbp-1);
                corr_time_dist_ana_gibbs_lbp = acorrtime(thin_dist_ana_gibbs_lbp);
                toc
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
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % RMSE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if node_marginals==1
                init = sqrt(mean( (dist_truth - emp_dist(initial_point))  .^2));
                error_hmc(q,1)=init;
                error_ana(q,1)=init;
                error_ana_lbp(q,1)=init;
                error_cmh(q,1)=init;
                error_cmh_lbp(q,1)=init;
                error_ana_gibbs=init;
                error_ana_gibbs_rb=init;
                error_ana_gibbs_rb_lbp(q,1)=init;
                error_sw(q,1)=init;

                for j=2:(number_samples/N)
                    error_hmc(q,j) = sqrt(mean(   (dist_truth - emp_dist(samples_hmc(:,2:j*N-1)))  .^2));
    %                 error_ana(q,j) = sqrt(mean(   (dist_truth - mean(cat(1,dist_ana(:,2:round(j*N*r_a-1))),2))   .^2));
                    error_ana_lbp(q,j) = sqrt(mean(   (dist_truth - mean(cat(1,dist_ana_lbp(:,2:round(j*N*r_a_lbp-1))),2))   .^2));
    %                 error_cmh(q,j) = sqrt(mean(  (dist_truth - emp_dist(samples_CMH(:,2:round(j*N*r_cmh-1))))  .^2));
                    error_cmh_lbp(q,j) = sqrt(mean(  (dist_truth - emp_dist(samples_CMH_lbp(:,2:round(j*N*r_cmh_lbp-1))))  .^2));
    %                 error_ana_gibbs(q,j) = sqrt(mean(  (dist_truth - mean(cat(1,dist_ana_gibbs(:,2:round(j*N*r_a_g-1))),2))  .^2));
    %                 error_ana_gibbs_rb(q,j) = sqrt(mean(  (dist_truth - mean(cat(1,dist_ana_gibbs_rb(:,2:round(j*N*r_a_g-1))),2))  .^2));
                    error_ana_gibbs_rb_lbp(q,j) = sqrt(mean(  (dist_truth - mean(cat(1,dist_ana_gibbs_lbp(:,2:round(j*N*r_a_g_lbp-1))),2))  .^2));
                    error_sw(q,j) = sqrt(mean(  (dist_truth - emp_dist(samples_sw(:,2:round(j*N*r_sw-1))))  .^2));
                end
                
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Pairwise marginal error
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if pair_marg==1
                err_pw_hmc(q,1) = sum(abs(empericalCounts(samples_hmc, a)' - true_marg));
                err_pw_ana(q,1) = sum(abs(empericalCounts(samples_ana, a)' - true_marg));
                err_pw_ana_lbp(q,1) = sum(abs(emp_count_ana_lbp - true_marg));
                err_pw_cmh(q,1) = sum(abs(empericalCounts(samples_CMH, a)' - true_marg));
                err_pw_cmh_lbp(q,1) = sum(abs(empericalCounts(samples_CMH_lbp, a)' - true_marg));
                err_pw_ana_gibbs(q,1) = sum(abs(empericalCounts(samples_ana_gibbs, a)' - true_marg));
                err_pw_ana_gibbs_rb(q,1) = sum(abs(emp_count_ana_gibbs - true_marg));
                err_pw_ana_gibbs_rb_lbp(q,1) = sum(abs(emp_count_ana_gibbs_lbp - true_marg));
                err_pw_sw(q,1) = sum(abs(empericalCounts(samples_sw, a)' - true_marg));
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ACtime for node-marginal estimates
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if act_marginals==1
                err_hmc = sum(corr_time_dist_hmc-corr_time_dist_ana_gibbs_lbp);
                err_ana = sum(corr_time_dist_hmc-corr_time_dist_ana_gibbs_lbp);
                err_ana_lbp = sum(corr_time_dist_ana_lbp-corr_time_dist_ana_gibbs_lbp);
                err_cmh = sum(corr_time_dist_cmh-corr_time_dist_ana_gibbs_lbp);
                err_cmh_lbp = sum(corr_time_dist_cmh_lbp-corr_time_dist_ana_gibbs_lbp);
                err_ana_gibbs = sum(corr_time_dist_ana_gibbs-corr_time_dist_ana_gibbs_lbp);
            end
            
        end
         
        if plot_log_liks==1
            figure
            plot(loglik_ana_gibbs_lbp(1:50),'b')
            hold on
            plot(log_lik_CMH_lbp(1:50),'g')
            hold on
            plot(loglik_hmc(1:50),'r')
            h=legend('AAG LBP', 'CMH LBP','HMC');
            set(h);
            xlabel('Iterations');
            ylabel('Log-likelihood');
        end
        
%         figure
%         semilogy(loglik_ana_gibbs_lbp(1:100),'b')
%         hold on
%         semilogy(log_lik_CMH_lbp(1:100),'g')
%         hold on
%         semilogy(loglik_hmc(1:100),'r')

%         hold on
%         semilogy(log_lik_CMH_lbp(1:end),'b')
%         hold on 
%         semilogy(loglik_sw(1:end),'k')
%         str2=sprintf('Logp:temp = %d, scale = %d', temp, scale);
%         title(str2)
%         saveas(gcf, str2, 'png')

        % Node Marginal Errors
        if  plot_marginals==1
            figure
            semilogy(0:length(error_ana_lbp)-1, mean(error_ana_lbp,1), '--', 'Color', [152,78,163]/255);
            hold on
            semilogy(0:length(error_ana_gibbs_rb_lbp)-1, mean(error_ana_gibbs_rb_lbp,1),'--', 'Color', [77,175,74]/255);
            hold on
            semilogy(0:length(error_hmc)-1, mean(error_hmc,1), '--', 'Color', [228,26,28]/255);
            hold on
            semilogy(0:length(error_cmh_lbp)-1, mean(error_cmh_lbp,1), '--', 'Color', [55,126,184]/255);
            hold on
            semilogy(0:length(error_sw)-1, mean(error_sw,1), '--', 'Color', [55,126,184]/255);
            hold on
            h=legend('RBS LBP', 'RBG LBP','HMC', 'CMH LBP', 'SW');
            set(h);
            xlabel('Function Evaluations');
            ylabel('RMSE (Marginals)');
            str=sprintf('Bias scale = %d', scale);
            title(str);
            set(gcf,'units','points','position',[10,10,800,800]);
            name = strcat('RMSE Temp', num2str(temp), 'Bias', num2str(scale), '.fig');
            path = '/Users/Jalaj/Documents/Github - LBPSS/New Outputs';
            savefig(gcf, fullfile(path, name))
         
        end
        
        % Pairwise node marginals
        if  plot_pair_marg ==1
            mean_err_pw = [mean(err_pw_hmc), mean(err_pw_ana), mean(err_pw_ana_lbp), mean(err_pw_cmh), mean(err_pw_cmh_lbp), mean(err_pw_ana_gibbs), mean(err_pw_ana_gibbs_rb), mean(err_pw_ana_gibbs_rb_lbp)];
            std_err_pw = [std(err_pw_hmc), std(err_pw_ana), std(err_pw_ana_lbp), std(err_pw_cmh), std(err_pw_cmh_lbp), std(err_pw_ana_gibbs), std(err_pw_ana_gibbs_rb), std(err_pw_ana_gibbs_rb_lbp)];

            h = bar(diag(mean_err_pw),'stacked');
            grid on
            l = cell(1,8);
            l{1}='HMC'; l{2}='AAS'; l{3}='AAS+LBP'; l{4}='CMH'; l{5}='CMH+LBP'; l{6}='AAG'; l{7}='AAG+RB'; l{8}='AAG+RB+LBP';    
            legend(h,l);
%             errorbar(1:9,mean_err_pw,std_err_pw,'.')
                        
        end
        
              
    end
    
end
            
        
        