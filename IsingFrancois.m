% Implementation of Different Algorithms for Ising Models
% Models:
% 1) 1D Ising Models
% 2) 2D Ising Models

% clear
% Parameters
d=10;
temp_vec=[-1*pi];
scale_vec = [20];

for u=1:1:length(temp_vec)
    
    temp = temp_vec(u);
    
    for v=1:1:length(scale_vec);
        
        scale = scale_vec(v);
        is1 = Ising1D_new(d,temp, scale);  % Create 1D Ising Object
        % is1 = Ising1D_rand_weight(d,temp);
        clique_size=2; %Clique size
        number_samples = 4000;
        num_examples = 1;
        initial_point = sign(normrnd(0,1,d,1));
        % initial_point = ones(d,1);
        
        % Algorithms:
        truth = true;
        ana_on_off = false;
        ana_gibbs_on_off =  false;
        exact_hmc_on_off = false;
        CMH_on_off = true;
        info_on_off = true;
        
        
        ground_truth = false;
        %%%%%%%%%%%%%%%%%%%%%%%%
        % Approx. ground truth
        %%%%%%%%%%%%%%%%%%%%%%%%
        %         if ground_truth == 1
        %             display('Getting samples to approximate grund truth')
        %             num_samp_truth = 100000;
        %             fn_evlas = clique_size*d*num_samp_truth;
        %             [samples_true, ~, ~] = CMH(is1, fn_evlas, clique_size, initial_point);
        %             dist_truth = emp_dist(samples_true);
        %
        %             [samples_true_1, ~, ~, ~]= analytic_slice_new( is1, fn_evlas_hmc, clique_size, initial_point);
        %             dist_truth_1 = emp_dist(samples_true_1);
        %             %     t = 1.5; T=t*pi;
        %             %     num_samp_truth = 50000;
        %             %     [samples_true, loglik_true, energy_true] = HMC_binary(is1,T,num_samp_truth, initial_point);
        %
        %         end
        %         dist_truth = 0.5*ones(d,1);
        
        
        
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
            % vector of probabilities of 2^d points, used for tv also 
            
            dist_truth = zeros(d,1);
            for j=0:2^d-1
                if mod(j,10000)==0
                end
                yy = (2*de2bi(j,d) - 1)';
                dist_truth = dist_truth + weight(1,j+1)*(yy >0);
                
            end
            
        end
        
        dist_LBP = get_LBP_marginals( -is1.bias , -is1.M);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Comparing Samplers
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        N=200; % Check error after every N equivalent iterations
        error_ana = zeros(num_examples, number_samples/N -1);
        error_ana_dist = zeros(num_examples, number_samples/N -1);
        error_ana_gibbs = zeros(num_examples, number_samples/N -1);
        error_ana_gibbs_dist = zeros(num_examples, number_samples/N -1);
        error_uss = zeros(num_examples, number_samples/N-1);
        error_uss_inout = zeros(num_examples, number_samples/N-1);
        error_hmc = zeros(num_examples, number_samples/N -1);
        error_cmh = zeros(num_examples, number_samples/N -1);
        
        error_ana_1 = zeros(num_examples, number_samples/N -1);
        error_ana_dist_1 = zeros(num_examples, number_samples/N -1);
        error_ana_gibbs_1 = zeros(num_examples, number_samples/N -1);
        error_ana_gibbs_dist_1 = zeros(num_examples, number_samples/N -1);
        error_uss_1 = zeros(num_examples, number_samples/N-1);
        error_uss_inout_1 = zeros(num_examples, number_samples/N-1);
        error_hmc_1 = zeros(num_examples, number_samples/N -1);
        error_cmh_1 = zeros(num_examples, number_samples/N -1);
        
%         error_ana_2 = zeros(num_examples, 1);
%         error_ana_dist_2 = zeros(num_examples, 1);
%         error_ana_gibbs_2 = zeros(num_examples, 1);
%         error_ana_gibbs_dist_2 = zeros(num_examples, 1);
%         error_hmc_2 = zeros(num_examples, 1);
%         error_cmh_2 = zeros(num_examples, 1);
        
        
        
        for q=1:num_examples
            
            q
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Metropolis with random flip proposal
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if CMH_on_off == 1
                tic
                disp('Coordinate MH')
                [samples_CMH, log_lik_CMH, nu_samples_cmh] = CMH(is1, fn_evlas_hmc, clique_size, initial_point);
                %[samples_CMH, log_lik_CMH, nu_samples_cmh] = CMH_LBP_DES(is1,  fn_evlas_hmc, clique_size, initial_point, dist_LBP);
                %[samples_CMH, log_lik_CMH, nu_samples_cmh, dist_cmh] = CMH_LBP_RB(is1,  fn_evlas_hmc, clique_size, initial_point, dist_LBP);
                %[samples_CMH, log_lik_CMH, nu_samples_cmh] = CMH_LBP(is1,  fn_evlas_hmc, clique_size, initial_point, dist_LBP);
                r_cmh = nu_samples_cmh/number_samples;
                macmh = mean(samples_CMH,1);
                toc
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Exact-HMC
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if exact_hmc_on_off == 1
                tic
                disp('Exact_HMC')
                t = 1.5; T=t*pi;
                [samples_hmc, loglik_hmc, energy_hmc] = HMC_binary(is1,T,number_samples, initial_point);
                mah = mean(samples_hmc,1);
                fn_evlas_hmc = number_samples*((is1.dim)*t + (is1.dim) + clique_size*(is1.dim));
                toc
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Analytic Slice Sampling
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if ana_on_off == 1
                tic
                disp('Analytic Slice Sampling')
                % [samples_ana, dist_ana, loglik_ana, fn_evals_ana, nu_samples_ana]= ussSampler(is1, 0, 1, fn_evlas_hmc, clique_size, initial_point);
                % More efficient way of doing slice sampling on a circle
                [samples_ana, dist_ana, loglik_ana, nu_samples_ana]= Stretched_analytic_slice_new( is1, fn_evlas_hmc, clique_size, info_on_off, initial_point, dist_LBP);%2*(dist_LBP>0.5)-1%analytic_slice_new( is1, fn_evlas_hmc, clique_size, initial_point);%
                
                r_a = nu_samples_ana/number_samples;
                mauss_ana = mean(samples_ana,1);
                toc
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Analytic Gibbs Sampling
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if ana_gibbs_on_off == 1
                tic
                disp('Analytic Gibbs Sampling')
                info_on_off = 1;
                %[samples_ana_gibbs, dist_ana_gibbs, loglik_ana_gibbs, nu_samples_ana_gibbs] = analytic_gibbs_new( is1, fn_evlas_hmc, clique_size, info_on_off, initial_point);
                [samples_ana_gibbs, dist_ana_gibbs, loglik_ana_gibbs, nu_samples_ana_gibbs] = Stretched_analytic_gibbs_new( is1, fn_evlas_hmc, clique_size, info_on_off, initial_point, dist_LBP);
                r_a_g = nu_samples_ana_gibbs/number_samples;
                mauss_ana_gibbs = mean(samples_ana_gibbs,1);
                toc
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Getting Errors
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % RMSE
            for j=2:number_samples/N
                
                %error_ana(q,j-1) = sqrt(mean(   (dist_truth - emp_dist(samples_ana(:,2:round(j*N*r_a-1)))) .^2));
                error_ana_dist(q,j-1) = sqrt(mean(   (dist_truth - mean(cat(1,dist_ana(:,2:round(j*N*r_a-1))),2))   .^2));
                error_ana_gibbs(q,j-1) = sqrt(mean(   (dist_truth - emp_dist(samples_ana_gibbs(:,2:round(j*N*r_a_g-1))))  .^2));
                if info_on_off ==1
                    error_ana_gibbs_dist(q,j-1) = sqrt(mean(  (dist_truth - mean(cat(1,dist_ana_gibbs(:,2:round(j*N*r_a_g-1))),2))  .^2));
                end
                %         error_uss(p,j-1) = sum(abs(dist_truth - emp_dist(samples_uss_line(:,1:floor(j*N*r_u)-10))));
                %         error_uss_inout(p,j-1) = sum(abs(dist_truth - emp_dist(samples_stepinout(1:round(j*N*r_u_inout)))));
                error_hmc(q,j-1) =sqrt(mean(   (dist_truth - emp_dist(samples_hmc(:,2:j*N-1)))  .^2));
                %error_cmh(q,j-1) = sqrt(mean(  (dist_truth - emp_dist(samples_CMH(:,2:round(j*N*r_cmh-1))))  .^2));
                error_cmh(q,j-1) = sqrt(mean(   (dist_truth - mean(cat(1,dist_cmh(:,2:round(j*N*r_cmh-1))),2))   .^2));
                
            end
            
            % Max Error
            for j=2:number_samples/N
                
                error_ana_1(q,j-1) = max(abs(dist_truth - emp_dist(samples_ana(:,2:round(j*N*r_a-1)))));
                error_ana_dist_1(q,j-1) = max(abs(dist_truth - mean(cat(1,dist_ana(:,2:round(j*N*r_a-1))),2)));
                error_ana_gibbs_1(q,j-1) = max(abs(dist_truth - emp_dist(samples_ana_gibbs(:,2:round(j*N*r_a_g-1)))));
                if info_on_off ==1
                    error_ana_gibbs_dist_1(q,j-1) = max(abs(dist_truth - mean(cat(1,dist_ana_gibbs(:,2:round(j*N*r_a_g-1))),2)));
                end
                %         error_uss(p,j-1) = sum(abs(dist_truth - emp_dist(samples_uss_line(:,1:floor(j*N*r_u)-10))));
                %         error_uss_inout(p,j-1) = sum(abs(dist_truth - emp_dist(samples_stepinout(1:round(j*N*r_u_inout)))));
                error_hmc_1(q,j-1) = max(abs(dist_truth - emp_dist(samples_hmc(:,2:j*N-1))));
                %error_cmh_1(q,j-1) = max(abs(dist_truth - emp_dist(samples_CMH(:,2:round(j*N*r_cmh-1)))));
                error_cmh_1(q,j-1) = max(abs(dist_truth - mean(cat(1,dist_cmh(:,2:round(j*N*r_cmh-1))),2)));
                
            end
            
%             % Total Variation
%                 error_ana_2(q,1) = total_var(is1, samples_ana, loglik_ana, weight);
%                 error_ana_dist_2(q,1) = tv_ana;
%                 error_ana_gibbs_2(q,1) = total_var(is1, samples_ana_gibbs, loglik_ana_gibbs', weight);
%                 error_ana_gibbs_dist_2(q,1) = tv_ana_gibbs;
%                 error_hmc_2(q,1) = total_var(is1, samples_hmc, loglik_hmc', weight);
%                 error_cmh_2(q,1) = total_var(is1, samples_CMH, log_lik_CMH, weight);
                            
            
        end
        
        
        % figure
        % semilogy(mauss_ana(1:end),'b')
        % hold on
        % semilogy(1:r_a:number_samples*r_a, mah(1:end),'r')
        % hold on
        % semilogy(1:r_a:number_samples*r_a, mabps(1:end), 'k')
        % title('Magnetization');
        % legend('Analytic','HMC', 'BPS')
        % xlabel('Iterations')
        % grid
        % box on
        
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % semilogyting relavant quantities
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % total Variation
        
%         display(['TV: Slice = ', num2str(mean(error_ana_2))])
%         display(['TV: Recycled Slice = ', num2str(mean(error_ana_dist_2))])
%         display(['TV: Gibbs = ', num2str(mean(error_ana_gibbs_2))])
%         display(['TV: RB Gibbs = ', num2str(mean(error_ana_gibbs_dist_2))])
%         display(['TV: HMC = ', num2str(mean(error_hmc_2))])
%         display(['TV: CMH = ', num2str(mean(error_cmh_2))])

    
        % Node Marginal Errors
        
        figure
%         semilogy(mean(error_ana,1),'b')
%         hold on
        semilogy(mean(error_ana_dist,1),'b')
        hold on
%         semilogy(mean(error_ana_gibbs,1),'g')
%         hold on
        semilogy(mean(error_ana_gibbs_dist,1),'g')
        hold on
        semilogy(mean(error_hmc,1),'r')
        hold on
        semilogy(mean(error_cmh,1),'k')
        h=legend('Recycled Slice', 'RB Gibbs', 'HMC', 'CMH');
        %set(h, 'Fontsize', 22)
        xlabel('Iterations')%, 'fontsize', 24)
        ylabel('RMSE (Marginals)')%, 'fontsize', 24)
        str=sprintf('Bias scale = %d', scale);
        title(str)%, 'fontsize', 24)
        filename = sprintf('Total Error Temp: %d and Scele: %d', temp, scale);
        saveas(gcf, filename, 'png')
        
        figure
        semilogy(mean(error_ana_dist_1,1),'b')
        hold on
        semilogy(mean(error_ana_gibbs_dist_1,1),'g')
        hold on
        semilogy(mean(error_hmc_1,1),'r')
        hold on
        semilogy(mean(error_cmh_1,1),'k')
        h=legend('Recycled Slice', 'RB Gibbs', 'HMC', 'CMH');
        %set(h, 'Fontsize', 22)
        xlabel('Iterations')%, 'fontsize', 24)
        ylabel('Max Error (Marginals)')%, 'fontsize', 24)
        title(str)%,'fontsize', 24)
        filename = sprintf('Max Error Temp: %d and Scele: %d', temp, scale);
        saveas(gcf, filename, 'png')

        
        % Log - likelihoods
        
        figure
        semilogy(loglik_ana(1:end),'b')
        hold on
        semilogy(loglik_ana_gibbs(1:end),'g')
        hold on
        semilogy(1:r_a:number_samples*r_a, loglik_hmc(1:end),'r')
        hold on
        % semilogy(1:r_cmh:number_samples*r_cmh, loglik_hmc(1:end),'r')
        % hold on
        semilogy(log_lik_CMH(1:end),'k')
        str2=sprintf('Logp:temp = %d, scale = %d', temp, scale);
        title(str2)
        saveas(gcf, str2, 'png')
        
        
    end
    
end

%% semilogyting Magnetism
% if fn_evals_uss_stepinout < fn_evals_uss
%
%     display('Stepping-In/Out works faster')
%     sc_ratio_stepinout = fn_evals_uss/fn_evals_uss_stepinout;
%     sc_ratio = fn_evlas_hmc/fn_evals_uss_stepinout;
%
% else
%
%     display('Only Stepping In works better')
%     sc_ratio_stepinout = fn_evals_uss_stepinout/fn_evals_uss;
%     sc_ratio = fn_evlas_hmc/fn_evals_uss;
%
% end

% % If stepping in stepping out works well
% subsemilogy(111)
% hold on
% semilogy(mauss_stepinout(1:end), 'k')
% semilogy(1:round(sc_ratio_stepinout):number_samples, mauss(1:ceil(number_samples/round(sc_ratio_stepinout))), 'b')
% semilogy(1:round(sc_ratio):number_samples, mah(1:ceil(number_samples/round(sc_ratio))), 'r')
% % semilogy (mah(1:end), 'r')
% title('Magnetization');
% legend('step In-Out','Uniform Slice Sampler', 'Exact-HMC')
% xlabel('Iterations')
% grid
% box on

% % If stepping in works well
% subsemilogy(111)
% hold on
% semilogy(1:round(sc_ratio_stepinout):number_samples, mauss_stepinout(1:ceil(number_samples/round(sc_ratio_stepinout))), 'g')
% semilogy(mauss(1:end), 'b')
% semilogy(1:round(sc_ratio):number_samples, mah(1:ceil(number_samples/round(sc_ratio))), 'r')
% semilogy (mah(1:end), 'r')
% title('Magnetization');
% legend('step In-Out','Uniform Slice Sampler', 'Exact-HMC')
% xlabel('Iterations')
% grid
% box on



%% This is for semilogyting acf for log-likelihoods for various samplers

% figure
% fig1 = subsemilogy(1,4,1);
% acf(loglik_ana',20);
% title('Analytic Slice Sampler')
% hold on
%
% fig2 = subsemilogy(1,4,2);
% acf(loglik_ana_gibbs',20);
% title('Analytic Gibbs Sampler')
% hold on
%
% % fig3 = subsemilogy(1,4,3);
% % acf(loglik_uss_stepinout',20);
% % title('Stepping In/Out')
% % hold on
%
% fig4 = subsemilogy(1,4,3);
% acf(loglik_hmc,20);
% title('HMC')
%
% fig5 = subsemilogy(1,4,4);
% acf(log_lik_CMH',20);
% title('Coordinate MH')



%% This semilogys the effective sample size from various samplers
% y=[62.5/(2*50),1, 500/(30*50)];
% bar(1, y(1),0.4,'r');
% hold on
% bar(2, y(2),0.4,'g');
% hold on
% bar(3, y(3),0.4,'b');
% legend('Uniform Slice Sampler','Step-in/out', 'Exact-HMC')
%
% colors = hsv(numel(y));
% for i = 1:numel(y)
%     bar(i, y(i), 'parent', 'facecolor', colors(i));
% end
%
% hb = bar(y);
% set(hb(1),'facecolor','g')
%
%
% bar(y(1),'b')
% hold on
% bar(y(2),'g'),bar(y(3),'y')



%% This semilogys the magnetism of points as the samplers move through the space
% subsemilogy(111)
% hold on
% semilogy(mauss_ana(1:end), 'k')
% semilogy(mabps(1:end),'g')
% % semilogy(mauss(1:end), 'b')
% % semilogy(mauss_stepinout(1:end), 'g')
% % semilogy(1:round(sc_ratio):number_samples,mah(1:ceil(number_samples/round(sc_ratio))), 'r')
% semilogy (mah(1:end), 'r')
% title('Magnetization');
% legend('Uniform Slice Sampler', 'Exact-HMC')
% xlabel('Iterations')
% grid
% box on


%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Uniform Slice Sampling
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     if uniform_SS_on_off == 1
%         tic
%         disp('Uniform Slice Sampling')
%         [samples_uss_line, dist_uss, loglik_uss, fn_evals_uss, nu_samples_uss] = ussSampler(is1, 1, 0, fn_evlas_hmc, clique_size);
%         r_u = nu_samples_uss/number_samples;
%         mauss = mean(samples_uss_line,1);
%         toc
%         %     slice sampling on a circle
%         %     [samples_uss_ess, loglik_uss_ess, fn_evals_uss_ess]= ussSampler(is1, number_samples, d, number_chains,0);
%         %     mauss_ess = mean(samples_uss_ess,1);
%     end
%
%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Uniform Slice Sampling with Stepping In, Stepping Out
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     if uss_stepinout_on_off == 1
%         disp('Uniform Slice Sampling with Stepping In and Stepping Out')
%         [samples_stepinout, loglik_uss_stepinout, fn_evals_uss_stepinout, nu_samples_inout]= uss_step_InOut(is1, d, fn_evlas_hmc, clique_size);
%         mauss_stepinout = mean(samples_stepinout,1);
%         r_u_inout = nu_samples_inout/number_samples;
%     end
