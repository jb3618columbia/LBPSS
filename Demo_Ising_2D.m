% Implementation of Different Algorithms for 1D and 2D Ising Models
% Parameters
% d=10;  % for 1D Ising models
d=16;   % for 2D Ising models, this creates a sqrt(d) X sqrt(d) lattice
temp_vec=[10*pi];
scale_vec = [1,3,5,7,10,20,50];
data_matrix = zeros(6, length(scale_vec)); % for heat map

    
for u=1:1:length(temp_vec)
    temp = temp_vec(u);
    
    for v=1:1:length(scale_vec);
        scale = scale_vec(v);
        % is1 = Ising1D_new(d,temp, scale);  % Create 1D Ising Object
        % clique_size=2; %Clique size
        is1 = Ising2D(sqrt(d),temp, scale);  % Create 2D Ising Object
        clique_size=2; %Clique size
        number_samples = 5000;
        num_examples = 20;
        initial_point = sign(normrnd(0,1,d,1));
        
        % Algorithms:
        truth = true;
        ana_on_off = true;
        ana_gibbs_on_off =  true;
        exact_hmc_on_off = true;
        CMH_on_off = true;
        info_on_off = true;
        % With LBP
        ana_lbp_on_off = true;
        ana_lbp_gibbs_on_off = true;
        
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
           
        dist_LBP = get_LBP_marginals( -is1.bias , -is1.M);
        
        % Error for various samplers        
        N=250; % Check error after every N equivalent iterations
        error_ana = zeros(num_examples, number_samples/N);
        error_ana_gibbs = zeros(num_examples, number_samples/N);
        error_hmc = zeros(num_examples, number_samples/N);
        error_cmh = zeros(num_examples, number_samples/N);
        error_ana_gibbs_lbp = zeros(num_examples, number_samples/N);
        error_ana_lbp = zeros(num_examples, number_samples/N);

        for q=1:num_examples
            q
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Exact-HMC
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if exact_hmc_on_off == 1
                tic
                disp('Exact_HMC')
                t = 1.5; T=t*pi;
                [samples_hmc, loglik_hmc, energy_hmc] = HMC_binary(is1,T,number_samples, initial_point);
                fn_evlas_hmc = number_samples*((is1.dim)*t + (is1.dim) + clique_size*(is1.dim));
                toc
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Analytic Slice Sampling
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ana_on_off == 1
                tic
                disp('Analytic Slice Sampling')
                [samples_ana, dist_ana, loglik_ana, nu_samples_ana]= analytic_slice_new( is1, fn_evlas_hmc, clique_size, initial_point);
                r_a = nu_samples_ana/number_samples;
                toc
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Analytic Slice Sampling with LBP
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ana_lbp_on_off == 1
                tic
                disp('Analytic Slice Sampling with LBP')
                [samples_ana_lbp, dist_ana_lbp, loglik_ana_lbp, nu_samples_ana_lbp]= Stretched_analytic_slice_new( is1, fn_evlas_hmc, clique_size, info_on_off, initial_point, dist_LBP);%2*(dist_LBP>0.5)-1%analytic_slice_new( is1, fn_evlas_hmc, clique_size, initial_point);%
                r_a_lbp = nu_samples_ana_lbp/number_samples;
                toc
            end
        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Coordiante MH
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if CMH_on_off == 1
                tic
                disp('Coordinate MH')
                [samples_CMH, log_lik_CMH, nu_samples_cmh] = CMH(is1, fn_evlas_hmc, clique_size, initial_point);
                r_cmh = nu_samples_cmh/number_samples;
                toc
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Analytic Gibbs Sampling
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ana_gibbs_on_off == 1
                tic
                disp('Analytic Gibbs Sampling')
                [samples_ana_gibbs, dist_ana_gibbs, loglik_ana_gibbs, nu_samples_ana_gibbs] = analytic_gibbs_new( is1, fn_evlas_hmc, clique_size, info_on_off, initial_point);
                r_a_g = nu_samples_ana_gibbs/number_samples;
                toc
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Analytic Gibbs Sampling with LBP
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ana_lbp_gibbs_on_off == 1
                tic
                disp('Analytic Gibbs Sampling with LBP')
                [samples_ana_gibbs_lbp, dist_ana_gibbs_lbp, loglik_ana_gibbs_lbp, nu_samples_ana_gibbs_lbp] = Stretched_analytic_gibbs_new( is1, fn_evlas_hmc, clique_size, info_on_off, initial_point, dist_LBP);
                r_a_g_lbp = nu_samples_ana_gibbs_lbp/number_samples;
                toc
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % RMSE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            error_ana(q,1) = sqrt(mean(   (dist_truth - emp_dist(initial_point))   .^2));
            error_ana_gibbs(q,1) = sqrt(mean(  (dist_truth - emp_dist(initial_point))  .^2));
            error_hmc(q,1) =sqrt(mean(   (dist_truth - emp_dist(initial_point))  .^2));
            error_cmh(q,1) = sqrt(mean(  (dist_truth - emp_dist(initial_point))  .^2));
            
            error_ana_lbp(q,1) = sqrt(mean(  (dist_truth - emp_dist(initial_point))  .^2));
            error_ana_gibbs_lbp(q,1) = sqrt(mean(  (dist_truth - emp_dist(initial_point))  .^2));
            
            for j=2:number_samples/N
                
                error_ana(q,j) = sqrt(mean(   (dist_truth - mean(cat(1,dist_ana(:,2:round(j*N*r_a-1))),2))   .^2));
                error_ana_gibbs(q,j) = sqrt(mean(  (dist_truth - mean(cat(1,dist_ana_gibbs(:,2:round(j*N*r_a_g-1))),2))  .^2));
                error_hmc(q,j) =sqrt(mean(   (dist_truth - emp_dist(samples_hmc(:,2:j*N-1)))  .^2));
                error_cmh(q,j) = sqrt(mean(  (dist_truth - emp_dist(samples_CMH(:,2:round(j*N*r_cmh-1))))  .^2));
                
                error_ana_lbp(q,j) = sqrt(mean(   (dist_truth - mean(cat(1,dist_ana_lbp(:,2:round(j*N*r_a-1))),2))   .^2));
                error_ana_gibbs_lbp(q,j) = sqrt(mean(  (dist_truth - mean(cat(1,dist_ana_gibbs_lbp(:,2:round(j*N*r_a_g-1))),2))  .^2));
            end
            
        end
        
        % For Heat Maps
         data_matrix(1, v) = mean(error_hmc(:,number_samples/N));
         data_matrix(2, v) = mean(error_cmh(:,number_samples/N));
         data_matrix(3, v) = mean(error_ana(:,number_samples/N));
         data_matrix(4, v) = mean(error_ana_lbp(:,number_samples/N));
         data_matrix(5, v) = mean(error_ana_gibbs(:,number_samples/N));
         data_matrix(6, v) = mean(error_ana_gibbs_lbp(:,number_samples/N));
        
        % Node Marginal Errors
%         figure
%         semilogy(0:length(error_ana)-1, mean(error_ana,1), 'Color', [152,78,163]/255);
%         hold on
%         semilogy(0:length(error_ana_gibbs)-1, mean(error_ana_gibbs,1), 'Color', [77,175,74]/255);
%         hold on
%         semilogy(0:length(error_hmc)-1, mean(error_hmc,1), 'Color', [228,26,28]/255);
%         hold on
%         semilogy(0:length(error_cmh)-1, mean(error_cmh,1),'Color', [55,126,184]/255);
%         hold on
%         semilogy(0:length(error_ana_lbp)-1, mean(error_ana_lbp,1), '--', 'Color', [152,78,163]/255);
%         hold on
%         semilogy(0:length(error_ana_gibbs_lbp)-1, mean(error_ana_gibbs_lbp,1),'--', 'Color', [77,175,74]/255);
%         h=legend('AAS', 'AAG', 'HMC', 'CMH','AAS LBP', 'AAG LBP');
%         set(h);
%         xlabel('Function Evaluations');
%         ylabel('RMSE (Marginals)');
%         str=sprintf('Bias scale = %d', scale);
%         title(str);
%         set(gcf,'units','points','position',[10,10,800,800]);
%         name = strcat('RMSE Temp', num2str(temp), 'Bias', num2str(scale), '.fig');
%         path = '/Users/Jalaj/Documents/Github - LBPSS/New Outputs';
%         savefig(gcf, fullfile(path, name))
        
        
%         figure
%         plot(loglik_ana(1:end),'b')
%         hold on
%         plot(loglik_ana_lbp(1:end),'m') 
%         hold on 
%         plot(loglik_ana_gibbs(1:end),'g')
%         hold on
%         plot(loglik_ana_gibbs_lbp(1:end),'c')
%         hold on
%         plot(1:r_a:number_samples*r_a, loglik_hmc(1:end),'r')
%         hold on
%         plot(log_lik_CMH(1:end),'k')
%         str2=sprintf('Logp:temp = %d, scale = %d', temp, scale);
%         title(str2)
%         saveas(gcf, str2, 'png')

         
         
        
    end
    
    colormap('hot')
    imagesc(data_matrix);
    xlabel('Bias');
    ylabel('Algorithms')
    title('RMSE (Marginals)');
    ax = gca;
    ax.YTickLabel = {'HMC','CMH','AAS','AAS LBP','AAG', 'AAG LBP'};
    ax.XTickLabel = {'', '1', '', '3','','5' '', '7','','10', '', '20','','50'};
    h=colorbar;

end
            
        
        
        
        
        
        