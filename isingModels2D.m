% Implementation of Different Algorithms for 2D Ising Models


% Parameters
% d=5; 
% temp=5*pi;
% clique_size=4; % Clique size
% is2 = Ising2D(d,temp);  % Create 2D Ising Object 

d=5; 
temp=5*pi;
k=d; % Clique size
is2 = Ising2D_rand_weight(d,temp);  % Create a random Ising 

number_samples = 1500;
num_examples = 5;

% Algorithms:
ana_on_off = true;
ana_gibbs_on_off = true;
exact_hmc_on_off = true;
uss_stepinout_on_off = false;
uniform_SS_on_off = false;
CMH_on_off = true;

ground_truth = false;
lbpSS_on_off = false;


%%%%%%%%%%%%%%%%%%%%%%%%
% Approx. ground truth
%%%%%%%%%%%%%%%%%%%%%%%%
if ground_truth == 1
    display('Getting samples to approximate grund truth')
    t = 1.5; T=t*pi;
    num_samp_truth = 10000;
    fn_evlas = num_samp_truth*( (is2.dim)*t + (is2.dim) + clique_size*(is2.dim));
    [samples_true, loglik_true, energy_true] = CMH(is2, fn_evlas, clique_size);
%     dist_truth = emp_dist(samples_true); 
end
dist_truth = 0.5*ones(d*d,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparing Samplers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=50; % Check error after every N equivalent iterations
error_ana = zeros(num_examples, number_samples/N -1);
error_ana_dist = zeros(num_examples, number_samples/N -1);
error_ana_gibbs = zeros(num_examples, number_samples/N -1);
error_ana_gibbs_dist = zeros(num_examples, number_samples/N -1);
error_uss = zeros(num_examples, number_samples/N-1);
error_uss_inout = zeros(num_examples, number_samples/N-1);
error_hmc = zeros(num_examples, number_samples/N -1);
error_bps = zeros(num_examples, number_samples/N -1);
error_cmh = zeros(num_examples, number_samples/N -1);


for p=1:num_examples
    
    p
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Exact-HMC
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if exact_hmc_on_off == 1
        tic
        disp('Exact_HMC')
        t = 1.5; T=t*pi;
        [samples_hmc, loglik_hmc, energy_hmc] = HMC_binary(is2,T,number_samples);
        mah = mean(samples_hmc,1);
        fn_evlas_hmc = number_samples*((is2.dim)*t + (is2.dim) + clique_size*(is2.dim));
        toc
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Uniform Slice Sampling
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if uniform_SS_on_off == 1
        tic
        disp('Uniform Slice Sampling')
        [samples_uss_line, dist_uss, loglik_uss, fn_evals_uss, nu_samples_uss] = ussSampler(is2, 1, 0, fn_evlas_hmc, clique_size);
        r_u = nu_samples_uss/number_samples;
        mauss = mean(samples_uss_line,1);
        toc
        %     slice sampling on a circle
        %     [samples_uss_ess, loglik_uss_ess, fn_evals_uss_ess]= ussSampler(is2, number_samples, d, number_chains,0);
        %     mauss_ess = mean(samples_uss_ess,1);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Analytic Slice Sampling
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ana_on_off == 1
        tic
        disp('Analytic Slice Sampling')
        [samples_ana, dist_ana, loglik_ana, fn_evals_ana, nu_samples_ana]= ussSampler(is2, 0, 1, fn_evlas_hmc, clique_size);
        r_a = nu_samples_ana/number_samples;
        mauss_ana = mean(samples_ana,1);
        toc
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Uniform Slice Sampling with Stepping In, Stepping Out
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if uss_stepinout_on_off == 1
        disp('Uniform Slice Sampling with Stepping In and Stepping Out')
        [samples_stepinout, loglik_uss_stepinout, fn_evals_uss_stepinout, nu_samples_inout]= uss_step_InOut(is2, fn_evlas_hmc, clique_size);
        mauss_stepinout = mean(samples_stepinout,1);
        r_u_inout = nu_samples_inout/number_samples;
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Coordinate MH
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if CMH_on_off == 1
        tic
        disp('Coordinate MH')
        [samples_CMH, log_lik_CMH, nu_samples_cmh] = CMH(is2, fn_evlas_hmc, clique_size);
        r_cmh = nu_samples_cmh/number_samples;
        macmh = mean(samples_CMH,1);
        toc
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Analytic Gibbs Sampling
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ana_gibbs_on_off == 1
        tic
        disp('Analytic Gibbs Sampling')
        info_on_off = 0;
        [samples_ana_gibbs, dist_ana_gibbs, loglik_ana_gibbs, nu_samples_ana_gibbs]= analytic_gibbs( is2, fn_evlas_hmc, clique_size ,info_on_off);
        r_a_g = nu_samples_ana_gibbs/number_samples;
        mauss_ana_gibbs = mean(samples_ana_gibbs,1);
        toc
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Getting Errors
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for j=2:number_samples/N
        
        error_ana(p,j-1) = sum(abs(dist_truth - emp_dist(samples_ana(:,1:round(j*N*r_a-1)))));
        error_ana_dist(p,j-1) = sum(abs(dist_truth - mean(cat(1,dist_ana(:,1:round(j*N*r_a-1))),2)));
        error_ana_gibbs(p,j-1) = sum(abs(dist_truth - emp_dist(samples_ana_gibbs(:,1:round(j*N*r_a_g-1)))));
        if info_on_off ==1
            error_ana_gibbs_dist(p,j-1) = sum(abs(dist_truth - mean(cat(1,dist_ana_gibbs(:,1:round(j*N*r_a_g-1))),2)));
        end
%         error_uss(p,j-1) = sum(abs(dist_truth - emp_dist(samples_uss_line(:,1:floor(j*N*r_u)-10))));
%         error_uss_inout(p,j-1) = sum(abs(dist_truth - emp_dist(samples_stepinout(1:round(j*N*r_u_inout)))));
        error_hmc(p,j-1) = sum(abs(dist_truth - emp_dist(samples_hmc(:,1:j*N-1))));
        error_cmh(p,j-1) = sum(abs(dist_truth - emp_dist(samples_CMH(:,1:round(j*N*r_cmh-1)))));
        
    end
    
end


% figure 
% plot(mauss_ana(1:end),'b')
% hold on
% plot(1:r_a:number_samples*r_a, mah(1:end),'r')
% hold on 
% plot(1:r_a:number_samples*r_a, mabps(1:end), 'k')
% title('Magnetization');
% legend('Analytic','HMC', 'BPS')
% xlabel('Iterations')
% grid 
% box on



figure
plot(mean(error_ana,1),'b')
hold on
plot(mean(error_ana_dist,1),'m')
hold on
plot(mean(error_ana_gibbs,1),'g')
hold on
% plot(mean(error_ana_gibbs_dist,1),'c')
% hold on
% plot(mean(error_uss,1),'c')
% hold on
% plot(error_uss_inout,'g')
% hold on 
plot(mean(error_hmc,1)+0.5,'r')
hold on
plot(mean(error_cmh,1),'k')
legend('Ana Slice', 'Recycled Ana Slice', 'Ana Gibbs', 'HMC', 'CMH')
xlabel('Iterations')
ylabel('Error (Node Marginals)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting relavant quantities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Plotting log likelihoods
figure
plot(loglik_ana(1:end),'b')
hold on
plot(loglik_ana_gibbs(1:end),'g')
hold on
plot(1:r_a:number_samples*r_a, loglik_hmc(1:end),'r')
hold on
plot(1:r_cmh:number_samples*r_cmh, loglik_hmc(1:end),'r')
hold on
plot(log_lik_CMH(1:end),'k')


%% This is for plotting acf for log-likelihoods for various samplers

% figure
% fig1 = subplot(1,4,1);
% acf(loglik_ana',20);
% title('Analytic Slice Sampler')
% hold on
% 
% fig2 = subplot(1,4,2);
% acf(loglik_ana_gibbs',20);
% title('Analytic Gibbs Sampler')
% hold on 
% 
% % fig3 = subplot(1,4,3);
% % acf(loglik_uss_stepinout',20);
% % title('Stepping In/Out')
% % hold on 
% 
% fig4 = subplot(1,4,3);
% acf(loglik_hmc,20);
% title('HMC')
% 
% fig5 = subplot(1,4,4);
% acf(log_lik_CMH',20);
% title('Coordinate MH')

     
%% This plots the magnetism of points as the samplers move through the space
% subplot(111)
% hold on 
% plot(mauss_ana(1:end), 'k')
% plot(mabps(1:end),'g')
% % plot(mauss(1:end), 'b')
% % plot(mauss_stepinout(1:end), 'g')
% % plot(1:round(sc_ratio):number_samples,mah(1:ceil(number_samples/round(sc_ratio))), 'r')
% plot (mah(1:end), 'r')
% title('Magnetization');
% legend('Uniform Slice Sampler', 'Exact-HMC')
% xlabel('Iterations')
% grid
% box on



