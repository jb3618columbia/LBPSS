% Implementation of Different Algorithms for Ising Models
% Models:
% 1) 1D Ising Models
% 2) 2D Ising Models

% clear
% Parameters
d=100; 
temp=5*pi;
is1 = Ising1D_new(d,temp);  % Create 1D Ising Object 
% is1 = Ising1D_rand_weight(d,temp); 
clique_size=2; %Clique size 
number_samples = 500;
num_examples = 1;
% initial_point = sign(normrnd(0,1,d,1));
initial_point = ones(d,1);

% Algorithms:
ana_on_off = true;
ana_gibbs_on_off = false;
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
    [samples_true, loglik_true, energy_true] = HMC_binary(is1,T,num_samp_truth, initial_point);
%     dist_truth = emp_dist(samples_true);
%     fn_evlas = 10000*2*d;
%     [samples_true, loglik_ana, fn_evals_ana, nu_samples_ana]= ussSampler(is1, d, 0, 1, fn_evlas);
%     dist_truth = emp_dist(samples_true);
%     
end
dist_truth = 0.5*ones(d,1);

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
        [samples_hmc, loglik_hmc, energy_hmc] = HMC_binary(is1,T,number_samples, initial_point);
        mah = mean(samples_hmc,1);
        fn_evlas_hmc = number_samples*((is1.dim)*t + (is1.dim) + clique_size*(is1.dim));
        toc
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Uniform Slice Sampling
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if uniform_SS_on_off == 1
        tic
        disp('Uniform Slice Sampling')
        [samples_uss_line, dist_uss, loglik_uss, fn_evals_uss, nu_samples_uss] = ussSampler(is1, 1, 0, fn_evlas_hmc, clique_size);
        r_u = nu_samples_uss/number_samples;
        mauss = mean(samples_uss_line,1);
        toc
        %     slice sampling on a circle
        %     [samples_uss_ess, loglik_uss_ess, fn_evals_uss_ess]= ussSampler(is1, number_samples, d, number_chains,0);
        %     mauss_ess = mean(samples_uss_ess,1);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Analytic Slice Sampling
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ana_on_off == 1
        tic
        disp('Analytic Slice Sampling')
        %         [samples_ana, dist_ana, loglik_ana, fn_evals_ana, nu_samples_ana]= ussSampler(is1, 0, 1, fn_evlas_hmc, clique_size, initial_point);
        [samples_ana, dist_ana, loglik_ana, nu_samples_ana]= analytic_slice_new( is1, fn_evlas_hmc, clique_size, initial_point);
        r_a = nu_samples_ana/number_samples;
        mauss_ana = mean(samples_ana,1);
        toc
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Uniform Slice Sampling with Stepping In, Stepping Out
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if uss_stepinout_on_off == 1
        disp('Uniform Slice Sampling with Stepping In and Stepping Out')
        [samples_stepinout, loglik_uss_stepinout, fn_evals_uss_stepinout, nu_samples_inout]= uss_step_InOut(is1, d, fn_evlas_hmc, clique_size);
        mauss_stepinout = mean(samples_stepinout,1);
        r_u_inout = nu_samples_inout/number_samples;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Metropolis with random flip proposal
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if CMH_on_off == 1
        tic
        disp('Coordinate MH')
        [samples_CMH, log_lik_CMH, nu_samples_cmh] = CMH(is1, fn_evlas_hmc, clique_size, initial_point);
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
        info_on_off = 1;
        [samples_ana_gibbs, dist_ana_gibbs, loglik_ana_gibbs, nu_samples_ana_gibbs]= analytic_gibbs( is1, fn_evlas_hmc, clique_size, info_on_off, initial_point);
        r_a_g = nu_samples_ana_gibbs/number_samples;
        mauss_ana_gibbs = mean(samples_ana_gibbs,1);
        toc
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Getting Errors
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for j=2:number_samples/N
        
        error_ana(p,j-1) = sum(abs(dist_truth - emp_dist(samples_ana(:,2:round(j*N*r_a-1)))));
        error_ana_dist(p,j-1) = sum(abs(dist_truth - mean(cat(1,dist_ana(:,2:round(j*N*r_a-1))),2)));
%         error_ana_gibbs(p,j-1) = sum(abs(dist_truth - emp_dist(samples_ana_gibbs(:,2:round(j*N*r_a_g-1)))));
%         if info_on_off ==1
%             error_ana_gibbs_dist(p,j-1) = sum(abs(dist_truth - mean(cat(1,dist_ana_gibbs(:,2:round(j*N*r_a_g-1))),2)));
%         end
%         error_uss(p,j-1) = sum(abs(dist_truth - emp_dist(samples_uss_line(:,1:floor(j*N*r_u)-10))));
%         error_uss_inout(p,j-1) = sum(abs(dist_truth - emp_dist(samples_stepinout(1:round(j*N*r_u_inout)))));
        error_hmc(p,j-1) = sum(abs(dist_truth - emp_dist(samples_hmc(:,2:j*N-1))));
        error_cmh(p,j-1) = sum(abs(dist_truth - emp_dist(samples_CMH(:,2:round(j*N*r_cmh-1)))));
        
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
% plot(mean(error_ana_gibbs,1),'g')
% hold on
% plot(mean(error_ana_gibbs_dist,1),'c')
% hold on
% plot(mean(error_uss,1),'c')
% hold on
% plot(error_uss_inout,'g')
% hold on 
% plot(mean(error_bps,1),'g')
% hold on 
plot(mean(error_hmc,1)+0.5,'r')
hold on
plot(mean(error_cmh,1),'k')
legend('Ana Slice', 'Recycled Ana Slice', 'Ana Gibbs', 'Recycled Ana Gibbs', 'HMC', 'CMH')
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
% plot(1:r_cmh:number_samples*r_cmh, loglik_hmc(1:end),'r')
% hold on
plot(log_lik_CMH(1:end),'k')




%% Plotting Magnetism
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
% subplot(111)
% hold on 
% plot(mauss_stepinout(1:end), 'k')
% plot(1:round(sc_ratio_stepinout):number_samples, mauss(1:ceil(number_samples/round(sc_ratio_stepinout))), 'b')
% plot(1:round(sc_ratio):number_samples, mah(1:ceil(number_samples/round(sc_ratio))), 'r')
% % plot (mah(1:end), 'r')
% title('Magnetization');
% legend('step In-Out','Uniform Slice Sampler', 'Exact-HMC')
% xlabel('Iterations')
% grid 
% box on 

% % If stepping in works well
% subplot(111)
% hold on
% plot(1:round(sc_ratio_stepinout):number_samples, mauss_stepinout(1:ceil(number_samples/round(sc_ratio_stepinout))), 'g')
% plot(mauss(1:end), 'b')
% plot(1:round(sc_ratio):number_samples, mah(1:ceil(number_samples/round(sc_ratio))), 'r')
% plot (mah(1:end), 'r')
% title('Magnetization');
% legend('step In-Out','Uniform Slice Sampler', 'Exact-HMC')
% xlabel('Iterations')
% grid
% box on



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



%% This plots the effective sample size from various samplers
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



