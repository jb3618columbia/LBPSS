% Implementation of Slice sampling for Ising Models

% To do:
% 
% Models:
% 1) 1D Ising Models
% 2) 2D Ising Models
% 
% Algorithms:
% 
% 1) Naive slice sampler (Uniform proposals)
% 2) LBP SS
% 
% Comparative Algorithms:
% 
% 1) Exact HMC 
% 2) Metropolis Sampler


% Parameters
d=400; 
temp=12*pi;

is1 = Ising1D(d,temp);  % Create 1D Ising Object 
number_samples = 1000;
number_chains = 1;  % Maybe used in Effective Sample Size Computations

uniform_SS_on_off = true;
uss_stepinout_on_off = true;
lbpSS_on_off = false;
exact_hmc_on_off = false;
metroBinary_on_off = false;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uniform Slice Sampling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if uniform_SS_on_off == 1
    
    
    tic
    disp('Uniform Slice Sampling')
    [samples_uss_line, loglik_uss, fn_evals_uss]= ussSampler(is1, number_samples, d, number_chains,1);
    mauss = mean(samples_uss_line,1);
    
    % slice sampling on a circle
%     [samples_uss_ess, loglik_uss_ess, fn_evals_uss_ess]= ussSampler(is1, number_samples, d, number_chains,0);
%     mauss_ess = mean(samples_uss_ess,1);  

    toc
    
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uniform Slice Sampling with Stepping In, Stepping Out 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if uss_stepinout_on_off == 1
    
    
    tic
    disp('Uniform Slice Sampling with Stepping In and Stepping Out')
    [samples_stepinout, loglik_uss_stepinout, fn_evals_uss_stepinout]= uss_step_InOut(is1, number_samples, d, number_chains,1);
    mauss_stepinout = mean(samples_stepinout,1);

%     [samples_uss_ess, loglik_uss_ess]= ussSampler(is1, number_samples, d, number_chains,0);
    toc
    
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exact-HMC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exact_hmc_on_off == 1
    
    tic
    disp('Exact_HMC')
    t = 12.5;
    T=t*pi;
    [samples_hmc, energy_hmc, loglik_hmc] = HMC_binary(is1,T,number_samples);
    mah = mean(samples_hmc,1);
    fn_evlas_hmc = number_samples*(3*d*t + d);
    toc
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Metropolis with random flip proposal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if metroBinary_on_off == 1
    
    disp('Metroplois')
    
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting relavant quantities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if fn_evals_uss_stepinout < fn_evals_uss
    
    display('Stepping-In/Out works faster')
    sc_ratio_stepinout = fn_evals_uss/fn_evals_uss_stepinout;
    sc_ratio = fn_evlas_hmc/fn_evals_uss_stepinout;
   
else 
   
    display('Only Stepping In works better')
    sc_ratio_stepinout = fn_evals_uss_stepinout/fn_evals_uss;
    sc_ratio = fn_evlas_hmc/fn_evals_uss;  
    
end

% If stepping in stepping out works well
subplot(111)
hold on 
plot(mauss_stepinout(1:end), 'k')
plot(1:round(sc_ratio_stepinout):number_samples, mauss(1:ceil(number_samples/round(sc_ratio_stepinout))), 'b')
plot(1:round(sc_ratio):number_samples, mah(1:ceil(number_samples/round(sc_ratio))), 'r')
% plot (mah(1:end), 'r')
title('Magnetization');
legend('step In-Out','Uniform Slice Sampler', 'Exact-HMC')
xlabel('Iterations')
grid 
box on 

% If stepping in works well
subplot(111)
hold on
plot(1:round(sc_ratio_stepinout):number_samples, mauss_stepinout(1:ceil(number_samples/round(sc_ratio_stepinout))), 'g')
plot(mauss(1:end), 'b')
plot(1:round(sc_ratio):number_samples, mah(1:ceil(number_samples/round(sc_ratio))), 'r')
plot (mah(1:end), 'r')
title('Magnetization');
legend('step In-Out','Uniform Slice Sampler', 'Exact-HMC')
xlabel('Iterations')
grid
box on

% loglik_uss = loglik_uss(1:round(sc_ratio):number_samples);
subplot(1,3,1)
acf(loglik_uss',20);
title('Uniform Slice Sampler')
hold on 

% loglik_uss_stepinout = loglik_uss_stepinout(1:round(sc_ratio):number_samples);
subplot(1,3,2)
acf(loglik_uss_stepinout',20);
title('Stepping In/Out')
hold on 

subplot(1,3,3)
acf(loglik_hmc,20);
title('HMC')

y=[62.5/(2*50),1, 500/(30*50)];
bar(1, y(1),0.4,'r'); 
hold on
bar(2, y(2),0.4,'g'); 
hold on 
bar(3, y(3),0.4,'b'); 
legend('Uniform Slice Sampler','Step-in/out', 'Exact-HMC')

colors = hsv(numel(y));
for i = 1:numel(y)
    bar(i, y(i), 'parent', 'facecolor', colors(i));
end

hb = bar(y);
set(hb(1),'facecolor','g') 


bar(y(1),'b')
hold on
bar(y(2),'g'),bar(y(3),'y')




% bwmap = [.97 .97 1 ; 0 0 0];
% colormap(bwmap);
% 
% 
% imagesc(samples_uss_line)
% title('USS');
% xlabel('Iteration')
% 
% imagesc(samples_hmc)
% title('HMC');
% xlabel('Iteration')
     

% subplot(111)
% hold on 
% plot(mauss(1:end), 'b')
% % plot(mauss_stepinout(1:end), 'g')
% plot(1:round(sc_ratio):number_samples,mah(1:ceil(number_samples/round(sc_ratio))), 'r')
% % plot (mah(1:end), 'r')
% title('Magnetization');
% legend('Uniform Slice Sampler', 'Exact-HMC')
% xlabel('Iterations')
% grid
% box on


% subplot(212)
% hold on 
% plot (mah(1:end))
% grid
% box off

