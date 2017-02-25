function [ ] = plotZerror(logZ_est_true, logZ_est_HMC, logZ_est_CMH, logZ_est_AAS, logZ_est_AAG)

figure
% plot(logZ_est_true, 'm')
% hold on
plot(abs(logZ_est_true - logZ_est_HMC), 'r')
hold on
plot(abs(logZ_est_true - logZ_est_CMH), 'b')
hold on
plot(abs(logZ_est_true - logZ_est_AAS), 'g')
hold on
plot(abs(logZ_est_true - logZ_est_AAG), 'k')
legend('HMC','CMH','AAS','AAG')
end

