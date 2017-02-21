function [ ] = plotRealdata()
% This is to plot the histigrams for samples obtained in outer MH loop.
% Not including AAS as we have too less samples here.

path = '/Users/Jalaj/Documents/Github - LBPSS/Real_Data/run1/';
lab = char(['tru','HMC', 'CMH', 'AAS', 'AAG']);
figure
for i=1:15
    subplot(3,5,i);
    for j=1:length(lab)/3
        fileName = sprintf('%s%s.mat', path, lab((j-1)*3 + 1:j*3));
        load(fileName);
    end
    pd_true = fitdist(samples_true(:,i),'Normal');
    y_true = pdf(pd_true, samples_true(:,i));
    pd_HMC = fitdist(samples_HMC(:,i),'Normal');
    y_HMC = pdf(pd_HMC, samples_HMC(:,i));
    pd_CMH = fitdist(samples_CMH(:,i),'Normal');
    y_CMH = pdf(pd_CMH, samples_CMH(:,i));
%     pd_AAS = fitdist(samples_AAS(:,i),'Normal');
%     y_AAS = pdf(pd_AAS, samples_AAS(:,i));
    pd_AAG = fitdist(samples_AAG(:,i),'Normal');
    y_AAG = pdf(pd_AAG, samples_AAG(:,i));
    
    plot(samples_true(:,i), y_true,'LineWidth',2)
    hold on
    plot(samples_HMC(:,i), y_HMC,'LineWidth',2)
    hold on
    plot(samples_CMH(:,i), y_CMH,'LineWidth',2)
    hold on
%     plot(samples_AAS(:,i), y_AAS,'LineWidth',2)
%     hold on
    plot(samples_AAG(:,i), y_AAG,'LineWidth',2)
    
    %     histogram(samples_true(:,i))
    %     hold on
    %     histogram(samples_HMC(:,i))
    %     hold on
    %     histogram(samples_CMH(:,i))
    %     hold on
    %     histogram(samples_AAS(:,i))
    %     hold on
    %     histogram(samples_AAG(:,i))
    
end
legend('True','HMC','CMH','AAG')

end

