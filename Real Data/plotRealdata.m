function [ ] = plotRealdata(samples_true, samples_HMC, samples_CMH, samples_AAS, samples_AAG)

figure
for i=1:15
    subplot(3,5,i);
    %     plot(samples_true(:,i), 'm')
    %     hold on
    %     plot(samples_HMC(:,i), 'r')
    %     hold on
    %     plot(samples_CMH(:,i), 'b')
    %     hold on
    %     plot(samples_AAS(:,i), 'g')
    %     hold on
    %     plot(samples_AAG(:,i), 'k')
    pd_true = fitdist(samples_true(:,i),'Normal');
    y_true = pdf(pd_true, samples_true(:,i));
    pd_HMC = fitdist(samples_HMC(:,i),'Normal');
    y_HMC = pdf(pd_HMC, samples_HMC(:,i));
    pd_CMH = fitdist(samples_CMH(:,i),'Normal');
    y_CMH = pdf(pd_CMH, samples_CMH(:,i));
    pd_AAS = fitdist(samples_AAS(:,i),'Normal');
    y_AAS = pdf(pd_AAS, samples_AAS(:,i));
    pd_AAG = fitdist(samples_AAG(:,i),'Normal');
    y_AAG = pdf(pd_AAG, samples_AAG(:,i));
    
    plot(samples_true(:,i), y_true,'LineWidth',2)
    hold on
    plot(samples_HMC(:,i), y_HMC,'LineWidth',2)
    hold on
    plot(samples_CMH(:,i), y_CMH,'LineWidth',2)
    hold on
    plot(samples_AAS(:,i), y_AAS,'LineWidth',2)
    hold on
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
end

