function [ ] = plotRealdata()
% This is to plot the histigrams for samples obtained in outer MH loop.
% Not including AAS as we have too less samples here.

path = '/Users/Jalaj/Documents/Github - LBPSS/Real_Data/fake_data_highdim/run1/';
lab = char(['HMC', 'CMH', 'AAS', 'AAG']);
% lab = char(['tru1','HMC1', 'CMH1', 'AAS1', 'AAG1']);

figure
for i=1:15
    subplot(3,5,i);
    for j=1:length(lab)/3
        fileName = sprintf('%s%s.mat', path, lab((j-1)*3 + 1:j*3));
        load(fileName);
    end
    i=i+30;
    filename = [path, 'W.mat'];
    load(filename);
    W = nonzeros(W)';
    x_grid = -4:0.1:4;
%     pd_true = fitdist(samples_true1(:,i),'Normal');
%     y_true = pdf(pd_true, x_grid);
    true = W(1, i);
    pd_HMC = fitdist(samples_HMC(:,i),'Normal');
    y_HMC = pdf(pd_HMC, x_grid);
    pd_CMH = fitdist(samples_CMH(:,i),'Normal');
    y_CMH = pdf(pd_CMH, x_grid);
    pd_AAS = fitdist(samples_AAS(:,i),'Normal');
    y_AAS = pdf(pd_AAS, x_grid);
    pd_AAG = fitdist(samples_AAG(:,i),'Normal');
    y_AAG = pdf(pd_AAG, x_grid);
    
    
    %     pd_true = fitdist(samples_true1(:,i),'Normal');
    %     y_true = pdf(pd_true, samples_true1(:,i));
    %     pd_HMC = fitdist(samples_HMC1(:,i),'Normal');
    %     y_HMC = pdf(pd_HMC, samples_HMC1(:,i));
    %     pd_CMH = fitdist(samples_CMH1(:,i),'Normal');
    %     y_CMH = pdf(pd_CMH, samples_CMH1(:,i));
    %     pd_AAS = fitdist(samples_AAS1(:,i),'Normal');
    %     y_AAS = pdf(pd_AAS, samples_AAS1(:,i));
    %     pd_AAG = fitdist(samples_AAG1(:,i),'Normal');
    %     y_AAG = pdf(pd_AAG, samples_AAG1(:,i));
%     plot(x_grid, y_true,'LineWidth',.4)
%     hold on
    gt = [true true];
    y = [0 0.5];
    plot(gt, y);
    hold on 
    plot(x_grid, y_HMC,'LineWidth',2)
    hold on
    plot(x_grid, y_CMH,'Color','green', 'LineWidth',.4)
    hold on
    plot(x_grid, y_AAS,'Color','red','Linestyle', '--','LineWidth',.4)
    ylim([0, 0.6])
    if i >=2
        set(gca,'YTickLabel',[])
    end
    %     hold on
    %     plot(x_grid, y_AAG,'LineWidth',2)
    
    % histogram(samples_true1(:,i))
    % hold on
    % histogram(samples_HMC1(:,i))
    % hold on
    % histogram(samples_CMH1(:,i))
    % hold on
    % histogram(samples_AAS1(:,i))
    % hold on
    % histogram(samples_AAG1(:,i))
    
end
lgd = legend('True','CMH','AAS');
lgd.Orientation = 'horizontal';
lgd.Location = 'south';
lgd.Position = [0.018,-0.45,1,1];
name = 'fit_samples';
saveas(1, [path, name], 'epsc')
saveas(1, [path, name], 'png')
close all

end

