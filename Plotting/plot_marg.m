function [ ] = plot_marg( mean_err_pw, std_err_pw, title_pass, temp, scale_bias, scale_corr)
% Draw error bar chart with means and standard deviations
% for both node

figure(2)
FigHandle = figure(2);
set(FigHandle, 'PaperPosition', [0 0 7 5]);
hold on
color_chart = color_code_marg();
ind = [1 3 5 7 8 10 11 12 14 15 16];
for kk = 1:length(ind)
    k = ind(kk);
    %e1 = errorbar(k,mean_err_pw(k),std_err_pw(k),'o');
    e1(kk) = plot([k; k], [mean_err_pw(kk)-std_err_pw(kk); mean_err_pw(kk)+std_err_pw(kk)], '-');
    e0(kk) = plot(k,mean_err_pw(kk),'o');
    set(e0(kk),'Color', color_chart(kk,:))
    set(e0(kk),'MarkerFaceColor',color_chart(kk,:))
    set(e0(kk),'MarkerEdgeColor',[0 0 0])
    set(e0(kk),'MarkerSize',5)
    set(e1(kk),'Color', color_chart(kk,:))
    set(e1(kk),'LineWidth',2)
end


title(title_pass);
label_y = 'RMSE';
ylabel(label_y)
box on
% Change the labels for the tick marks on the x-axis
dum_label = '';
xtl1 = '\begin{tabular}{c} MH-LBP \end{tabular}';
xtl2 = '\begin{tabular}{c} HMC \end{tabular}';
xtl3 = '\begin{tabular}{c} CMH \end{tabular}';
xtl4 = '\begin{tabular}{c} CMH \\+\\LBP\end{tabular}';
xtl5 = '\begin{tabular}{c} AAS \end{tabular}';
xtl6 = '\begin{tabular}{c} AAS \\+\\RB\end{tabular}';
xtl7 = '\begin{tabular}{c} AAS \\ + \\ RB\\ + \\ LBP\end{tabular}';
xtl8 = '\begin{tabular}{c} AAG \end{tabular}';
xtl9 = '\begin{tabular}{c} AAG\\+\\ RB \end{tabular}';
xtl10 = '\begin{tabular}{c} AAG \\ + \\ RB\\ + \\ LBP\end{tabular}';
Algorithms = {xtl1, xtl2, xtl3, xtl4, xtl5, xtl6, xtl7, xtl8, xtl9, xtl10};
legend_lab = {'LBP', 'MH-LBP','HMC','CMH','CMH + LBP','AAS','AAS+RB','AAS+RB+LBP','AAG','AAG+RB','AAG+RB+LBP'};
%set(gca, 'XTick', 1:11, 'XTickLabel', Algorithms, 'TickLabelInterpreter', 'latex');
set(gca,'XTick',[])
ylim([0 0.1+max(mean_err_pw)])
xlim([0 17])
% legend(e0,legend_lab,'Orientation','horizontal','Location','southoutside');
[legend_h,~,~,~] = columnlegend(6, e0, legend_lab, 2,'southoutside');
set(legend_h, 'position', [-0.005 -0.22 1.05 0.35]);
name = strcat('Node_marg','','Temp', num2str(round(temp)), 'Bias', num2str(scale_bias*100), 'Corr', num2str(scale_corr*100));
path = '/Users/Jalaj/Documents/Github - LBPSS/Outputs_after_NIPS/Nodemarginal/';
% saveas(1, [path, name], 'epsc')
saveas(2, [path, name], 'png')
close all

end
