function [ ] = plot_fn( mean_err_pw, std_err_pw, title_pass, act, temp, scale_bias, scale_corr)
% Draw error bar chart with means and standard deviations
% for both node

figure(1)
FigHandle = figure(1);
set(FigHandle, 'PaperPosition', [0 0 7 5]);
hold on
color_chart = color_code();
ind = [1 3 5 6 8 9 10 12 13 14];
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
if act==1
    label_y = 'ACT';
else
    label_y = 'RMSE';
end
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
legend_lab = {'MH-LBP','HMC','CMH','CMH + LBP','AAS','AAS+RB','AAS+RB+LBP','AAG','AAG+RB','AAG+RB+LBP'};
%set(gca, 'XTick', 1:11, 'XTickLabel', Algorithms, 'TickLabelInterpreter', 'latex');
set(gca,'XTick',[])
ylim([0 0.1+max(mean_err_pw)])
xlim([0 15])
% legend(e0,legend_lab,'Orientation','horizontal','Location','southoutside');
[legend_h,~,~,~] = columnlegend(5, e0, legend_lab, 2,'southoutside');
set(legend_h, 'position', [0.045 -0.20 0.95 0.35]);
if act ==1
    name = strcat('Act_mag,'',''Temp', num2str(round(temp)), 'Bias', num2str(scale_bias), 'Corr', num2str(scale_corr));
    path = 'C:\Users\Student.DESKTOP-GMAHVHB\Downloads\LBPSS-master\LBPSS-master\Outputs_after_NIPS\Actmag\';
else
    name = strcat('Pair_marg','','Temp', num2str(round(temp)), 'Bias', num2str(scale_bias), 'Corr', num2str(scale_corr));
    path = 'C:\Users\Student.DESKTOP-GMAHVHB\Downloads\LBPSS-master\LBPSS-master\Outputs_after_NIPS\Pairmarginal\';
end
% saveas(1, [path, name], 'epsc')
saveas(1, [path, name], 'png')
close all

% Save as Matlab figs
% name = strcat('Temp', num2str(temp), 'Bias', num2str(scale_bias), '.fig');
% path = '/Users/Jalaj/Documents/Github - LBPSS/Outputs_after_NIPS/Pairwise_marginals';
% savefig(gcf, fullfile(path, name))


end

