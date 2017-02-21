function [ ] = plotZrealData()
% Draw error bar chart with means and standard deviations
% for both node
path = '/Users/Jalaj/Documents/Github - LBPSS/Real_Data/run3/';
fileName = [path,'Z_error', '.mat'];
load(fileName);
data = Z_error_mean;

figure(1)
FigHandle = figure(1);
set(FigHandle, 'PaperPosition', [0 0 8.2 6]);
hold on
color_chart = color_code();
ind = [1 2 3 4];
for kk = 1:length(ind)
    k = ind(kk);
    %e1 = errorbar(k,mean_err_pw(k),std_err_pw(k),'o');
    e1(kk) = plot([k; k], [data(kk,1)-data(kk,2); data(kk,1)+data(kk,2)], '-');
    e0(kk) = plot(k,data(kk,1),'o');
    set(e0(kk),'Color', color_chart(kk,:))
    set(e0(kk),'MarkerFaceColor',color_chart(kk,:))
    set(e0(kk),'MarkerEdgeColor',[0 0 0])
    set(e0(kk),'MarkerSize',8)
    set(e1(kk),'Color', color_chart(kk,:))
    set(e1(kk),'LineWidth',2)
end


title('log Z ratio');
ylabel('Mean Error ')
box on
legend_lab = {'HMC','CMH','AAS + RB','AAG + RB'};
set(gca,'XTick',[])
set(gca,'FontSize', 18)
ylim([0 0.1+max(data(:,1))])
xlim([0 5])
% legend(e0,legend_lab,'Orientation','horizontal','Location','southoutside');
[legend_h,~,~,~] = columnlegend(5, e0, legend_lab, 2,'southoutside');
set(legend_h, 'position', [0.11 -0.20 1.05 0.35]);

name = 'Z_error';
saveas(1, [path, name], 'epsc')
% saveas(1, [path, name], 'png')
close all

end

