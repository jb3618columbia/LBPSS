function [ ] = plotZrealData()
% Draw error bar chart with means and standard deviations
% for both node
% path = '/Users/Jalaj/Documents/Github - LBPSS/Real_Data/data/run4/';
path = '/Users/Jalaj/Documents/Github - LBPSS/Real_Data/fake_data/run3/';
fileName = [path,'Z_error', '.mat'];
load(fileName);
data = Z_error_mean;

path = '/Users/Jalaj/Documents/Github - LBPSS/Real_Data/data/run4/';
fileName = [path,'Z_error', '.mat'];
load(fileName);
data1 = Z_error_mean;

figure(1)
FigHandle = figure(1);
set(FigHandle, 'PaperPosition', [0 0 7.5 4]);
hold on


color_chart = color_code();
ind = [1 2 3 4];
for kk = 1:length(ind)
    k = ind(kk);
    %e1 = errorbar(k,mean_err_pw(k),std_err_pw(k),'o');
    e1(kk) = plot([k; k], [data(kk,1)-data(kk,2); data(kk,1)+data(kk,2)], '-');
    e0(kk) = plotyy(k,data(kk,1),k, data1(kk,1));
    set(e0(kk),'Color', color_chart(kk,:))
    set(e0(kk),'MarkerFaceColor',color_chart(kk,:))
    set(e0(kk),'MarkerEdgeColor',[0 0 0])
    set(e0(kk),'MarkerSize',8)
    set(e1(kk),'Color', color_chart(kk,:))
    set(e1(kk),'LineWidth',2)
end

box on
set(gca,'XTick',[])
set(gca,'FontSize', 18)
ylim([1.5 0.4+max(data(:,1))])
xlim([0 5])

title('Approx. MH sampler');
ylabel('Approximation error')

legend_lab = {'HMC','CMH','AAS + RB','AAG + RB'};

% legend(e0,legend_lab,'Orientation','horizontal','Location','southoutside');
[legend_h,~,~,~] = columnlegend(5, e0, legend_lab, 2,'southoutside');
set(legend_h, 'position', [0.15 -0.25 0.95 0.35]);

name = 'Z_error';
saveas(1, [path, name], 'epsc')
saveas(1, [path, name], 'png')
close all

end

