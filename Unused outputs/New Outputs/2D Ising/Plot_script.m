
scale = 1;
width = 6.25*scale;
height = 5.25*scale;

set(findall(gcf,'-property','FontSize'),'FontSize', 24) 
set(findall(gcf,'Type','Line'),'LineWidth',2)
set(findall(gcf,'Type','Line'), 'LineStyle','-')
ax=gca;
ax.XTickLabel = {'0','1250','2500','3750','5000','2500','3000','3500','4000','4500','5000'};
set(gcf, 'PaperSize', [width height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 width height]);

% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperSize', [width height]);
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperPosition', [0 0 width height]);

set(gcf, 'renderer', 'painters');
print(gcf, '-dpdf', 'fig14.pdf');
%print(gcf, '-dpng', 'my-figure.png');
%print(gcf, '-depsc2', 'my-figure.eps');