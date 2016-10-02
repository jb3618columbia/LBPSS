function [ ] = plot_heatmap_data(data, legendText, scale_conn, plotsPath, data_code, bias_ind, on_off)

nAlg = size(data,1);
nStrength = size(data,3);
heatmapMatrix = zeros(nAlg, nStrength);

for k = 1:nStrength
    heatmapMatrix(:,k) = data(:,1,k);
end

figure(1)
colormap('hot')
imagesc(heatmapMatrix);
xlabel('Strength');
if on_off ==1
    % ylabel('Algorithm')
    title('abc')
    set(gca,'XTick',[1:nStrength], 'XTickLabel', strread(num2str(scale_conn),'%s'))
    set(gca, 'YTickLabel', legendText)
    set(gca,'Position',[0.15 0.13 0.8 0.8])
    c = colorbar;
    minVal = round(0,3);
    maxVal = round(max(max(heatmapMatrix)),3);
    nLabels = 5;
    c.Limits = [minVal, maxVal];
    cTickLabels = linspace(minVal,maxVal,nLabels);
    c.Ticks = cTickLabels;
else
    set(gca,'XTick',[1:nStrength], 'XTickLabel', strread(num2str(scale_conn),'%s'))
    set(gca, 'YTick', [])
    set(gca,'Position',[0.05 0.1 0.9 0.85])
end
set(gcf,'renderer','opengl');
fileName = [plotsPath, 'plot_', data_code, '_', num2str(bias_ind)];
fileNamePNG = [plotsPath, 'p_plot_', data_code, '_', num2str(bias_ind)];
saveas(1, fileName, 'epsc')
saveas(1, fileNamePNG, 'png')
close

end

