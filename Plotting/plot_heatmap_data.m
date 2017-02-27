function [ ] = plot_heatmap_data(data, legendText, scale_conn, plotsPath, data_code, bias_ind, on_off)

nAlg = size(data,1); % for err_pw
% nAlg = size(data,1)-1; % for err to remove the result for LBP
nStrength = size(data,3);

% Manipulate the data matrix construction here to cherry pick the values
% for connection strength
cherry_pick = 0;
if cherry_pick == 1
    scale_conn_new = [2,4,6,8,10,15]; % to cherry pick values
    scale_conn_new = round(scale_conn_new,2);
    arr = [2, 3, 4, 5, 6, 7];   % the corresponding indeces hard coded.
    heatmapMatrix = zeros(nAlg, length(scale_conn_new));
    for l=1:length(scale_conn_new)
        heatmapMatrix(:,l) = data(:,1,arr(l));
    end
    scale_conn = scale_conn_new;
else
    heatmapMatrix = zeros(nAlg, nStrength);
    for k = 1:nStrength
%         heatmapMatrix(:,k) = data(2:nAlg+1,1,k); % for err
        heatmapMatrix(:,k) = data(:,1,k); % for err_pw
    end
end


figure(1)
% colormap('hot')
colormap(autumn(50000));
brighten(0.6);
imagesc(heatmapMatrix);
xlabel('Strength (W)');
% xlabel('Bias scale (c)'); 

scale_conn = round(scale_conn,3);
if on_off ==1
    if strcmp(data_code, 'err')
        %         ylabel('Algorithm')
        title('RMSE(Node Marginals)')
        set(gca,'Position',[0.3 0.13 0.65 0.77])
    else
        title('RMSE(Pairwise Marginals)')
    end
    
    %     set(gca,'XTick',[1:nStrength], 'XTickLabel', strread(num2str(scale_conn),'%s'))
    set(gca,'XTick',[1:nStrength], 'XTickLabel', strread(num2str(scale_conn),'%s'))
    set(gca, 'YTickLabel', legendText)
    set(gca,'Position',[0.23 0.13 0.73 0.77])
    set(gca,'FontSize', 18)
    c = colorbar;
    minVal = round(0,3);
    maxVal = round(max(max(heatmapMatrix)),2);
    nLabels = 5;
    c.Limits = [minVal, maxVal];
    cTickLabels = round(linspace(minVal,maxVal,nLabels),2);
    c.Ticks = cTickLabels;
else
    %     set(gca,'XTick',[1:nStrength], 'XTickLabel', strread(num2str(scale_conn),'%s'))
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

