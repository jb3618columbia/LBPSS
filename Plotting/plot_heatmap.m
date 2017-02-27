function [ ] = plot_heatmap(bias_ind, scale_conn, data_code, legendText, plotsPath, dataPath, on_off)

%3D. Third dimension is correlation. In each dimension, column 1 is mean
%and column 2 is stdev.
data = load_data(bias_ind, scale_conn, data_code, dataPath);

%plot
if strcmp(data_code, 'err')
    %do nothing
    legendText(1) = [];
else
    legendText(1) = [];
end
plot_heatmap_data(data, legendText, scale_conn, plotsPath, data_code, bias_ind, on_off)

end

