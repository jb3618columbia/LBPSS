function [ ] = plot_heatmap_script()

on_off = 1;  % To remove the colorbar and yaxis labels
if on_off == 1
    plotsPath = '/Users/Jalaj/Documents/Github - LBPSS/Outputs_after_NIPS/plots/';
else
    plotsPath = '/Users/Jalaj/Documents/Github - LBPSS/Outputs_after_NIPS/plots2/';
end
scale_vec = linspace(0,0.15,4);
scale_conn = linspace(0.1,0.6,5);
data_code = [{'err'}, {'err_pw'}, {'act'}];
legendText = {'LBP','HMC','CMH','CMH + LBP','AAS','AAS+RB','AAS+RB+LBP','AAG','AAG+RB','AAG+RB+LBP'};
for j = 1:length(data_code)
    for bias_ind = 1:length(scale_vec)
        plot_heatmap(bias_ind, scale_conn, data_code{j}, legendText, plotsPath, on_off)
    end
end


end

