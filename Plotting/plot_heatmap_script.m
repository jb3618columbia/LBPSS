function [ ] = plot_heatmap_script()

on_off = 1;  % To remove the colorbar and yaxis labels
dataPath = '/Users/Jalaj/Documents/Github - LBPSS/Outputs_AISTATS/manipulated_data/fig_2/';
if on_off == 1
    plotsPath = '/Users/Jalaj/Documents/Github - LBPSS/Outputs_AISTATS/paper_plots_AISTATS/';
else
    plotsPath = '/Users/Jalaj/Documents/Github - LBPSS/Outputs_after_NIPS/plots2/';
end
% scale_vec = linspace(0,0.15,4);
% scale_conn = [linspace(0.2,0.5,5) , linspace(0.55,0.75,10)]; 
% for figure 1
scale_vec = 0.2;
scale_conn = [0.1:0.1:0.8];
% for figure 2
% scale_vec = 0
% scale_conn = [0.25:0.25:1, 2, 3, 4, 5]
% for figure 3
% scale_vec = 0.2;
% scale_conn = [0.2:0.2:1, 1.5, 2, 3, 4];

% Passing the entire connection strength, 
data_code = [{'err_pw'}];
% data_code = [{'err'}, {'err_pw'}];
% data_code = [{'err'}, {'err_pw'}, {'act'}];

legendText = {'LBP','HMC','CMH','CMH + LBP','AAS','AAS+RB','AAS+RB+LBP','AAG','AAG+RB','AAG+RB+ST', 'AAG+RB+LBP'};
for j = 1:length(data_code)
    for bias_ind = 1:length(scale_vec)
        plot_heatmap(bias_ind, scale_conn, data_code{j}, legendText, plotsPath, dataPath, on_off)
    end
end


end

