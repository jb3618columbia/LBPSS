function [data_out] = load_data(bias_ind, scale_conn, data_code, dataPath)

for w = 1:length(scale_conn)
%     path = '/Users/Jalaj/Documents/Github - LBPSS/Outputs_after_NIPS/data/';
    path = dataPath;
    fileName = [path, 'data_', data_code ,'_', num2str(bias_ind), num2str(w), '.mat'];
    load(fileName);
    data_out(:,:,w) = mean_std;
end
data_out
end

