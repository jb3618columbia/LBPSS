function [data_out] = load_data(bias_ind, scale_conn, data_code)

for w = 1:length(scale_conn)
    path = '/Users/Jalaj/Documents/Github - LBPSS/Outputs_after_NIPS/data/';
        fileName = [path, 'data_', data_code ,'_', num2str(bias_ind), num2str(w), '.mat'];
        load(fileName);
        data_out(:,:,w) = mean_std;
end

end

