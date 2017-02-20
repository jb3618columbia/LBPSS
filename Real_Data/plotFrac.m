function [ ] = plotFrac()

path = '/Users/Jalaj/Documents/Github - LBPSS/Real_Data/fake_data_highdim/run2/';
lab = char(['HMC', 'CMH', 'AAS', 'AAG']);
for j=1:length(lab)/3
    fileName = sprintf('%s%s.mat', path, lab((j-1)*3 + 1:j*3));
    load(fileName);
end

filename = [path, 'W', '.mat'];
load(filename);

percent=0.2;
frac_hmc = giveFrac(samples_HMC, W, percent);
frac_cmh = giveFrac(samples_CMH, W, percent);
frac_aas = giveFrac(samples_AAS, W, percent);
frac_aag = giveFrac(samples_AAG, W, percent);

% fileName = [path, 'frac_hmc','.mat'];
% load(fileName);
% fileName1 = [path, 'frac_cmh','.mat'];
% load(fileName1);
% fileName2 = [path, 'frac_aas','.mat'];
% load(fileName2);
% fileName3 = [path, 'frac_aag','.mat'];
% load(fileName3);


figure(1)
FigHandle = figure(1);
set(FigHandle, 'PaperPosition', [0 0 7 5]);
plot(frac_hmc, 'LineStyle', '--')
hold on
plot(frac_cmh, 'LineStyle', '--')
hold on
plot(frac_aas, 'LineStyle', '--');
hold on
plot(frac_aag, 'LineStyle', '--');
set(gca,'Position',[0.16 0.15 0.73 0.77])
set(gca,'FontSize', 18)

title('Approx. MH sampler');
xlabel('Parameters');
ylabel('fraction (f)');

lgd = legend('HMC','CMH','AAS','AAG');
lgd.Orientation = 'vertical';
% lgd.Location = 'north';
% lgd.Position = [0.018,-0.45,1,1];

name = 'frac_plot';
saveas(1, [path, name], 'epsc')
saveas(1, [path, name], 'png')
close all
