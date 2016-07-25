function IsingPlot(sigma, title2)

% ISINGPLOT(SIGMA, TITLE) plots the Ising state SIGMA and sets the 
%   title of the figure to TITLE.

%% Control the input arguments
if nargin<2
    title2 = '';
elseif (~ischar(title2))
    error('Second argument must be a string!');
end


%% Plot the Ising state

[row, col] = size(sigma);
plot_title = sprintf('%dx%d Ising model',row,col);

image((sigma+1)*128);

title(plot_title)
xlabel(title2,'Color','k')
set(gca,'YTickLabel',[],'XTickLabel',[],'XTick',[],'YTick',[]); 
axis square; colormap bone; drawnow;

