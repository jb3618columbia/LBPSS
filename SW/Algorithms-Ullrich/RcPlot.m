function RcPlot(omega, title2)
% RCPLOT(SIGMA, TITLE) plots the Random cluster state OMEGA and sets the 
%   title of the figure to TITLE.

[row, col] = size(omega);   % size of connection matrix
row2 = sqrt(row);           % rows in graph
col2 = sqrt(col);           % cols in graph

%% Control the input arguments
if nargin<2
    title2 = '';
elseif (~ischar(title2))
    error('Second argument must be a string!');
end


%% Initialize plot
clf reset      % close current figure
plot_title = sprintf('%dx%d Random-cluster model',sqrt(row),sqrt(col));
title(plot_title)

set(gca,'YTickLabel',[],'XTickLabel',[],'XTick',[],'YTick',[]); 
axis ij; axis([0.99 row2+.5 0.99 col2+.5]); axis square; 
set(gca,'XColor','w','YColor','w')
colormap bone; drawnow;

xlabel(title2,'Color','k')

%% Plot the Random cluster state
[v1_all v2_all] = find(omega);    % edges

% % without periodic boundary condition
% v1 = mod(E-1,row)+1; % definition of connected edges
% v2 = ceil(E/row);
% c1y = sqrt(row)-mod(v1-1,sqrt(row));    % definition of the coordinates
% c1x = ceil(v1/sqrt(row));     %   (origin in the upper left)
% c2y = sqrt(row)-mod(v2-1,sqrt(row));
% c2x = ceil(v2/sqrt(row));
% 
% coordx=[c1x'; c2x'];    % preparing for plot
% coordy=[c1y'; c2y'];
% 
% line(coordx,coordy,'Color','k')


% with periodic boundary condition
for i = 1:length(v1_all)
%for i =1:length(E)
%     v1 = mod(E(i)-1,row)+1; % definition of connected edges
%     v2 = ceil(E(i)/row);
    
    v1 = v1_all(i);
    v2 = v2_all(i);

    c1y = mod(v1-1,row2)+1;      % definition of the coordinates
    c1x = ceil(v1/row2);            %   (origin in the upper left)
    c2y = mod(v2-1,row2)+1;
    c2x = ceil(v2/row2);
    
    
    % right outer edge
    if ( (c1x==row2 && c2x==1) || (c1x==1 && c2x==row2) ) 
        
        line([row2 row2+1],[c1y c1y],'Color','k');
        
    % lower outer edge
    elseif ( (c1y==col2 && c2y==1) || (c1y==1 && c2y==col2) )
        
        line([c1x c1x],[col2 col2+1],'Color','k');
        
    % inner edges
    else 
        line([c1x c2x],[c1y c2y],'Color','k');
    end

end



