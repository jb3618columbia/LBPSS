function IsingSave(sigma, title);

% ISINGSAVE(SIGMA) saves the picture of the Ising state SIGMA. 
%   The name of the generated file will be TITLE.

% Control the input arguments
if nargin<2 
    [row, col] = size(sigma);
    time = round(clock);
    title = strcat('Ising_',num2str(row),'x',num2str(col),'_',date,'_',...
        num2str(time(4)),'-',num2str(time(5)),'-',num2str(time(6)))
elseif (~ischar(title))
    error('Second argument must be a string!');
end

% Write the data
imwrite((sigma+1)*128,strcat(title,'.jpg'),'jpg');