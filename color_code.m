function [ color_chart ] = color_code()
%Colors 

color_chart = ones(10,3);

color_chart(1,:) = [0,0,0];

color_chart(2,:) = [1,0,1];

color_chart(3,:) = [0,102,0]/255;
color_chart(4,:) = [0,204,0]/255;

color_chart(5,:) = [0,0,102]/255;
color_chart(6,:) = [0,0,204]/255;
color_chart(7,:) = [102,102,255]/255;

color_chart(8,:) = [153,0,0]/255;
color_chart(9,:) = [255,0,0]/255;
color_chart(10,:) = [255,153,153]/255;


end

