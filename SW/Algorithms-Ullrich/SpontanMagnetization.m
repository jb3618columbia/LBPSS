function result = SpontanMagnetization(steps, start, size)
% RESULT = SPONTANMAGNETIZATION(STEPS, START, SIZE) approximates the 
%   magnetization of an SIZExSIZE Ising model with start determined by
%   START and STEPS steps of the sweep heat bath algorithm algorithm.
%   e.g. result = SpontanMagnetization(1000,0);
%	RESULT is a vector with the results of the approximation.
%   START determine the start values 
%       (0 for given PW data (DEFAULT; data ...
%				must be in ".../SpontanMagnetization/initialStates/")
%        1 for deterministic ones 
%       -1 for random uniform distributed states)
%	The exactly distributed states are given by a run of the 
%   Propp-Wilson algorithm.
%   SIZE is lattice size of the model
%
%	The vector RESULT will be saved every 1000000 steps in the directory 
%	".../SpontanMagnetization/" as the file "SM-'SIZE'-'STEP'.mat".


%% Control the input arguments
if nargin<2
    start = 0;
    size = 250;
elseif nargin<3
    size = 250;
elseif nargin>3
    error('Too many input arguments!');
elseif nargin==0
    error('No input arguments!');
end

%% Preliminaries

format long;

seed_start = 12345;
RandStream.setDefaultStream(RandStream('mt19937ar','seed',seed_start));

%% Define the x-axis and variables
name = {'1.763' '1.322' '1.146'... 
        '1.058' '1.014' '0.970'... 
        '0.925' '0.908' '0.899'... 
        '0.890' '0.881' '0.877'... 
        '0.873' '0.868' '0.864'... 
        '0.855' '0.837' '0.793'... 
        '0.749' '0.705' '0.661'... 
        '0.573' '0.5'};
beta = [2*log(1+sqrt(2)) 1.5*log(1+sqrt(2)) 1.3*log(1+sqrt(2))... 
        1.2*log(1+sqrt(2)) 1.15*log(1+sqrt(2)) 1.1*log(1+sqrt(2))... 
        1.05*log(1+sqrt(2)) 1.03*log(1+sqrt(2)) 1.02*log(1+sqrt(2))... 
        1.01*log(1+sqrt(2)) log(1+sqrt(2)) .995*log(1+sqrt(2))... 
        .99*log(1+sqrt(2)) .985*log(1+sqrt(2)) .98*log(1+sqrt(2))... 
        .97*log(1+sqrt(2)) .95*log(1+sqrt(2)) .9*log(1+sqrt(2))... 
        .85*log(1+sqrt(2)) .8*log(1+sqrt(2)) 0.75*log(1+sqrt(2))... 
        0.65*log(1+sqrt(2)) .5];

%omega = sparse([],[],[],250^2,250^2,2*250^2); % empty set

state = zeros(size,size,length(name));
result = zeros(1,length(name)); % absolute magnetization

%% Load results from the Propp-Wilson algorithm
% state(:,:,1) = load('SpontanMagnetization/PW-250_0.970.mat');
% state(:,:,2) = load('SpontanMagnetization/PW-250_0.88.mat');
% state(:,:,3) = load('SpontanMagnetization/PW-250_0.837.mat');
% state(:,:,4) = load('SpontanMagnetization/PW-250_0.793.mat');
% state(:,:,5) = load('SpontanMagnetization/PW-250_0.661.mat');
% state(:,:,6) = load('SpontanMagnetization/PW-250_0.5.mat');
if start == -1
    for i=1:length(name)
        state(:,:,i) = (-1).^(round(rand(size)));
        result(i) = abs(sum(sum(state(:,:,i))))/(size^2);
    end
elseif start == 1
    for i=1:length(name)
        state(:,:,i) = ones(size);
        result(i) = abs(sum(sum(state(:,:,i))))/(size^2);
    end
else
    for i=1:length(name)
        state(:,:,i) = cell2mat(struct2cell(...
            load(strcat('SpontanMagnetization/initialStates/PW-',num2str(size),'_', name{i}, '.mat'))));
        result(i) = abs(sum(sum(state(:,:,i))))/(size^2);
    end
end
%% Create plot
% exact solution (N=\infty)
f = @(t) (1-(sinh(1./t)).^(-4)).^(1/8).*(t<1/log(1+sqrt(2)));
fplot(f,[0 2],'r')
hold on

% simulation data
data = [1 result];
h = plot([0 1./beta], data,'-*','YDataSource','data');
axis([0 2 -.1 1.1])
title(strcat('Spontaneous Magnetization of a ',...
    num2str(size), '\times', num2str(size), ' lattice'))
hold off


%% Exact solution for comparison
exact = [1 f(1./beta)];

%% Approximation
for k=1:steps

% 	if getappdata(h,'canceling')
% 		break
% 	end

	% One step for each value of beta
	for l=1:length(name)
% 		% One step for beta(l) (Swendsen-Wang)
% 		omega = IsingToRc(state(:,:,l),beta(l));
% 		state(:,:,l) = RcToIsing(omega);
%         % Calculation of the new estimate (k+1 summands)
% 		temp = result(l);
% 		result(l) = (k*temp + sum(sum(state(:,:,l))))/(k+1);
		
        % One step for beta(l) (Heat bath)
        [state(:,:,l) M] = IsingHeatbathSweepStep(state(:,:,l),beta(l),-1);
        result(l) = (k*result(l) + abs(M)/(size^2))/(k+1);
        
    end
	
    % Update the plot
        data = [1 result];
        refreshdata(h,'caller')
        drawnow;
		% plot([0 1./beta 2], [1 result 0],'-*')
		% axis([0 2 -.1 1.1])
        % pause(0.00001)
    
    % Comparison (Mean square error on the sampling points)
    xlabel(sprintf('step = %d          MSE = %g',k,(sum((data-exact).^2))^(1/2)))
    
    % Save data
    if ~mod(k,1000000)
        save(strcat('SpontanMagnetization/SM-',num2str(size),'-',num2str(k),'.mat'),'result')
    end
    
end

%close(h)

end