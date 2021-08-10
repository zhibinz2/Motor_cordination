function [x,y] = drawU1(steplength)
% Bimanual type U1
%   Detailed explanation goes here
x = [ones(1,steplength) [2:1:steplength-1] steplength*ones(1,steplength)];
y = [[1:1:steplength] steplength*ones(1,steplength-2) [steplength:-1:1]];
end

%% test
% steplength=200;
% [x,y] = drawU1(steplength);
% plot(x,y,'ro');
% set(gca, 'YDir', 'reverse');