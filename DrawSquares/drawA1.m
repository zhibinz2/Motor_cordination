function [x,y] = drawA1(steplength)
% Bimanual type A1
%   Detailed explanation goes here
x = [ones(1,steplength) [2:1:steplength-1] steplength*ones(1,steplength)];
y = [[steplength:-1:1] ones(1,steplength-2) [1:1:steplength]];
end

%% test
% steplength=200;
% [x,y] = drawA1(steplength);
% plot(x,y,'ro');
% set(gca, 'YDir', 'reverse');