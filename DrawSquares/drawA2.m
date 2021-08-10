function [x,y] = drawA2(steplength)
% Bimanual type A2
%   Detailed explanation goes here
x = [steplength*ones(1,steplength) [steplength-1:-1:2] ones(1,steplength)];
y = [[steplength:-1:1] ones(1,steplength-2) [1:1:steplength]];
end

%% test
% steplength=200;
% [x,y] = drawA2(steplength);
% plot(x,y,'ro');
% set(gca, 'YDir', 'reverse');