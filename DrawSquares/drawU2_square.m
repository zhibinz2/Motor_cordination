function [x,y] = drawU2(steplength)
% Bimanual type U2
%   Detailed explanation goes here
x = [steplength*ones(1,steplength) [steplength-1:-1:2] ones(1,steplength)];
y = [[1:1:steplength] steplength*ones(1,steplength-2) [steplength:-1:1]];
end

%% test
% steplength=200;
% [x,y] = drawU2(steplength);
% plot(x,y,'ro');
% set(gca, 'YDir', 'reverse');