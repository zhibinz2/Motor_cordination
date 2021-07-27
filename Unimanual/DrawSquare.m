function [x,y] = DrawSquare(steplength)
x = [ones(1,steplength) [2:steplength-1] repmat(steplength,1,steplength)];
y = [[steplength:-1:1] ones(1,steplength-1) [(1+1):steplength]];
end

%% test
% steplength=200;
% [x,y] = DrawSquare(steplength)
% plot(x,y,'ro');
% set(gca, 'YDir', 'reverse');