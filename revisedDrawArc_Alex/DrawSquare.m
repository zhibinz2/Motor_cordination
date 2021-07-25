function [x,y] = DrawSquare(steplength)
x = [ones(1,steplength) [2:steplength-1] repmat(steplength,1,steplength)];
y = [[steplength:-1:1] ones(1,steplength-1) [(1+1):steplength]];
end


% plot(x,y);
% set(gca, 'YDir', 'reverse');