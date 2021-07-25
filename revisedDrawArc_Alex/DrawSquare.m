function [x,y] = DrawSquare(length)
x = [ones(1,length) [2:length-1] repmat(length,1,length)];
y = [[length:-1:1] ones(1,length-1) [(1+1):length]];
end

      