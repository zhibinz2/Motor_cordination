function [xxx, yyy] = keepvalues(x,y,lowlim, uplim)
% keep only coordinates between lowlim and uplim in the data
% x and y are of the same length
xxx=x(logical((x>lowlim & x<uplim).*(y>lowlim & y<uplim)));
yyy=y(logical((x>lowlim & x<uplim).*(y>lowlim & y<uplim)));
end
% example:
% lowlim=0.5;uplim=1;
% x=[0.1 0.2 0.3 0.5 0.6 0.7 0.8 1];
% y=[0.4 0.5 0.6 0.7 0.8 0.9 1 1.1];
% (x>lowlim & x<uplim).*(y>lowlim & y<uplim)