function [poly0] = polycircle(x,y,r)
% Use polyshape to construct a circle

th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;

poly0 = polyshape(xunit(1:end-1),yunit(1:end-1));
end

