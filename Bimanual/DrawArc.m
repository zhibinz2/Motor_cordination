function [x,y] = DrawArc(radius,N)
% N is the number of dots 
circr  = @(Radius,rad_ang)  [radius*cos(rad_ang);  radius*sin(rad_ang)]; % Circle Function For Angles In Radians
r_angl = linspace(pi/4,3*pi/4,N);
xy_r   = circr(radius,r_angl);
x      = xy_r(1,:);
y      = xy_r(2,:);

end

      