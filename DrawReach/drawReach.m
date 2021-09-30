function [x,y] = drawReach(radius,rad_ang, xCenter, yCenter);
% Bimanual type A1
% radius = maximum radius
% N = how many pixiels
% rad_ang = phase, or angle in radian
% r_radius = range of the radius

reach  = @(radius,rad_ang)  [radius*cos(rad_ang);  radius*sin(rad_ang)]; % Reach Trace Function For radius in pixels
N      = radius;
r_radius = linspace(0,radius,N);
xy_r   = reach(r_radius,rad_ang);
x      = xy_r(1,:);
y      = xy_r(2,:);


x=x+xCenter;
y=y+yCenter;

y=flip(y);

end


%% test
% radius=320; rad_ang=pi/6; 
% [x,y] = drawReach(radius, rad_ang, xCenter, yCenter);
% plot(x,y,'go');
% hold on; 
% plot(x(end),y(end),'ro');
% set(gca, 'YDir', 'reverse');