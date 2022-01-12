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

y=flip(y);

x=x+xCenter;
y=y+yCenter-max(y);



end


%% test
radius=320; rad_ang=pi/6; 
[x,y] = drawReach(radius, rad_ang, xCenter, yCenter);
plot(x,y,'go');
xlim([0 1920]);ylim([0 1080]);
hold on; 
plot(xCenter, yCenter,'ro');
plot(x(end),y(end),'ro');
set(gca, 'YDir', 'reverse');