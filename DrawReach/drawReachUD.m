function [xU,yU,xD,yD] = drawReachUD(radius,rad_ang, xCenter, yCenter);
% Bimanual type A1
% radius = maximum radius
% N = how many pixiels
% rad_ang = phase, or angle in radian
% r_radius = range of the radius

reach  = @(radius,rad_ang)  [radius*cos(rad_ang);  radius*sin(rad_ang)]; % Reach Trace Function For radius in pixels
N      = radius;
r_radius = linspace(0,radius,N);
xy_r   = reach(r_radius,rad_ang);
xU      = xy_r(1,:);
yU      = xy_r(2,:);
xD      =-xU;
yD      =-yU;

% flip up side down for psychtoolbox coordinate system
yU=-yU;
yD=-yD;

% shift to the Screen Center
xU=xU+xCenter;
yU=yU+yCenter;

xD=xD+xCenter;
yD=yD+yCenter;

end


%% test

% close;figure;
% radius=320; rad_ang=pi/6; xCenter=1920/2; yCenter=1080/2;
% [xU,yU,xD,yD] = drawReachUD(radius, rad_ang, xCenter, yCenter);
% plot(xU,yU,'mo');hold on;plot(xD,yD,'bo')
% xlim([0 1920]);ylim([0 1080]);
% plot(xCenter, yCenter,'go');
% plot(xU(end),yU(end),'go');
% plot(xD(end),yD(end),'go');
% set(gca, 'YDir', 'reverse');