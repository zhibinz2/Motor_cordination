function [x,y] = drawA1(radius);
% Bimanual type A1
% radius, or steplength
% N = how many pixiels
% rand_ang = phase, or angle
% r_rangl= range of the angle

circr  = @(Radius,rad_ang)  [radius*cos(rad_ang);  radius*sin(rad_ang)]; % Circle Function For Angles In Radians
N      = 3 * radius;
r_angl = linspace(pi,2*pi,N);
xy_r   = circr(radius,r_angl);
x      = xy_r(1,:);
x      = x-min(x)+1; % shift into the corner of the screen
% x      = flip(x);
y      = xy_r(2,:);
y      = y-min(y)+1; % shift into the corner of the screen
end


%% test
% steplength=200; %radius=320
% [x,y] = drawA1(steplength);
% plot(x,y,'ro');
% hold on; plot(steplength,steplength,'go');
% plot(steplength,steplength,'go');plot(x(1),y(1),'bo');
% set(gca, 'YDir', 'reverse');