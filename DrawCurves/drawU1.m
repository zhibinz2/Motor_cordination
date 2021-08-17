function [x,y] = drawU1(radius)
% Bimanual type A1
% radius, or steplenth
% N = how many pixiels
% rand_ang = phase, or angle
% r_rangl= range of the angle

circr  = @(Radius,rad_ang)  [radius*cos(rad_ang);  radius*sin(rad_ang)]; % Circle Function For Angles In Radians
N      = 3 * radius;
r_angl = linspace(-0.25*pi,1.25*pi,N);
% r_angl = linspace(0.75*pi,2.25*pi,N);
xy_r   = circr(radius,r_angl);
x      = xy_r(1,:);
x      = x-min(x)+1; % shift into the corner of the screen
x      = flip(x);
y      = xy_r(2,:);
y      = y-min(y)+1; % shift into the corner of the screen
end

%% test
% steplength=200; %radius=200
% [x,y] = drawU1(steplength);
% plot(x,y,'ro');
% set(gca, 'YDir', 'reverse');