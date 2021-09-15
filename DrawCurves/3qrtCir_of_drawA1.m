function [x,y] = drawA1(radius)
% Bimanual type A1
% radius, or steplenth
% N = how many pixiels
% rand_ang = phase, or angle
% r_rangl= range of the angle

circr  = @(Radius,rad_ang)  [radius*cos(rad_ang);  radius*sin(rad_ang)]; % Circle Function For Angles In Radians
N      = 3 * radius;
% r_angl = linspace(-0.25*pi,1.25*pi,N);
r_angl = linspace(0.75*pi,2.25*pi,N);
xy_r   = circr(radius,r_angl);
x      = xy_r(1,:);
x      = x-min(x)+1; % shift into the corner of the screen
% x      = flip(x);
y      = xy_r(2,:);
y      = y-min(y)+1; % shift into the corner of the screen
end

% function [x,y] = drawA1(steplength);
% % Bimanual type A1
% %   Detailed explanation goes here
% x = [ones(1,steplength) [2:1:steplength-1] steplength*ones(1,steplength)];
% y = [[steplength:-1:1] ones(1,steplength-2) [1:1:steplength]];
% end

%% test
% steplength=200; %radius=200
% [x,y] = drawA1(steplength);
% plot(x,y,'ro');
% hold on;
% plot(steplength,steplength,'go')
% set(gca, 'YDir', 'reverse');