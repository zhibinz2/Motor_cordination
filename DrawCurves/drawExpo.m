function [x,y] = drawExpo(radius, xCenter, yCenter);
% Bimanual type A1
% radius = maximum radius
% N = how many pixiels
% rad_ang = phase, or angle in radian
% r_radius = range of the radius

DirectDistance=sqrt((radius.^2)/2);
N      = radius;
x=linspace(0,DirectDistance,N) ;
ROOTexp=nthroot(DirectDistance+1,DirectDistance);
y=ROOTexp.^x-1;plot(x,y,'r');

y=flip(y);

x=x+xCenter;
y=y+yCenter-max(y);
end


%%
% close;figure;
% radius=320; rad_ang=pi/6; xCenter=1920/2; yCenter=1080/2;
% [x,y] = drawExpo(radius, xCenter, yCenter);
% plot(x,y,'go');
% xlim([0 1920]);ylim([0 1080]);
% hold on; 
% plot(xCenter, yCenter,'ro');
% plot(x(end),y(end),'ro');
% set(gca, 'YDir', 'reverse');

%%
% DirectDistance=sqrt((radius.^2)/2);
% close;figure;N=radius;
% xx=linspace(0,DirectDistance,N) ;
% plot(xx,xx,'k--');
% hold on;
% xline(0,'k--');yline(0,'k--');
% xline(DirectDistance/2,'k--');yline(DirectDistance/2,'k--');
% xlim([-radius radius]);ylim([-radius radius]);
% 
% ROOTexp=nthroot(DirectDistance+1,DirectDistance);
% yy=ROOTexp.^xx-1;plot(xx,yy,'r');
% % examine
% yy=ROOTexp.^DirectDistance-1
% 
% % example: log9(x) = log(x) / log(9)
% yy=log(xx+1)/log(ROOTexp);plot(xx,yy,'b');
% % examine
% yy=log(DirectDistance+1)/log(ROOTexp)


%%
% conditions = [0 2*atan(1/4) 2*atan(1/2) 2*atan(1)];
% 
% for i=1:length(conditions)
%     
%     rad_ang=conditions(i); 
% 
%     reach  = @(radius,rad_ang)  [radius*cos(rad_ang);  radius*sin(rad_ang)]; % Reach Trace Function For radius in pixels
%     N      = radius;
%     r_radius = linspace(0,radius,N);
%     xy_r   = reach(r_radius,rad_ang);
%     x      = xy_r(1,:);
%     y      = xy_r(2,:);
%     plot(x,y,'g');
% end

%%
% [x,y] = drawReach(radius, rad_ang, xCenter, yCenter);
% plot(x,y,'go');
% xlim([0 1920]);ylim([0 1080]);
% hold on; 
% plot(xCenter, yCenter,'ro');
% plot(x(end),y(end),'ro');
% set(gca, 'YDir', 'reverse');

%% plot different log and expo
% close;figure;
% xx=linspace(0,DirectDistance,N) ;
% plot(xx,xx,'k--');
% hold on;
% xline(0,'k--');yline(0,'k--');
% xline(DirectDistance/2,'k--');yline(DirectDistance/2,'k--');
% xlim([-radius radius]);ylim([-radius radius]);
% for ROOT=0.1:0.1:1
%     yy=ROOT.^xx-1;plot(xx,yy,'r');hold on;
% end
% for ROOT=1:0.1:3
%     yy=log(xx+1)/log(ROOT);plot(xx,yy,'b');hold on
% end
