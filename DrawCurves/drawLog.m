function [x,y] = drawLog(radius, xCenter, yCenter);
% Bimanual type A1
% radius = maximum radius
% N = how many pixiels
% rad_ang = phase, or angle in radian
% r_radius = range of the radius

DirectDistance=sqrt((radius.^2)/2);
N      = radius;
x=linspace(0,DirectDistance,N) ;
ROOTexp=nthroot(DirectDistance+1,DirectDistance);
y=y=log(x+1)/log(ROOTexp);


y=flip(y);

x=x+xCenter;
y=y+yCenter-max(y);
end


%%
% close;figure;
% radius=320; rad_ang=pi/6; xCenter=1920/2; yCenter=1080/2;
% [x,y] = drawLog(radius, xCenter, yCenter);
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
