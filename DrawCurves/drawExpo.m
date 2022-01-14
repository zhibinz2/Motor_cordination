DirectDistance=sqrt((radius.^2)/2);
close;figure;N=radius;
xx=linspace(0,DirectDistance,N) ;
plot(xx,xx,'k--');
hold on;
xline(0,'k--');yline(0,'k--');
xline(DirectDistance/2,'k--');yline(DirectDistance/2,'k--');
xlim([-radius radius]);ylim([-radius radius]);

ROOTexp=nthroot(DirectDistance+1,DirectDistance);
yy=ROOTexp.^xx-1;plot(xx,yy,'r');
% examine
yy=ROOTexp.^DirectDistance-1

% example: log9(x) = log(x) / log(9)
yy=log(xx+1)/log(ROOTexp);plot(xx,yy,'b');
% examine
yy=log(DirectDistance+1)/log(ROOTexp)


%%
conditions = [0 2*atan(1/4) 2*atan(1/2) 2*atan(1)];

for i=1:length(conditions)
    
    rad_ang=conditions(i); 

    reach  = @(radius,rad_ang)  [radius*cos(rad_ang);  radius*sin(rad_ang)]; % Reach Trace Function For radius in pixels
    N      = radius;
    r_radius = linspace(0,radius,N);
    xy_r   = reach(r_radius,rad_ang);
    x      = xy_r(1,:);
    y      = xy_r(2,:);
    plot(x,y,'g');
end


%%



% [x,y] = drawReach(radius, rad_ang, xCenter, yCenter);
% plot(x,y,'go');
% xlim([0 1920]);ylim([0 1080]);
% hold on; 
% plot(xCenter, yCenter,'ro');
% plot(x(end),y(end),'ro');
% set(gca, 'YDir', 'reverse');