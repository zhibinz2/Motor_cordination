angle=linspace(0,pi,100);
ytan=tan(angle);
plot(angle,ytan);ylabel('tangent value');xlabel('angle');

angles=linspace(0,pi/2,100);
ytan=tan(angles);
plot(angle,ytan);ylabel('tangent value');xlabel('angle');ylim([0 5]);

atan2(1)

atan(1)


%%
XL=-100; XR=100;
XL=-2000; XR=1;
XL=-1; XR=2000;
XL=-500; XR=1500;

d=sqrt(XL.^2+XR.^2)
theta=2*atan(-XL/XR)
XJ=d*cos(theta);
YJ=d*sin(theta);
plot(XJ,YJ,'go',XL,0,'ro',XR,0,'bo');xlim([-2000 2000]);ylim([-2000 2000]);
