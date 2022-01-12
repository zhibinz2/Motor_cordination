N=radius;
xx=linspace(0,radius,N) ;
yy=2.^(xx);plot(xx,yy,'r');
xline(0,'k--');yline(0,'k--');
xline(1,'k--');yline(1,'k--');
xlim([-1 3]);ylim([-1 3]);
hold on;plot(xx,xx,'k--');

yy=log10(xx);plot(xx,yy,'b');
