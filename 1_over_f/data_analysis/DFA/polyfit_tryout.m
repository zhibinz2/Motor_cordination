x = 1:50;
y = -0.3*x + 2*randn(1,50);
p = polyfit(x,y,1)

f = polyval(p,x);
plot(x,y,'o',x,f,'-')
legend('data','linear fit')
