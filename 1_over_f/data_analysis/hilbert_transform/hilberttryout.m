%% https://www.mathworks.com/help/signal/ug/hilbert-transform.html
t = 0:1/1024:1;
x = sin(2*pi*60*t);
figure;
plot(t(1:50),x(1:50),'r');hold on;

y = hilbert(x);

plot(t(1:50),real(y(1:50)))
hold on
plot(t(1:50),imag(y(1:50)))
hold off
axis([0 0.05 -1.1 2])
legend('Real Part','Imaginary Part')
%% https://www.mathworks.com/help/signal/ug/envelope-extraction-using-the-analytic-signal.html
t = 0:1e-4:0.1;
x = (1+cos(2*pi*50*t)).*cos(2*pi*1000*t);

plot(t,x)
xlim([0 0.04]);

clf;
y = hilbert(x);
env = abs(y);
plot_param = {'Color', [0.6 0.1 0.2],'Linewidth',2}; 

plot(t,x)
hold on
plot(t,[-1]*env,plot_param{:})
plot(t,[1]*env,plot_param{:})
plot(t,[-1;1]*env,plot_param{:})
hold off
xlim([0 0.04])
title('Hilbert Envelope')


