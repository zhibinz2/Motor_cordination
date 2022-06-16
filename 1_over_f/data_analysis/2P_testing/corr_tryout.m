x = 0.1*rand(3,1); 
y = 0.1*rand(3,1);
z1 = corr(0.6+x,0.6+y)
z2 = corr(0.6+x,0.7+y)


subplot(2,1,1);
plot(0.6+x,'r');
hold on;
plot(0.6+y,'b');
hold off;
legend('0.6+x','0.6+y');
title(['z1 ' num2str(z1)]);
subplot(2,1,2);
plot(0.6+x,'r');
hold on;
plot(0.7+y,'b');
hold off;
legend('0.6+x','0.7+y');
title(['z1 ' num2str(z1)]);

%%
figure;
plot(rand(10,1),'r');hold on; plot(0.5+rand(10,1),'b');hold off;
%%
t = [0.01:0.01:1];
% two sin waves
f1 = 2
x1 = sin(2*pi*f1*t);

f2 = 2
x2 = sin(2*pi*f2*t);
phase_shift=pi;
offset=3;
x2 = offset + sin(2*pi*f2*t+phase_shift);

figure
clf
subplot(3,1,1);

hold off
plot(t,x1,'r');
hold on
plot(t,x2,'b');
legend('x1','x2');
title(['phase-shift = ' num2str(phase_shift) '    Corr = ' num2str(corr(x1',x2')) ...
    '   offset = ' num2str(offset)]);
hold off

subplot(3,1,2);

subplot(3,1,3);
% [c,lags] = xcorr(x1,x2,100,'coeff');
% stem(lags,c);xline(0,'r');