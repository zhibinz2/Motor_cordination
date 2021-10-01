x=linspace(0,2*pi,100);
y=sin(x);
subplot(1,3,1);
plot(x,y,'r');
hold on;
% y2=sin(x+0.25*pi);% y2 is leading by 25 samples or value of 3.14/4=0.78 seconds
y2=circshift(y,-10);% blue is leading by 10 samples
plot(x,y2,'b');


[r,lags]=xcorr(y,y2);
subplot(1,3,2);
plot(lags,r);

subplot(1,3,3)
[pks,locs] = findpeaks(r,lags,'MinPeakProminence',40,'Annotate','extents');
findpeaks(r,lags,'MinPeakProminence',40,'Annotate','extents');
hold off;

%%
n = 0:15;
x = 0.84.^n;
plot(x,'r');hold on;
y = circshift(x,5);
plot(y,'b');hold off;
[c,lags] = xcorr(x,y);
stem(lags,c)



  N = 16;
  n = 0:N-1;
  a = 0.92;
  b = 0.92;
  xa = a.^n;
  xb = b.^(n+1);
  [r,lags] = xcorr(xa,xb);
  stem(lags,r)