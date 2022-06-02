x1=zeros(1,101);
x1([51 53 56])=1;
x2=zeros(1,101);
x2([46 49 51])=1;
xcorr(x1,x2,101);
[r,lags]=xcorr(x1,x2,101,'normalized');
plot(lags,r);xlabel('lags');ylabel('Xcorr');