x1=zeros(1,101);
x1([51 53 56])=1;
x2=zeros(1,101);
x2([46 49 51])=1;
xcorr(x1,x2,101);
[r,lags]=xcorr(x1,x2,101,'normalized');
plot(lags,r);xlabel('lags');ylabel('Xcorr');

%%
n = 0:15;
x = 0.84.^n; plot(x); plot(x);
y = circshift(x,5)-0.5; plot(y);
[c,lags] = xcorr(x,y,15,'normalized');
[c,lags] = xcorr(x,y,15,'coeff');
stem(lags,c)

%% EMG
[r,lags]=xcorr(EMG_L, BP_L, 1000,'normalized');
plot(lags,r);xlabel('time [ms]');ylabel('Xcorr');
plot(EMG_L);
EMG_L0=EMG_L-100;plot(EMG_L0);
[r,lags]=xcorr(EMG_L0, BP_L, 1000);
plot(lags,r);xlabel('lags');ylabel('Xcorr');
[r,lags]=xcorr(EMG_L0, BP_L,'normalized');
plot(lags,r);xlabel('lags');ylabel('Xcorr');
[r,lags]=xcorr(EMG_L0, BP_L,500,'normalized');
plot(lags,r);xlabel('lags');ylabel('Xcorr');

%% EEG
for j = 1:32
[y(:,j),lags] = xcorr(BPCondi2L(1:end),EEGCondi2L(1:end,j),2000,'coeff');
end
plot(lags,y);

figure;
for j = 1:32
[y(:,j),lags] = xcorr(EEGCondi2L(1:end,j),BPCondi2L(1:end),2000,'coeff');
end
plotx(lags,y);

%% hist of BP interval
cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing
win=20;
y=smoothing(Calinterval(BPCondi2L'),win);hist(y)
y=Calinterval(BPCondi2L');hist(y./2);xlabel('BP-intervals (ms)');ylabel('number');
