x1=zeros(1,101);
x1([51 53 56])=1;
x2=zeros(1,101);
x2([46 49 51])=1;
xcorr(x1,x2,101);
[r,lags]=xcorr(x1,x2,101,'normalized');
plot(lags,r);xlabel('lags');ylabel('Xcorr');

%%
n = 0:15;
x = 0.84.^n; plot(x); figure;plot(x);title('x');
y = circshift(x,5)-0.5; figure;plot(y);title('y');
[c,lags] = xcorr(x,y,15,'normalized');
[c,lags] = xcorr(x,y,15,'coeff');
figure;stem(lags,c)

%%
t = [0.01:0.01:1];
% two sin waves
f1 = 2
x1 = sin(2*pi*f1*t);

f2 = 2
x2 = sin(2*pi*f2*t);
phase_shift=pi*0.1
x2 = sin(2*pi*f2*t+phase_shift);

figure
clf
subplot(2,1,1);
plot(t,x1,'r');
hold on
plot(t,x2,'b');title(['phase-shift = ' num2str(phase_shift)]);
legend('x1','x2');

subplot(2,1,2);
[c,lags] = xcorr(x1,x2,100,'coeff');
stem(lags,c);xline(0,'r');

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
