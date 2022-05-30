%% Create two sinusoidal signal
open /home/zhibin/Downloads/HNL/FFT 2019/Demo code/20191212/20191212code.m
open /home/zhibin/Downloads/HNL/FFT 2019/Demo code/20191212/code20191212.m
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
subplot(4,1,1);
plot(t,x1,'r');
hold on
plot(t,x2,'b');title(['phase-shift = ' num2str(phase_shift)]);
legend('x1','x2');
% plot(t,x1.*x2,'k')
sum(x1.*x2)
% sum to 0, when f 2 and 3

%% Xspetra
cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing

%% Method1: multply by conjugate (work!)
maxfreq = 5;
sr=100;
T=1;
df=1/T;% sr/(sr*T);
nbins = ceil(maxfreq/df) + 1;
freqs = [0:(nbins-1)]*df;
fcoef1 = fft(ndetrend(x1',1),[],1)/length(x1);
fcoef2 = fft(ndetrend(x2',1),[],1)/length(x2);

cprod=[];
for i=1:nbins
    cprod(i)=fcoef1(i)*conj(fcoef2(i))/length(x1);
end

% look at all frequency
subplot(4,1,2);plot(freqs,cprod);title('cprod');xline(2,'m');
subplot(4,1,3);plot(freqs,abs(cprod)); title('abs(cprod)');xline(2,'m');% cross spectra
subplot(4,1,4);plot(freqs,angle(cprod)); title('angle(cprod)');% phase difference
xline(2,'m');yline(-phase_shift,'b');

%% two channels into data_trials (not working)
time_series1=x1'; 
time_series2=x2'; 
data=cat(3,time_series1,time_series2); % treat the two time series as two trials in the data

maxfreq = 5;
goodepochs = [1:size(data,3)];
win = [1:size(data,1)];

rate=100;
df = rate/length(win);
nbins = ceil(maxfreq/df) + 1;
freqs = [0:(nbins-1)]*df;
fcoef = fft(ndetrend(data(win,:,:),1),[],1)/length(win);
eppow = abs(squeeze(mean(fcoef(1:nbins,:,goodepochs),3))).^2; % averaged across all good trials (power spectra)
figure;plot(freqs,eppow);% plot just one channel
pow = squeeze(var(fcoef(1:nbins,:,goodepochs),[],3)); % variance across all good trials (power spectra)
figure;plot(freqs,pow);% plot just one channel

cprod=[];corr=[];
for k = 1:nbins % for each frequency
   sf = corrcoef(transpose(squeeze(fcoef(k,:,goodepochs)))); % coherence (correlation coefficients)
   corr(k,:,:) = sf; % coherence  (correlation coefficients)
%  coh(k,:,:) = abs(sf).^2; % coherence  (correlation coefficients)
%  phase(k,:,:) = angle(sf);
  cprod(k,:,:) = cov(transpose(squeeze(fcoef(k,:,goodepochs))));% covarience (cross spectra)
end

% look at one frequency
f= 2; % Hz 
ind=round(f/df)+1;
figure;imagesc(abs(squeeze(cprod(ind,:,:)))); xlabel('channels');ylabel('channels'); colorbar;% cross spetra
figure;imagesc(angle(squeeze(cprod(ind,:,:)))); xlabel('channels');ylabel('channels'); colorbar;% phase differences

% look at all frequency
figure;plot(freqs,cprod);
figure;plot(freqs,abs(cprod)); % cross spectra
figure;plot(freqs,angle(cprod)); % phase difference
figure;plot(freqs,corr);