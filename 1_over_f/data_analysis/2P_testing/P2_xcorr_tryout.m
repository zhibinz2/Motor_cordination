%% Raw data
timeL;samplesL;TRIGGERindL;srL;channels_infoL;
timeR;samplesR;TRIGGERindR;srR;channels_infoR;
% code for the analysis
cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing
% load data
clear;close all;
cd /ssd/zhibin/1overf/20220515_2P
cd /ssd/zhibin/1overf/20220517_2P
cd /ssd/zhibin/1overf/20220518_2P
%% Segment data process
open /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing/P2_EEG_processing.m

%% load segmented data
clear
cd /ssd/zhibin/1overf/20220515_2P/Segmented_data

%% Xspectra Method1: multply by conjugate (work!)
open /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing/Xspectra_tryout.m

% two channels
x1=EEGCondi2L(:,5); % time_series1=BPCondi2L;
x2=EEGCondi2R(:,5); 

maxfreq = 5;
sr=2000; % unique([srL srR]);
T=size(EEGCondi2R,1)/sr; % in second
df=1/T;
nbins = ceil(maxfreq/df) + 1;
freqs = [0:(nbins-1)]*df;

fcoef1 = fft(ndetrend(x1,1),[],1)/length(x1);
fcoef2 = fft(ndetrend(x2,1),[],1)/length(x2);

cprod=[];
for i=1:nbins
    cprod(i)=fcoef1(i)*conj(fcoef2(i))/length(x1);
end

% look at all frequency
figure;
subplot(4,1,1);plot(x1,'r');hold on;plot(x2,'b');title('two time series');
subplot(4,1,2);plot(freqs,cprod);title('cprod');
subplot(4,1,3);plot(freqs,abs(cprod)); title('abs(cprod)');% cross spectra
subplot(4,1,4);plot(freqs,angle(cprod)); title('angle(cprod)');% phase difference

% two channels
x1=EEGCondi2L; % time_series1=BPCondi2L;
x2=EEGCondi2R; 

maxfreq = 5;
sr=2000; % unique([srL srR]);
SamplesLength=size(x2,1);
T=SamplesLength/sr; % in second
df=1/T;
nbins = ceil(maxfreq/df) + 1;
freqs = [0:(nbins-1)]*df;

fcoef1 = fft(ndetrend(x1,1),[],1)/SamplesLength;
fcoef2 = fft(ndetrend(x2,1),[],1)/SamplesLength;

cprod=[];
for i=1:nbins
    cprod(i,:)=fcoef1(i,:).*conj(fcoef2(i,:))/length(x1);
end

% look at all frequency
figure;
subplot(4,1,1);plot(x1);hold on;plot(x2);title('two time series');
subplot(4,1,2);plot(freqs,cprod);title('cprod');
subplot(4,1,3);plot(freqs,abs(cprod)); title('abs(cprod)');% cross spectra
subplot(4,1,4);plot(freqs,angle(cprod)); title('angle(cprod)');% phase difference

%% cross spectra (allspectra)
open allspectra.m

% all channels
time_series1=EEGCondi2L; % time_series1=BPCondi2L;
time_series2=EEGCondi2R; % time_series2=FBCondi2R; % time_series2=BPCondi2L;
data=cat(3,time_series1,time_series2); % treat the two time series as two trials in the data

maxfreq = 50;
goodepochs = [1:size(data,3)];
win = [1:size(data,1)];

df = rate/length(win);
nbins = ceil(maxfreq/df) + 1;
freqs = [0:(nbins-1)]*df;
fcoef = fft(ndetrend(data(win,:,:),1),[],1)/length(win);
eppow = abs(squeeze(mean(fcoef(1:nbins,:,goodepochs),3))).^2; % averaged across all good trials
plot(freqs,eppow(:,5));% plot just one channel
pow = squeeze(var(fcoef(1:nbins,:,goodepochs),[],3)); % variance
plot(freqs,pow(:,5));% plot just one channel

cprod=[];corr=[];
for k = 1:nbins % for each frequency
   sf = corrcoef(transpose(squeeze(fcoef(k,:,goodepochs)))); % coherence
%  coh(k,:,:) = abs(sf).^2;
%  phase(k,:,:) = angle(sf);
  cprod(k,:,:) = cov(transpose(squeeze(fcoef(k,:,goodepochs))));% cross spectra
  corr(k,:,:) = sf; % coherence
end

f= 10; % Hz % looking at one frequency
ind=round(f/df)+1;
imagesc(abs(squeeze(cprod(ind,:,:)))); xlabel('channels');ylabel('channels'); colorbar;% cross spetra
imagesc(angle(squeeze(cprod(ind,:,:)))); xlabel('channels');ylabel('channels'); colorbar;% phase differences

% label all channels
xticks([1:32]);
xticklabels(AllchanNames); % all 32 channels
xtickangle(90);
yticks([1:32]);
yticklabels(AllchanNames) % all 32 channels


% two channels
time_series1=EEGCondi2L(:,5); % time_series1=BPCondi2L;
time_series2=EEGCondi2R(:,5); % time_series2=FBCondi2R; % time_series2=BPCondi2L;
data=cat(3,time_series1,time_series2); % treat the two time series as two trials in the data

maxfreq = 50;
goodepochs = [1:size(data,3)];
win = [1:size(data,1)];

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
   sf = corrcoef(transpose(squeeze(fcoef(k,:,goodepochs)))); % coherence
   corr(k,:,:) = sf; % coherence
%  coh(k,:,:) = abs(sf).^2; % coherence (
%  phase(k,:,:) = angle(sf);
  cprod(k,:,:) = cov(transpose(squeeze(fcoef(k,:,goodepochs))));% covarience, cross spectra
  
end

% look at one frequency
f= 10; % Hz 
ind=round(f/df)+1;
figure;imagesc(abs(squeeze(cprod(ind,:,:)))); xlabel('channels');ylabel('channels'); colorbar;% cross spetra
figure;imagesc(angle(squeeze(cprod(ind,:,:)))); xlabel('channels');ylabel('channels'); colorbar;% phase differences

% look at all frequency
figure;plot(freqs,cprod);
figure;plot(freqs,abs(cprod)); % cross spectra
figure;plot(freqs,angle(cprod)); % phase difference
figure;plot(freqs,corr);

%% autocorr
% https://www.mathworks.com/help/econ/autocorr.html

% autocorr on boolean (not meaningful)
time_series=FBCondi2L;
[acf,lags] = autocorr(time_series);
figure; autocorr(time_series,500) % number of lags

% autocorr on intervals
time_series=Calinterval(FBCondi2L');
figure; autocorr(time_series,500) % number of lags
title('autocorr on intervals');

% autocorr on errors
time_series1=find(FBCondi2R==1);% plot(time_series1);
time_series2=find(BPCondi2R==1);% plot(time_series2);
% TimeLength=min([length(time_series1) length(time_series2)]);
% TimeLength=max([length(time_series1) length(time_series2)]);
TimeLength=length(time_series1);
Error=[];
% find the time difference with the closest botton press 
for i=1:TimeLength % i=StartStim2 % syncopation has one press less than condition 1 and 3
    [minValue,closetIndex]=min(abs(time_series2-time_series1(i))); % closetIndex in BottonPressTime
    Error(i)=time_series2(closetIndex)-time_series1(i);
end

% sliding window of 10 errors
function [ErrorSmooth] = smoothError(Error,win);
ErrorSmooth=[];
for i=1:length(Error)-20
    ErrorSmooth(i)=mean(Error(i:i+9));
end
end

% plot
figure;
subplot(3,1,1);
plot(Error/2000,'r.');title('error');xlabel('taps');ylabel('timing error (s)');
subplot(3,1,2);
autocorr(Error,250);% number of lags
title('autocorr of error');
subplot(3,1,3);
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/DFA
[D,Alpha1]=DFA_main(Error);
title(['DFA of error']);
hold off

%% xcorr
cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing

% xcorr on boolean
time_series1=FBCondi2R';% plot(time_series1);
time_series2=BPCondi2R';% plot(time_series2);
TimeLength=min([length(time_series1) length(time_series2)]);
[r,lags]=xcorr(time_series1(1:TimeLength), time_series2(1:TimeLength), 1000, 'normalized');
figure;plot(lags./2,r);% devided by 2 so that the unit is ms
xlabel('time [ms]');ylabel('cross correlation');title('xcorr on boolean');

% xcorr on intervals
time_series1=Calinterval(FBCondi2L');% plot(time_series1);
time_series2=Calinterval(BPCondi2L');% plot(time_series2);
TimeLength=min([length(time_series1) length(time_series2)]);
[r,lags]=xcorr(time_series1(1:TimeLength), time_series2(1:TimeLength), 1000,'normalized');
figure;plot(lags./2,r);xlabel('time [ms]');ylabel('cross correlation');title('xcorr on intervals');



