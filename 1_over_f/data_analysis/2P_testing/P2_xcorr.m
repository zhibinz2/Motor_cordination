%% Raw data
timeL;samplesL;TRIGGERindL;srL;channels_infoL;
timeR;samplesR;TRIGGERindR;srR;channels_infoR;
% code for the analysis
cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing
% data
clear;close all;
cd /ssd/zhibin/1overf/20220515_2P
cd /ssd/zhibin/1overf/20220517_2P
cd /ssd/zhibin/1overf/20220518_2P
%% Segmented data
% Segment EEG
% mixedsigL=mixedsigL';
% mixedsigR=mixedsigR';
% PacersL=PacerTimeIndL([1 2 3 4 5 34 35 36 37 66 67 68 69 98 99 100 101 130 131 132 133 162 163 164 165 194]);
% PacersR=PacerTimeIndR([1 2 3 4 5 34 35 36 37 66 67 68 69 98 99 100 101 130 131 132 133 162 163 164 165 194]);
EEGOpenEyeRestingL=mixedsigL(PacersL(1):PacersL(2),:);
EEGCloseEyeRestingL=mixedsigL(PacersL(3):PacersL(4),:);
EEGCondi1L=mixedsigL(PacersL(6):PacersL(7),:);
EEGResting2L=mixedsigL(PacersL(7):PacersL(8),:);
EEGCondi2L=mixedsigL(PacersL(10):PacersL(11),:);
EEGResting3L=mixedsigL(PacersL(11):PacersL(12),:);
EEGCondi3L=mixedsigL(PacersL(14):PacersL(15),:);
EEGResting4L=mixedsigL(PacersL(15):PacersL(16),:);
EEGCondi4L=mixedsigL(PacersL(18):PacersL(19),:);
EEGResting5L=mixedsigL(PacersL(19):PacersL(20),:);
EEGCondi5L=mixedsigL(PacersL(22):PacersL(23),:);
EEGResting6L=mixedsigL(PacersL(23):PacersL(24),:);
EEGCondi6L=mixedsigL(PacersL(26):end,:);

EEGOpenEyeRestingR=mixedsigR(PacersR(1):PacersR(2),:);
EEGCloseEyeRestingR=mixedsigR(PacersR(3):PacersR(4),:);
EEGCondi1R=mixedsigR(PacersR(6):PacersR(7),:);
EEGResting2R=mixedsigR(PacersR(7):PacersR(8),:);
EEGCondi2R=mixedsigR(PacersR(10):PacersR(11),:);
EEGResting3R=mixedsigR(PacersR(11):PacersR(12),:);
EEGCondi3R=mixedsigR(PacersR(14):PacersR(15),:);
EEGResting4R=mixedsigR(PacersR(15):PacersR(16),:);
EEGCondi4R=mixedsigR(PacersR(18):PacersR(19),:);
EEGResting5R=mixedsigR(PacersR(19):PacersR(20),:);
EEGCondi5R=mixedsigR(PacersR(22):PacersR(23),:);
EEGResting6R=mixedsigR(PacersR(23):PacersR(24),:);
EEGCondi6R=mixedsigR(PacersR(26):end,:);


%% segment button presses
BottonPresTimeL01; BottonPresTimeR01;

BPOpenEyeRestingL=BottonPresTimeL01(PacersL(1):PacersL(2),:);
BPCloseEyeRestingL=BottonPresTimeL01(PacersL(3):PacersL(4),:);
BPCondi1L=BottonPresTimeL01(PacersL(6):PacersL(7),:); plot(BPCondi1L,'r');
BPResting2L=BottonPresTimeL01(PacersL(7):PacersL(8),:);
BPCondi2L=BottonPresTimeL01(PacersL(10):PacersL(11),:);
BPResting3L=BottonPresTimeL01(PacersL(11):PacersL(12),:);
BPCondi3L=BottonPresTimeL01(PacersL(14):PacersL(15),:);
BPResting4L=BottonPresTimeL01(PacersL(15):PacersL(16),:);
BPCondi4L=BottonPresTimeL01(PacersL(18):PacersL(19),:);
BPResting5L=BottonPresTimeL01(PacersL(19):PacersL(20),:);
BPCondi5L=BottonPresTimeL01(PacersL(22):PacersL(23),:);
BPResting6L=BottonPresTimeL01(PacersL(23):PacersL(24),:);
BPCondi6L=BottonPresTimeL01(PacersL(26):end,:);

BPOpenEyeRestingR=BottonPresTimeR01(PacersR(1):PacersR(2),:);
BPCloseEyeRestingR=BottonPresTimeR01(PacersR(3):PacersR(4),:);
BPCondi1R=BottonPresTimeR01(PacersR(6):PacersR(7),:);plot(BPCondi1R,'b');
BPResting2R=BottonPresTimeR01(PacersR(7):PacersR(8),:);
BPCondi2R=BottonPresTimeR01(PacersR(10):PacersR(11),:);
BPResting3R=BottonPresTimeR01(PacersR(11):PacersR(12),:);
BPCondi3R=BottonPresTimeR01(PacersR(14):PacersR(15),:);
BPResting4R=BottonPresTimeR01(PacersR(15):PacersR(16),:);
BPCondi4R=BottonPresTimeR01(PacersR(18):PacersR(19),:);
BPResting5R=BottonPresTimeR01(PacersR(19):PacersR(20),:);
BPCondi5R=BottonPresTimeR01(PacersR(22):PacersR(23),:);
BPResting6R=BottonPresTimeR01(PacersR(23):PacersR(24),:);
BPCondi6R=BottonPresTimeR01(PacersR(26):end,:);

%

%% segment feedbacks from the other 
FeedbTimeL01; FeedbTimeR01;

FBOpenEyeRestingL=FeedbTimeL01(PacersL(1):PacersL(2),:);
FBCloseEyeRestingL=FeedbTimeL01(PacersL(3):PacersL(4),:);
FBCondi1L=FeedbTimeL01(PacersL(6):PacersL(7),:); plot(FBCondi1L,'b');
FBResting2L=FeedbTimeL01(PacersL(7):PacersL(8),:);
FBCondi2L=FeedbTimeL01(PacersL(10):PacersL(11),:);
FBResting3L=FeedbTimeL01(PacersL(11):PacersL(12),:);
FBCondi3L=FeedbTimeL01(PacersL(14):PacersL(15),:);
FBResting4L=FeedbTimeL01(PacersL(15):PacersL(16),:);
FBCondi4L=FeedbTimeL01(PacersL(18):PacersL(19),:);
FBResting5L=FeedbTimeL01(PacersL(19):PacersL(20),:);
FBCondi5L=FeedbTimeL01(PacersL(22):PacersL(23),:);
FBResting6L=FeedbTimeL01(PacersL(23):PacersL(24),:);
FBCondi6L=FeedbTimeL01(PacersL(26):end,:);

FBOpenEyeRestingR=FeedbTimeR01(PacersR(1):PacersR(2),:);
FBCloseEyeRestingR=FeedbTimeR01(PacersR(3):PacersR(4),:);
FBCondi1R=FeedbTimeR01(PacersR(6):PacersR(7),:); plot(FBCondi1R,'r');
FBResting2R=FeedbTimeR01(PacersR(7):PacersR(8),:);
FBCondi2R=FeedbTimeR01(PacersR(10):PacersR(11),:);
FBResting3R=FeedbTimeR01(PacersR(11):PacersR(12),:);
FBCondi3R=FeedbTimeR01(PacersR(14):PacersR(15),:);
FBResting4R=FeedbTimeR01(PacersR(15):PacersR(16),:);
FBCondi4R=FeedbTimeR01(PacersR(18):PacersR(19),:);
FBResting5R=FeedbTimeR01(PacersR(19):PacersR(20),:);
FBCondi5R=FeedbTimeR01(PacersR(22):PacersR(23),:);
FBResting6R=FeedbTimeR01(PacersR(23):PacersR(24),:);
FBCondi6R=FeedbTimeR01(PacersR(26):end,:);

%% %% segment EMG
filtered_EMGL; filtered_EMGR;

EMGOpenEyeRestingL=filtered_EMGL(PacersL(1):PacersL(2),:);
EMGCloseEyeRestingL=filtered_EMGL(PacersL(3):PacersL(4),:);
EMGCondi1L=filtered_EMGL(PacersL(6):PacersL(7),:);
EMGResting2L=filtered_EMGL(PacersL(7):PacersL(8),:);
EMGCondi2L=filtered_EMGL(PacersL(10):PacersL(11),:);
EMGResting3L=filtered_EMGL(PacersL(11):PacersL(12),:);
EMGCondi3L=filtered_EMGL(PacersL(14):PacersL(15),:);
EMGResting4L=filtered_EMGL(PacersL(15):PacersL(16),:);
EMGCondi4L=filtered_EMGL(PacersL(18):PacersL(19),:);
EMGResting5L=filtered_EMGL(PacersL(19):PacersL(20),:);
EMGCondi5L=filtered_EMGL(PacersL(22):PacersL(23),:);
EMGResting6L=filtered_EMGL(PacersL(23):PacersL(24),:);
EMGCondi6L=filtered_EMGL(PacersL(26):end,:);

EMGOpenEyeRestingR=filtered_EMGR(PacersR(1):PacersR(2),:);
EMGCloseEyeRestingR=filtered_EMGR(PacersR(3):PacersR(4),:);
EMGCondi1R=filtered_EMGR(PacersR(6):PacersR(7),:);
EMGResting2R=filtered_EMGR(PacersR(7):PacersR(8),:);
EMGCondi2R=filtered_EMGR(PacersR(10):PacersR(11),:);
EMGResting3R=filtered_EMGR(PacersR(11):PacersR(12),:);
EMGCondi3R=filtered_EMGR(PacersR(14):PacersR(15),:);
EMGResting4R=filtered_EMGR(PacersR(15):PacersR(16),:);
EMGCondi4R=filtered_EMGR(PacersR(18):PacersR(19),:);
EMGResting5R=filtered_EMGR(PacersR(19):PacersR(20),:);
EMGCondi5R=filtered_EMGR(PacersR(22):PacersR(23),:);
EMGResting6R=filtered_EMGR(PacersR(23):PacersR(24),:);
EMGCondi6R=filtered_EMGR(PacersR(26):end,:);

%% load data
clear /ssd/zhibin/1overf/20220515_2P/Segmented_data


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
% sliding window of errors

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
figure;plot(lags./2,r);xlabel('time [ms]');ylabel('cross correlation');title('xcorr on boolean');

% xcorr on intervals
time_series1=Calinterval(FBCondi2L');% plot(time_series1);
time_series2=Calinterval(BPCondi2L');% plot(time_series2);
TimeLength=min([length(time_series1) length(time_series2)]);
[r,lags]=xcorr(time_series1(1:TimeLength), time_series2(1:TimeLength), 1000,'normalized');
figure;plot(lags./2,r);xlabel('time [ms]');ylabel('cross correlation');title('xcorr on intervals');



