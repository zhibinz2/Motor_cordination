%% load xdf from LSL
% cd C:\Users\NIRS\Documents\zhibin\2022032801
cd /home/zhibin/Documents/Artinis_NIRS/zhibin/20220410
fileList = dir('*.xdf');
fileList.name
[streams,fileheader] = load_xdf(fileList.name);
[streams,fileheader] = load_xdf('20220415.xdf');

info=streams{1,2}.info;
info.name % to comfirm it is the stream of OxySoft
time_series=streams{1,2}.time_series;
time_stamps=streams{1,2}.time_stamps; num2str(time_stamps)
segments=streams{1,2}.segments;
srate=info.effective_srate;
clock_offsets=info.clock_offsets;

% open LSL mapping and find the oxy doxy channels
open 20220415.xlsx

figure('units','normalized','outerposition',[0 0 1 0.6]);
plotx(time_stamps',time_series(17:2:59,:)','r');
hold on;
plotx(time_stamps',time_series(18:2:60,:)','b');
xlabel('time-stamps from LSL')
xlabel('time')

% examine the time of ambient light noise
signal3=sort(time_series(3,:));
find(time_series(3,:)==signal3(end))
% examine werid timestamps (don't use LSL time_stamps)
plot(time_stamps','r.') % time stamps values is high for a short time first then dropped

% event marker stream from LSL
info=streams{1,1}.info;
info.name % confirm that it is OxySoft Event Marker
time_series2=streams{1,1}.time_series; % noting here
time_stamps2=streams{1,1}.time_stamps; % noting here

% plot Markers on top the signal (not work on event markers from LSL in OxySoft Event Markers)
hold on;
for i=1:length(time_stamps2)
    xline(time_stamps2(i),'m');
end
hold off;

% event marker stream from LSL
info=streams{1,3}.info;
info.name % confirm that it is MyMarkerStream
time_series3=streams{1,3}.time_series; % C C A B A B A B ...
time_stamps3=streams{1,3}.time_stamps; % 
num2str(time_stamps3)

% plot Markers on top the signal (work on event markers from LSL in MyMarkerStream)
hold on;
for i=1:length(time_stamps3)
    xline(time_stamps3(i),'m');
end
hold off;

%% replace the start and end artifacts with each channel's mean
% open LSL mapping and find the oxy doxy channels
open 20220415.xlsx

NIRS44chan=time_series(17:60,:);

find(time_stamps==time_stamps3(1)) % no return, need to find the closest estimate of time for the events
% find the cloest event indices in time_stamps
EventIndices=zeros(length(time_stamps3),1);
for i=1:length(time_stamps3) % i=StartStim2 % syncopation has one press less than condition 1 and 3
    [minValue,closetIndex]=min(abs(time_stamps-time_stamps3(i))); % closetIndex in BottonPressTime
    EventIndices(i)=closetIndex;
end

figure('units','normalized','outerposition',[0 0 1 0.3]);
plotx(NIRS44chan(1:2:43,:)','r');
hold on;
plotx(NIRS44chan(2:2:44,:)','b');
xlabel('time-samples from LSL')
hold off;

[x, ~] = ginput(2); % read two mouse clicks on the plot % x were index, y were real values
% get the proximate index
string(x)
ind1=round(x(1))
ind2=round(x(2))

% replace the start and end artifacts with each channel's mean
Mean44chan=mean(NIRS44chan(:,ind1:ind2),2);
NIRS44chan(:,1:ind1)=Mean44chan*ones(1,ind1);
NIRS44chan(:,ind2:end)=Mean44chan*ones(1,length(time_stamps)-ind2+1);
plotx(NIRS44chan(1:2:43,:)','r');
hold on;
plotx(NIRS44chan(2:2:44,:)','b');
hold off;

% detrend and band-pass filter (0.05-0.8 Hz)
% detrend
detrend44chan=detrend(NIRS44chan',1);
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(detrend44chan(:,1:2:43),'r');hold on;
plot(detrend44chan(:,2:2:44),'b');
title('detrend oxyval & dxyval'); 

% high pass (no paddings needed? for now)
Hd = makefilter(Fs,0.05,0.04,6,20,0); 
 % tweak the last line (doplot) of makefilter.m with 'freqz(Hd,0.01:0.01:60,sr)' in order to show finer details
% Hd = makefilter(Fs,0.05,0.04,6,20,1); 
which filtfilt.m -all
% cd 'C:\Program Files\MATLAB\R2021b\toolbox\signal\signal\'
cd /usr/local/MATLAB/R2022a/toolbox/signal/signal/
which filtfilt.m
flt44chan1=filtfilthd(Hd,double(detrend44chan));
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(flt44chan1(:,1:2:43),'r');hold on;
plot(flt44chan1(:,2:2:44),'b');

% low pass (this will create short edge artifact, no padding for now)
Hd = makefilter(Fs,0.8,0.9,6,20,0);
% Hd = makefilter(Fs,0.8,0.9,6,20,1); 
which filtfilt.m -all
cd 'C:\Program Files\MATLAB\R2021b\toolbox\signal\signal\'
which filtfilt.m
flt44chan2=filtfilthd(Hd,flt44chan1);
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(flt44chan2(:,1:2:43),'r');hold on;
plot(flt44chan2(:,2:2:44),'b');

%% organize into time x trials
% Pick the event onset 
diff(EventIndices)/srate
numTrials=36;
numChans=44;
trial_length=round(srate*40); % 10s baseline + 10s stimulus + 20 s rest
TrialTime=([1:trial_length]*(1/srate)-10*srate*(1/srate));
data_trials=zeros(trial_length,numChans,numTrials);
for i=1:numTrials % i=numTrials
    data_trials(:,:,i)=flt44chan2((EventIndices(2+2*i-1)-round(srate*10)):(EventIndices(2+2*i-1)+round(srate*30)-2),:);
end

% examine one trial
figure;plotx(TrialTime,data_trials(:,:,20));
% examine from trial to trial
figure;
for i=1:numTrials
    plot(TrialTime,data_trials(:,1:2:43,i),'r',TrialTime,data_trials(:,2:2:44,i),'b');
    pause(0.5);
end

% baseline normalization (no need. it will make doxy values big)
normalized_trials=zeros(trial_length,numChans,numTrials);
for i=1:numTrials
    baselineMean=mean(abs(data_trials(1:round(Fs*10),:,i)),1);
    normalized_trials(:,:,i)=data_trials(:,:,i)./(baselineMean.*ones(trial_length,1));
end
% examine one trial
figure;plot(TrialTime,normalized_trials(:,:,10));
% examine trial to trial
figure;
for i=1:numTrials
    subplot(2,1,1); plot(TrialTime,data_trials(:,1:2:43,i),'r',TrialTime,data_trials(:,2:2:44,i),'b');
    subplot(2,1,2); plot(TrialTime,normalized_trials(:,1:2:43,i),'r',TrialTime,normalized_trials(:,2:2:44,i),'b');pause(2);
end

%% plot left and right vs condition 
% load stimulus data
cd C:\Users\NIRS\Documents\zhibin\2022032801
cd /home/zhibin/Documents/Artinis_NIRS/zhibin/20220415
load('2022032801.mat'); allPerm;
allPerm=allPerm(1:numTrials);
trialsL=find(allPerm==1);
trialsR=find(allPerm==2);

% separate left and right stimulus condition
data_trialsL=data_trials(:,:,trialsL);
data_trialsR=data_trials(:,:,trialsR);

% channels
oxychanL=1:2:21;
dxychanL=2:2:22;
oxychanR=23:2:43;
dxychanR=24:2:44;

% Only the occipital lobe
oxychanL=[1 5 7 9 11 17];
dxychanL=[2 6 8 10 12 18];
oxychanR=[23 27 29 31 33 35];
dxychanR=[24 28 30 32 34 36];

% examine trial to trial
figure;
for i=1:numTrials
    subplot(2,1,1); 
    plot(TrialTime,data_trials(:,oxychanL,i),'r',TrialTime,data_trials(:,dxychanL,i),'b');
    title(['Left Hemisphere condition ' num2str(allPerm(i))]);xline(0,'m');ylim([-2 2]);
    subplot(2,1,2); 
    plot(TrialTime,data_trials(:,oxychanR,i),'r',TrialTime,data_trials(:,dxychanR,i),'b');
    title(['Right Hemisphere condition ' num2str(allPerm(i))]);xline(0,'m');ylim([-2 2]);
    pause(2);
end




