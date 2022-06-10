%% tryout with calculation
A=[0 0 0 0 1 1 0 0 1 1 1 1 0 0 0 0 1 1 1 1 1 0 0 0 ];
pressInd=find([0 diff(A)]==1);
pressIntervals=[pressInd(1) diff(pressInd)];

%% load TMSi data
clear;close all;
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing;
cd /ssd/zhibin/1overf/20220515_2P
cd /ssd/zhibin/1overf/20220517_2P
cd /ssd/zhibin/1overf/20220518_2P
cd /ssd/zhibin/1overf/20220609_2P

[EEGfileNameL]=uigetfile('*.Poly5');% select the Left player EEG
Path_filenameL=[pwd '/' EEGfileNameL];
[timeL,samplesL,TRIGGERindL,srL,channels_infoL] = LoadTMSi(Path_filenameL);

[EEGfileNameR]=uigetfile('*.Poly5');% select the Right player EEG
Path_filenameR=[pwd '/' EEGfileNameR];
[timeR,samplesR,TRIGGERindR,srR,channels_infoR] = LoadTMSi(Path_filenameR);

%% Extract Trigger channel info
TriggersL=samplesL(TRIGGERindL,:)';unique(TriggersL)
TriggersR=samplesR(TRIGGERindR,:)';unique(TriggersR)

%% extract time points for key presses (maybe create a function to do the same for pacers,feedbacks)
% Left player
PresIndL=unique([find(TriggersL == 239);find(TriggersL == 235);find(TriggersL == 111 )]); % extract Index of real key presses in the values
plot(PresIndL,ones(1,length(PresIndL)),'ro'); % look at the above Index (one press produced several indices)
% ifi=1/100;waitframes=2;sr=2000;
% determine a threshold of numbers of frames in the button press interval
threshold = 1*1.5; % interval larger than 1 consecutive samples
% button presses recorded via DP port into stimtracker is sampled at 2000 Hz
% feedback via photocell to stimtracker has minial interval every other waitframes (minmal of 1000 Hz)
BottonPresTimeIndL=PresIndL(find([0 diff(PresIndL')]>threshold)); % exact index of key press onset in datatimes (reduce several indices into one)
% create a time series that assign botton presses as 1, all other as 0
BottonPresTimeL01=zeros(size(timeL));
BottonPresTimeL01(BottonPresTimeIndL)=1;
plot(BottonPresTimeL01,'r.')
sum(BottonPresTimeL01) % num of button presses


% Right player
PresIndR=unique([find(TriggersR == 239);find(TriggersR == 235);find(TriggersR == 111 )]); % extract Index of real key presses in the values
plot(PresIndR,ones(1,length(PresIndR)),'bo'); % look at the above Index (one press produced several indices)
% ifi=1/100;waitframes=2;sr=2000;
% determine a threshold of numbers of frames in the button press interval
threshold = 1*1.5; % interval larger than 1 consecutive samples
% button presses recorded via DP port into stimtracker is sampled at 2000 Hz
% feedback via photocell to stimtracker has minial interval every other waitframes (minmal of 1000 Hz)
BottonPresTimeInd=PresIndR(find([0 diff(PresIndR')]>threshold)); % exact index of key press onset in datatimes (reduce several indices into one)
% create a time series that assign botton presses as 1, all other as 0
BottonPresTimeR01=zeros(size(timeR));
BottonPresTimeR01(BottonPresTimeInd)=1;
plot(BottonPresTimeR01,'b.')
numPresR=sum(BottonPresTimeR01) % num of button presses


%% extract time points for feedbacks from other
% Left Player Recording
FeedbIndL=unique([find(TriggersL == 251);find(TriggersL == 235)]); % extract Index of real key presses in the values
plot(FeedbIndL,ones(1,length(FeedbIndL)),'bo'); % look at the above Index (one press produced several indices)
% ifi=1/100;waitframes=2;sr=2000;
% determine a threshold of numbers of frames in the button press interval
threshold = ifi*waitframes*srL*1.5; % number of samples within a interval larger than the waitframes
% button presses recorded via DP port into stimtracker is sampled at 2000 Hz
% feedback via photocell to stimtracker has minial interval every other waitframes (minmal of 1000 Hz)
FeedbTimeIndL=FeedbIndL(find([0 diff(FeedbIndL')]>threshold)); % exact index of key press onset in datatimes (reduce several indices into one)
% create a time series that assign botton presses as 1, all other as 0
FeedbTimeL01=zeros(size(timeL));
FeedbTimeL01(FeedbTimeIndL)=1;
plot(FeedbTimeL01,'b.')
numFbL=sum(FeedbTimeL01) % num of feedbacks from the other plus 4+2x5 resting photocells

% Right Player Recording
FeedbIndR=unique([find(TriggersR == 251);find(TriggersR == 235)]); % extract Index of real key presses in the values
plot(FeedbIndR,ones(1,length(FeedbIndR)),'ro'); % look at the above Index (one press produced several indices)
% ifi=1/100;waitframes=2;sr=2000;
% determine a threshold of numbers of frames in the button press interval
threshold = ifi*waitframes*srR*1.5; % number of samples within a interval larger than the waitframes
% button presses recorded via DP port into stimtracker is sampled at 2000 Hz
% feedback via photocell to stimtracker has minial interval every other waitframes (minmal of 1000 Hz)
FeedbTimeIndR=FeedbIndR(find([0 diff(FeedbIndR')]>threshold)); % exact index of key press onset in datatimes (reduce several indices into one)
% create a time series that assign botton presses as 1, all other as 0
FeedbTimeR01=zeros(size(timeR));
FeedbTimeR01(FeedbTimeIndR)=1;
plot(FeedbTimeR01,'r.')
numFbR=sum(FeedbTimeR01) % num of feedbacks from the other plus 4+2x5 resting photocells

%% extract time points for pacers (using this as index to segment data into conditions)
% Left player
% PacerIndL=unique([find(TriggersL == 127)]); % extract Index of real key presses in the values
% PacerIndL=unique([find(TriggersL == 127);find(TriggersL == 111)); % extract Index of real key presses in the values
PacerIndL=unique([find(TriggersL == 127);find(TriggersL == 111);find(TriggersL == 123)]); % extract Index of real key presses in the values
plot(PacerIndL,ones(1,length(PacerIndL)),'ro'); % look at the above Index (one press produced several indices)
% Pacer at 2Hz has intervals of 500ms (1000 samples)
% Pacer at 3Hz has interval of 333ms (666.6 samples)
% determine a threshold of numbers of samples in pacer interval
threshold = ifi*srL*0.75; % number of samples within a interval larger than the screen refresh rate
% button presses recorded via DP port into stimtracker is sampled at 2000 Hz
% feedback via photocell to stimtracker has minial interval every other waitframes (minmal of 1000 Hz)
PacerTimeIndL=PacerIndL(find([diff([0 PacerIndL'])]>threshold)); % exact index of key press onset in datatimes (reduce several indices into one)
% 2+2+30*6+2*5 = 194 % 2+30*6
% create a time series that assign botton presses as 1, all other as 0
PacerTimeL01=zeros(size(timeL));
PacerTimeL01(PacerTimeIndL)=1;
plot(PacerTimeL01,'r');
sum(PacerTimeL01) % check if == 194
% xlim([-1 length(PacerTimeL01)]);

% Right player
PacerIndR=unique([find(TriggersR == 127);find(TriggersR == 111)]); % extract Index of real key presses in the values
plot(PacerIndR,ones(1,length(PacerIndR)),'bo'); % look at the above Index (one press produced several indices)
% Pacer at 2Hz has intervals of 500ms (1000 samples)
% Pacer at 3Hz has interval of 333ms (666.6 samples)
% determine a threshold of numbers of samples in pacer interval
threshold = ifi*srR*0.75; % number of samples within a interval larger than the screen refresh rate
% button presses recorded via DP port into stimtracker is sampled at 2000 Hz
% feedback via photocell to stimtracker has minial interval every other waitframes (minmal of 1000 Hz)
PacerTimeIndR=PacerIndR(find([diff([0 PacerIndR'])]>threshold)); % exact index of key press onset in datatimes (reduce several indices into one)
% 2+2+30*6+2*5 = 194
% create a time series that assign botton presses as 1, all other as 0
PacerTimeR01=zeros(size(timeR));
PacerTimeR01(PacerTimeIndR)=1;
plot(PacerTimeR01,'b');
sum(PacerTimeR01) % check if == 194

% Plot not alligned
figure; % for 6 conditions
PacersL=PacerTimeIndL([1 2 3 4 5 34 35 36 37 66 67 68 69 98 99 100 101 130 131 132 133 162 163 164 165 194]);
plot(PacersL, ones(length(PacersL)),'r.');
hold on;
PacersR=PacerTimeIndR([1 2 3 4 5 34 35 36 37 66 67 68 69 98 99 100 101 130 131 132 133 162 163 164 165 194]);
plot(PacersR, ones(length(PacersR)),'b.');

figure; % for 5 conditions
PacersL=PacerTimeIndL([1 2 3 4 5 34 35 36 37 66 67 68 69 98 99 100 101 130 131 132 133 162]);
plot(PacersL, ones(length(PacersL)),'r.');
hold on;
PacersR=PacerTimeIndR([1 2 3 4 5 34 35 36 37 66 67 68 69 98 99 100 101 130 131 132 133 162]);
plot(PacersR, ones(length(PacersR)),'b.');

% make shift for missming pacers (optional, only needed for 20220517_2P)
ShiftPacer=PacerTimeIndL(end)-PacerTimeIndR(end);
PacersR=PacerTimeIndR([1 2 3 4 5 34 35 36 37 66 67 68 69 98 99 100 101 130 131 132 133 162 163 164 165 194]);
PacersL=PacersR+ShiftPacer;
figure;
plot(PacersL,ones(1,length(PacersL)),'r.');hold on;ylim([0 3]);
plot(PacerTimeIndL,2*ones(1,length(PacerTimeIndL)),'r.');

% Plot alligned
figure;
plot(PacersL-PacersL(1), ones(length(PacersL)),'r.');
hold on;
plot(PacersR-PacersR(1), ones(length(PacersR)),'b.');

%% Extract feedbacks from light detector ISO (skip for now)
% look for the second ISO aux channel for the photocell 
ISOauxindL=find(channels_infoL.labels=='ISO aux');
ISOauxindR=find(channels_infoR.labels=='ISO aux');


% examine the Left
figure;subplot(2,1,1);
plot(samplesL(ISOauxindL(1),:),'r'); title('ISOauxind(1)')%  ISO aux = analog
subplot(2,1,2);
plot(samplesL(ISOauxindL(2),:),'b'); title('ISOauxind(2)')
% select a good one
Photocell_L=samplesL(ISOauxindL(1),:)'.*-1;
Photocell_L=samplesL(ISOauxindL(2),:)';

% for photocell
% view the time course of photocell signals
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(Photocell_L);xlabel('time');ylabel('photocell signal');
% plot EEG on top
% hold on; plot(time,samples(2:33,:)'); % it zoom out the phtocell amplitude, too small to see
% click and select the start and end point for peak extraction
[x, y] = ginput(2); % read two mouse clicks on the plot % x were index, y were real values
Startpoint=round(x(1));Endpoint=round(x(2)); % Startpoint and Endpoint are sample number or index in time
hold on;xline(x(1),'r');xline(x(2),'r');hold off;

% replace the beginning and end with baseline value
Photocell(1:Startpoint)=mean(y);Photocell(Endpoint:end)=mean(y); % plot(Photocell');
plot(time,Photocell,'b'); 

% Examine peaks detection in analog1
Halfhigh1=3/4*(max(Photocell_L)-min(Photocell_L)); % value for 'MinPeakProminence'
% Check if need to adjust the Halfhigh cutoff
close;figure;
findpeaks(Photocell_L,1:length(timeL),'MinPeakProminence',Halfhigh1,'Annotate','extents');
yline(Halfhigh1,'m','MinPeakProminence');

% locate the trial sessions % pks=value of the peak % locs=time of the peak
[pksL,locsL] = findpeaks(Photocell_L,1:length(timeL),'MinPeakProminence',Halfhigh1,'Annotate','extents');

% examine pks and locs (both are values of analog1data and datatimes, not indices)
% i=1;
% find(Photocell==pks(i)) % return the index of peak in time
% find(datatimes==locs(i)) % return the same "index" in datatimes (ie, i=1, index=4226)
% so, infer from above, beacause the above "find values" for the same "index" works
% pks are values in analog1data
% locs are values in datatimes ((ie, i=1, value=4225))

% figure;plot(locs,'bo');
% figure;bar(locs);
% figure;plot(pks,'bo');
% figure;bar(pks);
% figure;plot(locs,pks,'bo');ylim([min(analog1data) max(analog1data)]);

% locsDurations=diff(locs);% time durations between peaks
% close;figure;
% plot(locsDurations,'ro');% look at distribution of these durations   
% xlabel('each peak');ylabel('locsDurations = Between-peaks time duration (ms)');
% Lowcutoff=5*mean(locsDurations);% cutoff standard of the between-peak duarations to separate between trials
% Highcutoff=180*mean(locsDurations);
% hold on; yline(Lowcutoff,'m--','lowcut');hold off; % examine the cutoff line
% hold on; yline(Highcutoff,'m--','highcut');hold off;



% examine the Right
figure;subplot(2,1,1);
plot(samplesR(ISOauxindR(1),:),'r'); title('ISOauxind(1)')%  ISO aux = analog
subplot(2,1,2);
plot(samplesR(ISOauxindR(2),:),'b'); title('ISOauxind(2)')
% select a good one
Photocell_R=samplesR(ISOauxindR(1),:)'.*-1;
Photocell_R=samplesR(ISOauxindR(2),:)';

% for photocell
% view the time course of photocell signals
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(Photocell_R);xlabel('time');ylabel('photocell signal');
% plot EEG on top
% hold on; plot(time,samples(2:33,:)'); % it zoom out the phtocell amplitude, too small to see
% click and select the start and end point for peak extraction
[x, y] = ginput(2); % read two mouse clicks on the plot % x were index, y were real values
Startpoint=round(x(1));Endpoint=round(x(2)); % Startpoint and Endpoint are sample number or index in time
hold on;xline(x(1),'r');xline(x(2),'r');hold off;

% replace the beginning and end with baseline value
Photocell(1:Startpoint)=mean(y);Photocell(Endpoint:end)=mean(y); % plot(Photocell');
plot(time,Photocell,'b'); 

% Examine peaks detection in analog1
Halfhigh1=3/4*(max(Photocell_R)-min(Photocell_R)); % value for 'MinPeakProminence'
% Check if need to adjust the Halfhigh cutoff
close;figure;
findpeaks(Photocell_R,1:length(timeR),'MinPeakProminence',Halfhigh1,'Annotate','extents');
yline(Halfhigh1,'m','MinPeakProminence');

% locate the trial sessions % pks=value of the peak % locs=time of the peak
[pksR,locsR] = findpeaks(Photocell_R,1:length(timeR),'MinPeakProminence',Halfhigh1,'Annotate','extents');

% examine pks and locs (both are values of analog1data and datatimes, not indices)
% i=1;
% find(Photocell==pks(i)) % return the index of peak in time
% find(datatimes==locs(i)) % return the same "index" in datatimes (ie, i=1, index=4226)
% so, infer from above, beacause the above "find values" for the same "index" works
% pks are values in analog1data
% locs are values in datatimes ((ie, i=1, value=4225))

% figure;plot(locs,'bo');
% figure;bar(locs);
% figure;plot(pks,'bo');
% figure;bar(pks);
% figure;plot(locs,pks,'bo');ylim([min(analog1data) max(analog1data)]);

% locsDurations=diff(locs);% time durations between peaks
% close;figure;
% plot(locsDurations,'ro');% look at distribution of these durations   
% xlabel('each peak');ylabel('locsDurations = Between-peaks time duration (ms)');
% Lowcutoff=5*mean(locsDurations);% cutoff standard of the between-peak duarations to separate between trials
% Highcutoff=180*mean(locsDurations);
% hold on; yline(Lowcutoff,'m--','lowcut');hold off; % examine the cutoff line
% hold on; yline(Highcutoff,'m--','highcut');hold off;