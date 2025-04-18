%% load behaviral data
clear;
close all;

% add load data functions
cd /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing
addpath /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing;
addpath(genpath('/home/zhibinz2/Documents/GitHub/TMSi_SAGA_matlab'))

% data directories
cd /ssd/zhibin/1overf/20220515_2P/Segmented_data/Plots % Kenta and Jenny: Synch 2 - 3 Hz
cd /ssd/zhibin/1overf/20220517_2P/Segmented_data/Plots % Hiro and Kenta: Synch 2 - 3 Hz
cd /ssd/zhibin/1overf/20220518_2P/Segmented_data/Plots % Jenny and Hiro: Synch 2 - 3 Hz
cd /ssd/zhibin/1overf/20220609_2P/Segmented_data/Plots % Patty and Hojjat: Synch 1.3Hz 600taps
cd /ssd/zhibin/1overf/20220610_2P/Segmented_data/Plots % Patty and Hojjat: Syncopate 1.3 Hz 600taps
cd /ssd/zhibin/1overf/20220713_2P/Segmented_data/Plots % Patty and Shane: Synch 1.3Hz 230 taps
cd /ssd/zhibin/1overf/20220721_2P/Segmented_data/Plots % Patty and Hojjat: Syncopate 1.3 Hz 230 taps
cd /ssd/zhibin/1overf/20220804_2P % Wing and Mariya: Synch 1.3Hz 230 taps
cd /ssd/zhibin/1overf/20220808_2P % Wing and Mariya: Synco 1.3Hz 230 taps
cd /ssd/zhibin/1overf/20220810_2P % Jordan and Issac: Synch 1.3Hz 230 taps
cd /ssd/zhibin/1overf/20220811_2P % Jordan and Issac: Synco 1.3Hz 230 taps
cd /ssd/zhibin/1overf/20220815_2P % Ivy and Rajat: Synch 1.3Hz 230 taps
cd /ssd/zhibin/1overf/20220816_2P % Ivy and Rajat: Synco 1.3Hz 230 taps
cd /ssd/zhibin/1overf/20221003_2P % Cassie and Harlin: Synch 1.3Hz 230 taps
cd /ssd/zhibin/1overf/2022100401_2P % Cassie and Harlin: Synco 1.3Hz 230 taps
cd /ssd/zhibin/1overf/2022100402_2P % Cassie and Shiri: Synch 1.3Hz 230 taps
cd /ssd/zhibin/1overf/20221005_2P % Cassie and Shiri: Synco 1.3Hz 230 taps


load 20220713.mat
load 20220721.mat
load 2022080400.mat
load 20221003.mat

%% load TMSi data
[EEGfileNameL]=uigetfile('*.Poly5');% select the Left player EEG
Path_filenameL=[pwd '/' EEGfileNameL];
[timeL,samplesL,TRIGGERindL,srL,channels_infoL,labelsL] = LoadTMSi(Path_filenameL);

[EEGfileNameR]=uigetfile('*.Poly5');% select the Right player EEG
Path_filenameR=[pwd '/' EEGfileNameR];
[timeR,samplesR,TRIGGERindR,srR,channels_infoR,labelsR] = LoadTMSi(Path_filenameR);

%% Extract Trigger channel info
TriggersL=[];TriggersR=[];valuesL=[];valuesR=[];
TriggersL=samplesL(TRIGGERindL,:)';valuesL=unique(TriggersL)'
TriggersR=samplesR(TRIGGERindR,:)';valuesR=unique(TriggersR)'

% Check if values are as expected
([valuesL] ==[0 119 127 183 191 247 255])
([valuesR] ==[0 119 127 183 191 247 255])

% Examine Trigger signals
% Method1: this method might normalize it to the other monitor
% figure('units','normalized','outerposition',[0 0 1 0.3]); 
% Method2: get screen size and set the figure size in pixels
ScreenSize=get(0,'MonitorPositions');
FigureXpixels=ScreenSize(3);FigureYpixels=ScreenSize(4);
figure('units','pixels','position',[0 0 FigureXpixels/2 FigureYpixels/4]);
subplot(2,1,1);
plot(TriggersL,'r');
hold on;
yline(63,'k','63 = light sensor 1+2'); 
plot(find(TriggersL==63),63*ones(length(find(TriggersL==63)),1),'mo'); 
yline(192,'k','192 = psychtoolbox display starts'); 
plot(find(TriggersL==192),192*ones(length(find(TriggersL==192)),1),'mo'); 
hold off;
subplot(2,1,2);
plot(TriggersR,'b');
hold on; 
yline(192,'k','192 = psychtoolbox display starts'); 
plot(find(TriggersR==192),192*ones(length(find(TriggersR==192)),1),'mo'); 
hold off;

% troubleshooting
hold on;
yline(95,'k','95 = light sensor 2 + key 6 '); 
plot(find(TriggersR==95),95*ones(length(find(TriggersR==95)),1),'mo'); % only happen in trial 1, to be discarded anyway
yline(247,'k','247 = middle key');
plot(find(TriggersR==247),247*ones(length(find(TriggersR==247)),1),'mo'); 
yline(223,'k','223 = key 6'); % R20220810 Issac pushes a lot of key 6 in trial 1 then move to key 7 until trial 5
plot(find(TriggersR==223),223*ones(length(find(TriggersR==223)),1),'mo'); % just in trial 1, to be discarded anyway
yline(239,'k','239 = key 5'); % see how many key 5 pushes there
plot(find(TriggersR==239),239*ones(length(find(TriggersR==239)),1),'mo'); % just one additional wrong push (can be ignored)
yline(239,'k','239 = key 5'); % see how many key 5 pushes there
plot(find(TriggersL==239),239*ones(length(find(TriggersL==239)),1),'ko'); % just one additional wrong push (can be ignored)
yline(251,'k','251 = key 3'); % see how many key 3 pushes there
plot(find(TriggersL==251),251*ones(length(find(TriggersL==251)),1),'ko'); % hold off;
% Xidon2 - setting
% Line 0: Key 1  (output = 255-2^0 = 255-1 =254)
% Line 1: Key 2  (output = 255-2^1 = 255-2 =253)
% Line 2: Key 3  (output = 255-2^2 = 255-4 =251)
% Line 3: Key 4  (output = 255-2^3 = 255-8 =247) middle key
% Line 4: Key 5  (output = 255-2^4 = 255-16 =239)
% Line 5: Key 6  (output = 255-2^5 = 255-32 =223)
% Line 6: Light Sensor 1 (Back light sensor) (output = 255-2^6 = 255-64 = 191)
% Line 7: Light Sensor 2 (White light sensor) (output = 255-2^7 = 255-128 = 127)
% 119=255-2^3-2^7 % middle key + light sensor 2
% 183=255-2^3-2^6 % middle key + light sensor 1

% 192=255-63=255-32-16-8-2 % only happen when psychtoolbox display starts
% 63 = 255-192=255-(128+64) % light senosor 1 + 2
% 95 = 255-160 = 255-(128+32) % light sensor 2 + key 6
% 223 = key 6
% 239 = key 5

% For 2022100401_2P, trial 7 lost feedback from L, restrated
hold on;
% pick cut off timepints and examine
cut1_L=3110520; cut2_L=4165560;
cut1_R=3096780; cut2_R=4152250;
subplot(2,1,1);
xline(cut1_L,'k');xline(cut2_L,'k');
subplot(2,1,2);
xline(cut1_R,'k');xline(cut2_R,'k');
% cut off middle section and reconnect the data
timeL=timeL([1:cut1_L cut2_L:end]);
samplesL=samplesL(:,[1:cut1_L cut2_L:end]);
timeR=timeR([1:cut1_R cut2_R:end]);
samplesR=samplesR(:,[1:cut1_R cut2_R:end]);
% Run from the beginning of this section again 

%% extract time points for key presses (maybe create a function to do the same for pacers,feedbacks)
% Left player
% % Key presses (all 7 keys same value)
% PresIndL=unique([find(TriggersL == 239);...
%     find(TriggersL == 235);...
%     find(TriggersL == 111 )]); % extract Index of real key presses in the values
% Key presses (6 different key values)
PresIndL=unique([find(TriggersL == 247);...
    find(TriggersL == 119);...
    find(TriggersL == 183 )]); 
% for 2022100402_2P
PresIndL=unique([find(TriggersL == 247);find(TriggersL == 239);find(TriggersL == 251);... 
    find(TriggersL == 119);...
    find(TriggersL == 183 )]); 

figure('units','normalized','outerposition',[0 0 1 0.3]);
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
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(BottonPresTimeL01,'r.')
numPresL=sum(BottonPresTimeL01) % num of button presses
% 12 blocks *230 taps = 2760 (synchronization) 
% subject might ommit a few hundred taps for syncopation

% Right player
% % Key presses (all 7 keys same value)
% PresIndR=unique([find(TriggersR == 239);...
%     find(TriggersR == 235);...
%     find(TriggersR == 111 )]); % extract Index of real key presses in the values
% Key presses (6 different key values)
PresIndR=unique([find(TriggersR == 247);...
    find(TriggersR == 119);...
    find(TriggersR == 183 )]); 
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(PresIndR,ones(1,length(PresIndR)),'bo'); % look at the above Index (one press produced several indices)
% ifi=1/100;waitframes=2;sr=2000;
% determine a threshold of numbers of frames in the button press interval
threshold = 1*1.5; % interval larger than 1 consecutive samples
% button presses recorded via DP port into stimtracker is sampled at 2000 Hz
% feedback via photocell to stimtracker has minial interval every other waitframes (minmal of 1000 Hz)
BottonPresTimeIndR=PresIndR(find([0 diff(PresIndR')]>threshold)); % exact index of key press onset in datatimes (reduce several indices into one)
% create a time series that assign botton presses as 1, all other as 0
BottonPresTimeR01=zeros(size(timeR));
BottonPresTimeR01(BottonPresTimeIndR)=1;
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(BottonPresTimeR01,'b.')
numPresR=sum(BottonPresTimeR01) % num of button presses
% 12 blocks *230 taps = 2760 (synchronization)
% subject might ommit a few hundred taps for syncopation

%% extract time points for feedbacks from other
% Left Player Recording
% previous ANT_Neuro default
% FeedbIndL=unique([find(TriggersL == 251);find(TriggersL == 235)]); % extract Index of real key presses in the values
% Light senor 1(black,top)&2(white,bottom) on output line 6&7
% Line 6: Light Sensor 1 (Back light sensor) (output = 255-2^6=191)
% Line 7: Light Sensor 2 (White light sensor) (output = 255-2^7=127)
% 119=255-2^3-2^7 % middle key + light sensor 2
% 183=255-2^3-2^6 % middle key + light sensor 1
FeedbIndL=unique([find(TriggersL == 127);find(TriggersL == 119)]);
FeedbIndL=unique([find(TriggersL == 127);find(TriggersL == 119);find(TriggersL == 63)]); % for 20220808_2P

figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(FeedbIndL,ones(1,length(FeedbIndL)),'ro'); % look at the above Index (one press produced several indices)
% ifi=1/100;waitframes=2;sr=2000;
% determine a threshold of numbers of frames in the button press interval
threshold = ifi*waitframes*srL*1.5; % number of samples within a interval larger than the waitframes
% button presses recorded via DP port into stimtracker is sampled at 2000 Hz
% feedback via photocell to stimtracker has minial interval every other waitframes (minmal of 1000 Hz)
FeedbTimeIndL=FeedbIndL(find([0 diff(FeedbIndL')]>threshold)); % exact index of key press onset in datatimes (reduce several indices into one)
% create a time series that assign botton presses as 1, all other as 0
FeedbTimeL01=zeros(size(timeL));
FeedbTimeL01(FeedbTimeIndL)=1;
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(FeedbTimeL01,'r.')
numFbL=sum(FeedbTimeL01) % num of feedbacks from the other plus 4+2x5 resting photocells 
% 12 blocks *230 taps = 2760

% Right Player Recording
% previous ANT_Neuro default
% FeedbIndR=unique([find(TriggersR == 251);find(TriggersR == 235)]); % extract Index of real key presses in the values
% Light senor 1(black,top)&2(white,bottom) on output line 6&7
% Line 6: Light Sensor 1 (Back light sensor) (output = 255-2^6=191)
% Line 7: Light Sensor 2 (White light sensor) (output = 255-2^7=127)
% 119=255-2^3-2^7 % middle key + light sensor 2
% 183=255-2^3-2^6 % middle key + light sensor 1
FeedbIndR=unique([find(TriggersR == 127);find(TriggersR == 119)]);

figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(FeedbIndR,ones(1,length(FeedbIndR)),'bo'); % look at the above Index (one press produced several indices)
% ifi=1/100;waitframes=2;sr=2000;
% determine a threshold of numbers of frames in the button press interval
threshold = ifi*waitframes*srR*1.5; % number of samples within a interval larger than the waitframes
% button presses recorded via DP port into stimtracker is sampled at 2000 Hz
% feedback via photocell to stimtracker has minial interval every other waitframes (minmal of 1000 Hz)
FeedbTimeIndR=FeedbIndR(find([0 diff(FeedbIndR')]>threshold)); % exact index of key press onset in datatimes (reduce several indices into one)
% create a time series that assign botton presses as 1, all other as 0
FeedbTimeR01=zeros(size(timeR));
FeedbTimeR01(FeedbTimeIndR)=1;
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(FeedbTimeR01,'b.')
numFbR=sum(FeedbTimeR01) % num of feedbacks from the other plus 4+2x5 resting photocells
% 12 blocks *230 taps = 2760

%% extract time points for pacers - (using this as index to segment data into conditions)
% Left player
% previous ANT_Neuro default
% PacerIndL=unique([find(TriggersL == 127)]); % extract Index of real key presses in the values
% PacerIndL=unique([find(TriggersL == 127);find(TriggersL == 111)); % extract Index of real key presses in the values
% PacerIndL=unique([find(TriggersL == 127);find(TriggersL == 111);find(TriggersL == 123)]); % extract Index of real key presses in the values
% Light senor 1(black,top)&2(white,bottom) on output line 6&7
% Line 6: Light Sensor 1 (Back light sensor) (output = 255-2^6=191)
% Line 7: Light Sensor 2 (White light sensor) (output = 255-2^7=127)
% 119=255-2^3-2^7 % middle key + light sensor 2
% 183=255-2^3-2^6 % middle key + light sensor 1
PacerIndL=unique([find(TriggersL == 191);find(TriggersL == 183)]);
PacerIndL=unique([find(TriggersL == 191);find(TriggersL == 183);find(TriggersL == 63)]); % for 20220808_2P

% figure('units','normalized','outerposition',[0 0 1 0.3]);
canvas(0.5,0.125);
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
% figure('units','normalized','outerposition',[0 0 1 0.3]);
canvas(0.5,0.125);
plot(PacerTimeL01,'r');
numPacerL=sum(PacerTimeL01) % check if == 194 / 386 for synchronization % 4+2*11+30*3+230*9=2186 for syncopation
% xlim([-1 length(PacerTimeL01)]);

% Right player
% previous ANT_Neuro default
% PacerIndR=unique([find(TriggersR == 127);find(TriggersR == 111)]); % extract Index of real key presses in the values
% Light senor 1(black,top)&2(white,bottom) on output line 6&7
% Line 6: Light Sensor 1 (Back light sensor) (output = 255-2^6=191)
% Line 7: Light Sensor 2 (White light sensor) (output = 255-2^7=127)
% 119=255-2^3-2^7 % middle key + light sensor 2
% 183=255-2^3-2^6 % middle key + light sensor 1
PacerIndR=unique([find(TriggersR == 191);find(TriggersR == 183)]);
% For 20220810_2P, it counted 517. (30*17=510 + 2 + 5)

% figure('units','normalized','outerposition',[0 0 1 0.3]);
canvas(0.5,0.125);
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
% figure('units','normalized','outerposition',[0 0 1 0.3]);
canvas(0.5,0.125);
plot(PacerTimeR01,'b');
numPacerR=sum(PacerTimeR01) % check if == 194 / 386 for synchronization % 4+2*11+30*3+230*9=2186 for syncopation
% For 20220810_2P, it counted 552. (30*17=510 + 2*15 + 4 + 8)
%% Select segment points  (synchronization)
% % Plot selected segment points (not alligned)
% figure('units','normalized','outerposition',[0 0 1 0.3]); % for 6 conditions (synchronization)
% PacersL=PacerTimeIndL([1 2 3 4 5 34 35 36 37 66 67 68 69 98 99 100 101 130 131 132 133 162 163 164 165 194]);
% plot(PacersL, ones(length(PacersL)),'r.');
% hold on;
% PacersR=PacerTimeIndR([1 2 3 4 5 34 35 36 37 66 67 68 69 98 99 100 101 130 131 132 133 162 163 164 165 194]);
% plot(PacersR, ones(length(PacersR)),'b.');

% figure('units','normalized','outerposition',[0 0 1 0.3]); % for 5 conditions (synchronization)
% PacersL=PacerTimeIndL([1 2 3 4 5 34 35 36 37 66 67 68 69 98 99 100 101 130 131 132 133 162]);
% plot(PacersL, ones(length(PacersL)),'r.');
% hold on;
% PacersR=PacerTimeIndR([1 2 3 4 5 34 35 36 37 66 67 68 69 98 99 100 101 130 131 132 133 162]);
% plot(PacersR, ones(length(PacersR)),'b.');

figure('units','normalized','outerposition',[0 0 1 0.3]); % for 12 blocks (synchronization)
PacersL=PacerTimeIndL([1:2 repelem([2+[32*([1:12]-1)]],4)+repmat([1 2 3 32],1,12)]);
plot(PacersL, ones(length(PacersL)),'r.');
hold on;
PacersR=PacerTimeIndR([1:2 repelem([2+[32*([1:12]-1)]],4)+repmat([1 2 3 32],1,12)]);
plot(PacersR, 0.5*ones(length(PacersR)),'b.');
ylim([0 1.5])

%% make shift for missming pacers (optional) 
% examine first and then pick one pacer from both recordings to calculate shift
canvas(0.5,0.25);
plot(PacerTimeIndL,1.1*ones(1,length(PacerTimeIndL)),'ro');hold on;
plot(PacerTimeIndR,ones(1,length(PacerTimeIndR)),'bo');
ylim([0 3]);

% Make up for missing pacers on the left (for 20220517_2P) (synchronization)
ShiftPacer=PacerTimeIndL(end)-PacerTimeIndR(end);
PacersR=PacerTimeIndR([1 2 3 4 5 34 35 36 37 66 67 68 69 98 99 100 101 130 131 132 133 162 163 164 165 194]);
PacersL=PacersR+ShiftPacer;
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(PacersL,ones(1,length(PacersL)),'r.');hold on;ylim([0 3]);
plot(PacerTimeIndL,2*ones(1,length(PacerTimeIndL)),'r.');

% Make up for missing pacers on the left (for 20220713_2P) (synchronization 12 blocks)
ShiftPacer=PacerTimeIndL(end)-PacerTimeIndR(end);
PacerTimeIndL=PacerTimeIndR+ShiftPacer;  % save("20220713workspace.mat","PacerTimeIndL",'-append');
PacersR=PacerTimeIndR([1:2 repelem([2+[32*([1:12]-1)]],4)+repmat([1 2 3 32],1,12)]);
PacersL=PacersR+ShiftPacer;
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(PacersL,ones(1,length(PacersL)),'r.');hold on;ylim([0 3]);
plot(PacerTimeIndL,2*ones(1,length(PacerTimeIndL)),'r.');

% Plot alligned (synchronization)
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(PacersL-PacersL(1), ones(length(PacersL)),'r.');
hold on;
plot(PacersR-PacersR(1), ones(length(PacersR)),'b.');

% Plot alligned (for syncopation 20220610_2P)
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(PacerTimeIndL-PacerTimeIndL(1), ones(length(PacerTimeIndL)),'r.');
hold on;
plot(PacerTimeIndR-PacerTimeIndR(1), ones(length(PacerTimeIndR)),'b.');

% Plot alligned (for syncopation 20220721_2P)
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(PacerTimeIndL-PacerTimeIndL(1), ones(length(PacerTimeIndL)),'r.');
hold on;
plot(PacerTimeIndR-PacerTimeIndR(2), 0.9*ones(length(PacerTimeIndR)),'b.');
ylim([0 2]);
legend('alligned-PacerTimeIndL', 'alligned-PacerTimeIndR');
% Make up for missing pacers on the left (for 20220721_2P) (sycopation 12 blocks)
ShiftPacer=PacerTimeIndL(1)-PacerTimeIndR(2);
PacerTimeIndL=PacerTimeIndR+ShiftPacer;
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(PacerTimeIndL,1.6*ones(1,length(PacerTimeIndL)),'r.');hold on;ylim([0 3]);
plot(PacerTimeIndR,1.4*ones(1,length(PacerTimeIndR)),'b.');
legend('restored PacerTimeIndL', 'PacerTimeIndR');
% Plot alligned again (for syncopation 20220721_2P)
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(PacerTimeIndL-PacerTimeIndL(1), ones(length(PacerTimeIndL)),'r.');
hold on;
plot(PacerTimeIndR-PacerTimeIndR(1), 0.9*ones(length(PacerTimeIndR)),'b.');
ylim([0 2]);
legend('alligned-PacerTimeIndL', 'alligned-PacerTimeIndR');

% plot pacers alligned (for 20220804_2P and 20220808_2P)
PacerTimeIndL;% might contain missing photocells, need to be recovered from R, should have about 30*12
PacerTimeIndR;% should have 4+2*11+30*12+4+2+30*2+2+1=455, alligned to last photocell in L
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(PacerTimeIndL-PacerTimeIndL(end), ones(1,length(PacerTimeIndL)),'r.');
hold on;
plot(PacerTimeIndR-PacerTimeIndR(end), 0.9*ones(1,length(PacerTimeIndR)),'b.');
ylim([0 2]);
legend('alligned-PacerTimeIndL', 'alligned-PacerTimeIndR');
% Make up for missing pacers on the left (for 20220804_2P and 20220808_2P)
ShiftPacer=PacerTimeIndL(end)-PacerTimeIndR(end);
PacerTimeIndL=PacerTimeIndR+ShiftPacer;
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(PacerTimeIndL,1.6*ones(1,length(PacerTimeIndL)),'r.');hold on;ylim([0 3]);
plot(PacerTimeIndR,1.4*ones(1,length(PacerTimeIndR)),'b.');
xline(0,'y');
legend('restored PacerTimeIndL', 'PacerTimeIndR');
% Plot alligned again (for 20220804_2P and 20220808_2P)
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(PacerTimeIndL-PacerTimeIndL(1), ones(length(PacerTimeIndL)),'r.');
hold on;
plot(PacerTimeIndR-PacerTimeIndR(1), 0.9*ones(length(PacerTimeIndR)),'b.');
ylim([0 2]);
legend('alligned-PacerTimeIndL', 'alligned-PacerTimeIndR');

% plot pacers alligned (for synchronization 20220810_2P)
PacerTimeIndL;% contain missing photocells, need to be recovered from R, should have about 30*12
PacerTimeIndR;% should have 4+2*4+30*5+4+4+30*12+2*11=552, alligned to last photocell in L
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot([PacerTimeIndL-PacerTimeIndL(end)], ones(1,length(PacerTimeIndL)),'r.');
hold on;
plot([PacerTimeIndR-PacerTimeIndR(end)], 0.9*ones(1,length(PacerTimeIndR)),'b.');
ylim([0 2]);
legend('alligned-PacerTimeIndL', 'alligned-PacerTimeIndR');
% Make up for missing pacers on the left (for synchronization 20220810_2P)
ShiftPacer=PacerTimeIndL(end)-PacerTimeIndR(end);
PacerTimeIndL=PacerTimeIndR+ShiftPacer;
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(PacerTimeIndL,1.6*ones(1,length(PacerTimeIndL)),'r.');hold on;ylim([0 3]);
plot(PacerTimeIndR,1.4*ones(1,length(PacerTimeIndR)),'b.');
xline(0,'y');
legend('restored PacerTimeIndL', 'PacerTimeIndR');
% Plot alligned again (for synchronization 20220810_2P)
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(PacerTimeIndL-PacerTimeIndL(1), ones(length(PacerTimeIndL)),'r.');
hold on;
plot(PacerTimeIndR-PacerTimeIndR(1), 0.9*ones(length(PacerTimeIndR)),'b.');
ylim([0 2]);
legend('alligned-PacerTimeIndL', 'alligned-PacerTimeIndR');


%% examine time points for pacers with segmentation purpose - synchronization -20220804_2P (12 blocks)
validIndL=find(PacerTimeIndL>0); % first pacer start from 67
validPacerTimeIndL=PacerTimeIndL(validIndL(1):end); % start from frist valid pacer 
validPacerTimeIndR=PacerTimeIndR(validIndL(1):end); % first corresponding pacer start from 67
% should both have 2*10+30*10+4+1+2+30+2+30+1=390 elements
% but both have 389, missing one photocell for the resting interval start before trial 8 (no idea why)
% make up for that one missing photocell 
value_diff=validPacerTimeIndL(2)-validPacerTimeIndL(1) % sample length of the first resting interval =60003
validPacerTimeIndL(34)-validPacerTimeIndL(33) % sample length of the 2nd resting interval =60003
validPacerTimeIndL(66)-validPacerTimeIndL(65) % sample length of the 3rd resting interval =60002
(validPacerTimeIndL(2)-validPacerTimeIndL(1))/2000 % varify that it is about 30 seconds
Ind_missing=2*(7-2)+30*(7-2)+1 % index for the resting itnerval end before trial 8
value_missingL=validPacerTimeIndL(Ind_missing)-value_diff; % value for the resting itnerval start before trial 8
value_missingR=validPacerTimeIndR(Ind_missing)-value_diff; % value for the resting itnerval start before trial 8
finalPacerTimeIndL=[validPacerTimeIndL(1:Ind_missing-1);value_missingL;validPacerTimeIndL(Ind_missing:end)];  % now 390 elements
finalPacerTimeIndR=[validPacerTimeIndR(1:Ind_missing-1);value_missingR;validPacerTimeIndR(Ind_missing:end)];  % now 390 elements
% plot final PacerTimeInd to make sure they look right
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(finalPacerTimeIndL,1.6*ones(1,length(finalPacerTimeIndL)),'r.');hold on;ylim([0 3]);
plot(finalPacerTimeIndR,1.4*ones(1,length(finalPacerTimeIndR)),'b.');
xline(0,'y');
legend('final PacerTimeIndL', 'final PacerTimeIndR');

%% examine time points for pacers with segmentation purpose - synchronization -20220810_2P (12 blocks)
% valid pacers are first 4 + first pacer for trial 1 (start from 4+5*30+4*2+4+4+1=171) to the end
validPacerTimeIndL=PacerTimeIndL([1:4 171:end]); % start from frist valid pacer 
validPacerTimeIndR=PacerTimeIndR([1:4 171:end]); % first corresponding pacer start from 67
% should both have 4+30*12+11*2=386 elements
% value for the resting itnerval start before trial 8
finalPacerTimeIndL=validPacerTimeIndL;  % now 386 elements
finalPacerTimeIndR=validPacerTimeIndR;  % now 386 elements
% plot final PacerTimeInd to make sure they look right
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(finalPacerTimeIndL,1.6*ones(1,length(finalPacerTimeIndL)),'r.');hold on;ylim([0 3]);
plot(finalPacerTimeIndR,1.4*ones(1,length(finalPacerTimeIndR)),'b.');
xline(0,'y');
legend('final PacerTimeIndL', 'final PacerTimeIndR');

%% Select segment points (20220804_2P)
figure('units','normalized','outerposition',[0 0 1 0.3]); % for 12 blocks (synchronization)
PacersL=finalPacerTimeIndL([[repelem([32*([1:10]-1)],4)+repmat([1 2 3 32],1,10)] 321:324 328 (328+29):(328+29+2+1) (360+29)]);
plot(PacersL, ones(length(PacersL)),'r.');
hold on;
PacersR=finalPacerTimeIndR([[repelem([32*([1:10]-1)],4)+repmat([1 2 3 32],1,10)] 321:324 328 (328+29):(328+29+2+1) (360+29)]);
plot(PacersR, 1.1*ones(length(PacersR)),'b.');
ylim([0 2]);

%% Select segment points (20220810_2P)
figure('units','normalized','outerposition',[0 0 1 0.3]); % for 12 blocks (synchronization)
PacersL=finalPacerTimeIndL([1:2 [repelem([32*([1:12]-1)],4)+repmat([3 4 5 34],1,12)] ]); % should have 4+2*12+2*11=50 elements
plot(PacersL, ones(length(PacersL)),'r.');
hold on;
PacersR=finalPacerTimeIndR([1:2 [repelem([32*([1:12]-1)],4)+repmat([3 4 5 34],1,12)] ]); % should have 4+2*12+2*11=50 elements
plot(PacersR, 1.1*ones(length(PacersR)),'b.');
ylim([0 2]);


%% extract time points for pacers with segmentation purpose - syncopation -20220610_2P (4 blocks)
% PacerTimeIndL; PacerTimeIndR; % index in datatime
% PacerTimeL01; PacerTimeR01;
% figure; plot(PacerTimeL01,'b'); figure; plot(PacerTimeR01,'b');
% L player pacer
thresholdL=2000*3; % after longer than 3 sec interval
SegPacerIndL=[1 find([0 diff(PacerTimeIndL')]>thresholdL)]; % should be 4+3*2+4*1=14 values for 4 blocks (first five should be 1 2 3 4 5)
SegtimeIndL=PacerTimeIndL(SegPacerIndL); % index in data time
% R player pacer
thresholdR=2000*3; % after longer than 3 sec interval
SegPacerIndR=[1 find([0 diff(PacerTimeIndR')]>thresholdR)]; % should be 14 values (first five should be 1 2 3 4 5)
SegtimeIndR=PacerTimeIndR(SegPacerIndR); % index in data time

% plot to examine
ScreenSize=get(0,'MonitorPositions');
FigureXpixels=ScreenSize(3);FigureYpixels=ScreenSize(4);
figure('units','pixels','position',[0 0 FigureXpixels/2 FigureYpixels/4]);
subplot(2,1,1);
plot(SegtimeIndL,ones(1,length(SegtimeIndL)),'ro');
subplot(2,1,2);
plot(SegtimeIndR,ones(1,length(SegtimeIndR)),'bo');

%% extract time points for pacers with segmentation purpose - syncopation -(12 blocks)
% Work for 20220721_2P; 20220816_2P;2022100401_2P;20221005_2P;

% L player pacer
% % thresholdL=2000*3; % after longer than 3 sec interval; for 20220721_2P;
thresholdL=2000*2.5; % after longer than 2.5 sec interval; for 20220816_2P;
SegPacerIndL=[1 find([0 diff(PacerTimeIndL')]>thresholdL)]; % should be 4+2*11+12=38 values (first five should be 1 2 3 4 5)
SegPacerIndL(1)=[];% remove the first one for 2022100401_2P
SegtimeIndL=PacerTimeIndL(SegPacerIndL); % index in data time
% R player pacer
% thresholdR=2000*3; % after longer than 3 sec interval; for 20220721_2P;
thresholdR=2000*2.5; % after longer than 3 sec interval; for 20220816_2P;
SegPacerIndR=[1 find([0 diff(PacerTimeIndR')]>thresholdR)]; % should be 38 values (first five should be 1 2 3 4 5)
SegPacerIndR(1)=[];% remove the first one for 2022100401_2P
SegtimeIndR=PacerTimeIndR(SegPacerIndR); % index in data time

% plot to examine
ScreenSize=get(0,'MonitorPositions');
FigureXpixels=ScreenSize(3);FigureYpixels=ScreenSize(4);
figure('units','pixels','position',[0 0 FigureXpixels/2 FigureYpixels/4]);
plot(PacerTimeIndL-PacerTimeIndL(1),3*ones(1,length(PacerTimeIndL)),'ro');hold on;
plot(SegtimeIndL-PacerTimeIndL(1),3.1*ones(1,length(SegtimeIndL)),'ro');
plot(PacerTimeIndR-PacerTimeIndR(1),2*ones(1,length(PacerTimeIndR)),'bo');hold on;
plot(SegtimeIndR-PacerTimeIndR(1),2.1*ones(1,length(SegtimeIndR)),'bo');
ylim([1 4])

%% extract time points for pacers with segmentation purpose - syncopation -(12 blocks)(20220808_2P)
% Work for 20220808_2P;

% L player pacer
thresholdL=2000*3; % after longer than 6 sec interval; for 20220808_2P;
SegPacerIndL=[1 find([0 diff(PacerTimeIndL')]>thresholdL)]; % should be 2*11+12+1+3=38 values
SegtimeIndL=PacerTimeIndL(SegPacerIndL); % index in data time
% R player pacer
thresholdR=2000*3; % after longer than 3 sec interval; for 20220808_2P;
SegPacerIndR=[1 find([0 diff(PacerTimeIndR')]>thresholdR)]; % should be 2*11+12+1+3=38 values
SegtimeIndR=PacerTimeIndR(SegPacerIndR); % index in data time

% plot to examine
ScreenSize=get(0,'MonitorPositions');
FigureXpixels=ScreenSize(3);FigureYpixels=ScreenSize(4);
figure('units','pixels','position',[0 0 FigureXpixels/2 FigureYpixels/4]);
subplot(2,1,1);
plot(PacerTimeIndL,ones(1,length(PacerTimeIndL)),'ro');
hold on;
plot(SegtimeIndL,1.1*ones(1,length(SegtimeIndL)),'ro');
ylim([0 2]);
subplot(2,1,2);
plot(PacerTimeIndR,ones(1,length(PacerTimeIndR)),'bo');
hold on;
plot(SegtimeIndR,1.1*ones(1,length(SegtimeIndR)),'bo');
ylim([0 2]);

% remove the 17th and 36th pacer from both sides
SegPacerIndL([17 36])=[];
SegPacerIndR([17 36])=[];
SegtimeIndL([17 36])=[];
SegtimeIndR([17 36])=[];

% plot to examine again 
% the last two pacers mark the start and end of resting eye eeg period 

%% Extract feedbacks from light detector ISO (skip for now)
% look for the second ISO aux channel for the photocell 
ISOauxindL=find(channels_infoL.labels=='ISO aux');
ISOauxindR=find(channels_infoR.labels=='ISO aux');

% examine the Left
figure('units','normalized','outerposition',[0 0 1 0.3]);
subplot(2,1,1);
plot(samplesL(ISOauxindL(1),:),'r'); title('ISOauxind(1)')%  ISO aux = analog
subplot(2,1,2);
plot(samplesL(ISOauxindL(2),:),'b'); title('ISOauxind(2)')
% select a good one
Photocell_L=samplesL(ISOauxindL(1),:)';
Photocell_L=samplesL(ISOauxindL(1),:)'.*-1;
Photocell_L=samplesL(ISOauxindL(2),:)';

% for photocell
% view the time course of photocell signals
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(Photocell_L,'r');xlabel('time');ylabel('photocell signal');
% plot EEG on top
% hold on; plot(time,samples(2:33,:)'); % it zoom out the phtocell amplitude, too small to see

% click and select the start and end point for peak extraction (optional)
[x, y] = ginput(2); % read two mouse clicks on the plot % x were index, y were real values
Startpoint=round(x(1));Endpoint=round(x(2)); % Startpoint and Endpoint are sample number or index in time
hold on;xline(x(1),'r');xline(x(2),'r');hold off;

% replace the beginning and end with baseline value (optional)
Photocell(1:Startpoint)=mean(y);Photocell(Endpoint:end)=mean(y); % plot(Photocell');
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(time,Photocell,'b'); 

% Examine peaks detection in analog1
Halfhigh1=3/4*(max(Photocell_L)-min(Photocell_L)); % value for 'MinPeakProminence'
% Check if need to adjust the Halfhigh cutoff
close;figure('units','normalized','outerposition',[0 0 1 0.3]);
findpeaks(Photocell_L,1:length(timeL),'MinPeakProminence',Halfhigh1,'Annotate','extents');
yline(Halfhigh1,'m','MinPeakProminence');

% locate the trial sessions % pks=value of the peak % locs=time of the peak
[pksL,locsL] = findpeaks(Photocell_L,1:length(timeL),'MinPeakProminence',Halfhigh1,'Annotate','extents');
% 12 blocks *230 taps = 2760

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
figure('units','normalized','outerposition',[0 0 1 0.3]);
subplot(2,1,1);
plot(samplesR(ISOauxindR(1),:),'r'); title('ISOauxind(1)')%  ISO aux = analog
subplot(2,1,2);
plot(samplesR(ISOauxindR(2),:),'b'); title('ISOauxind(2)')
% select a good one
Photocell_R=samplesR(ISOauxindR(1),:)';
Photocell_R=samplesR(ISOauxindR(1),:)'.*-1;
Photocell_R=samplesR(ISOauxindR(2),:)';

% for photocell
% view the time course of photocell signals
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(Photocell_R);xlabel('time');ylabel('photocell signal');
% plot EEG on top
% hold on; plot(time,samples(2:33,:)'); % it zoom out the phtocell amplitude, too small to see

% click and select the start and end point for peak extraction (optional)
[x, y] = ginput(2); % read two mouse clicks on the plot % x were index, y were real values
Startpoint=round(x(1));Endpoint=round(x(2)); % Startpoint and Endpoint are sample number or index in time
hold on;xline(x(1),'r');xline(x(2),'r');hold off;

% replace the beginning and end with baseline value (optional)
Photocell(1:Startpoint)=mean(y);Photocell(Endpoint:end)=mean(y); % plot(Photocell');
plot(time,Photocell,'b'); 

% Examine peaks detection in analog1
Halfhigh1=3/4*(max(Photocell_R)-min(Photocell_R)); % value for 'MinPeakProminence'
% Check if need to adjust the Halfhigh cutoff
close;figure('units','normalized','outerposition',[0 0 1 0.3]);
findpeaks(Photocell_R,1:length(timeR),'MinPeakProminence',Halfhigh1,'Annotate','extents');
yline(Halfhigh1,'m','MinPeakProminence');

% locate the trial sessions % pks=value of the peak % locs=time of the peak
[pksR,locsR] = findpeaks(Photocell_R,1:length(timeR),'MinPeakProminence',Halfhigh1,'Annotate','extents');
% 12 blocks *230 taps = 2760

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


%% Organize/Save all extracted variables
% extracted experimental data
allPerm;conditionNames;
ifi;

% extracted recording variables
timeL;samplesL;TRIGGERindL;srL;channels_infoL;labelsL;
timeR;samplesR;TRIGGERindR;srR;channels_infoR;labelsR;

% extracted behaviral variables
% BP from trigger chan (Response Pad)
BottonPresTimeIndL;BottonPresTimeL01;numPresL;
BottonPresTimeIndR;BottonPresTimeR01;numPresR;
% FB from trigger chan (light sensor 2)
FeedbTimeIndL;FeedbTimeL01;numFbL;
FeedbTimeIndR;FeedbTimeR01;numFbR;
% Pacers and Makers (light sensor 1)
PacerTimeIndL;PacerTimeL01;numPacerL;
PacerTimeIndR;PacerTimeR01;numPacerR;
% FB from light detector (ISO aux chan)
locsL; % Index in datatime; should be same as FeedbTimeIndL
locsR; % Index in datatime; should be same as FeedbTimeIndR

%% Compare FB from light sensor 2 and from ISO aux
close;figure('units','normalized','outerposition',[0 0 1 0.6]);
subplot(2,1,1);
plot(FeedbTimeIndL,ones(length(FeedbTimeIndL),1),'r.');
hold on;
plot(locsL,ones(length(locsL),1),'g.');
plot(locsL,ones(length(locsL),1).*1.1,'g.');
plot(FeedbTimeIndL,ones(length(FeedbTimeIndL),1).*1.1,'r.');
ylim([0 2]);legend('stimstracker','detector','detector','stimstracker','Location','southeast');
subplot(2,1,2);
plot(FeedbTimeIndR,ones(length(FeedbTimeIndR),1),'r.');
hold on;
plot(locsR,ones(length(locsR),1),'g.');
plot(locsR,ones(length(locsR),1).*1.1,'g.');
plot(FeedbTimeIndR,ones(length(FeedbTimeIndR),1).*1.1,'r.');
ylim([0 2]);legend('stimstracker','detector','detector','stimstracker','Location','southeast');


