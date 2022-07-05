%% load TMSi data
% open /home/zhibin/Documents/GitHub/TMSi_SAGA_matlab/SAGA_interface/testdatatryout2022.m

clear;close all;

% Step 1: Open Poly5 file.
addpath /home/zhibin/Documents/GitHub/TMSi_SAGA_matlab % hnlb
% cd /ssd/zhibin/1overf/20220324
cd /ssd/zhibin/1overf/20220331
cd /ssd/zhibin/1overf/2022041401/
% d = TMSiSAGA.Poly5.read('/ssd/zhibin/1overf/20220324/20220324.Poly5');
% d = TMSiSAGA.Poly5.read('/ssd/zhibin/1overf/20220331/RESTING_ec_2mim-20220331T161816.DATA.Poly5');
d = TMSiSAGA.Poly5.read('/ssd/zhibin/1overf/20220331/experiment_PILOT-20220331T162236.DATA.Poly5');
d = TMSiSAGA.Poly5.read('/ssd/zhibin/1overf/20220408/JACK20220408-20220408T160904.DATA.Poly5');
d = TMSiSAGA.Poly5.read('/ssd/zhibin/1overf/2022041401/2020041401Patty-20220414T145516.DATA.Poly5');
d = TMSiSAGA.Poly5.read('/ssd/zhibin/1overf/2022041402/Shane2022041402-20220414T160550.DATA.Poly5');

samples=d.samples;
sr=d.sample_rate;
channels=d.channels;

% Organize the channels

% unit_names=cell(length(channels),1);
% altertive_names=cell(length(channels),1);
% for i=1:length(channels)
%     altertive_names(i)={channels{i}.alternative_name};
%     unit_names(i)={channels{i}.unit_name};
% end

numbers=num2str([1:length(channels)]');
labels=strings(length(channels),1);
units=strings(length(channels),1);
for i=1:length(channels)
    labels(i)=channels{i}.alternative_name;
    units(i)=channels{i}.unit_name;
end
channels_info=table(numbers,labels,units)

% Create time stamps
num2str(d.time)
time=[1/sr:1/sr:d.time]';

% Plot channels of Key presses, photocells, EMG
% look for TRIGGERS channel;
TRIGGERind=find(labels=='TRIGGERS');
% plot(samples(54,:),'ro'); % trigger = key presses = 223 (255-2^5)
% plot(samples(53,:),'ro'); % trigger = key presses = 223 (255-2^5)
% plot(samples(46,:),'ro'); % trigger = key presses = 223 (255-2^5)
plot(samples(TRIGGERind,:),'ro'); % trigger = key presses = 223 (255-2^5)
% unique(samples(54,:))
% unique(samples(53,:))
% unique(samples(46,:))
unique(samples(TRIGGERind,:))
hold on; 
% look for the second ISO aux channel for the photocell 
ISOauxind=find(labels=='ISO aux');
% plot(samples(38,:),'b'); %  ISO aux = analog
% examine
figure;subplot(2,1,1);
plot(samples(ISOauxind(1),:),'r'); title('ISOauxind(1)')%  ISO aux = analog
subplot(2,1,2);
plot(samples(ISOauxind(2),:),'b'); title('ISOauxind(2)')
% select a good one
ISOauxind=ISOauxind(1);
ISOauxind=ISOauxind(2);
% look for BIP 01 channel for the EMG channel
BIP01ind=find(labels=='BIP 01');
% plot(samples(34,:),'k'); % EMG channel
% plot(samples(33,:),'k'); % EMG channel
plot(samples(BIP01ind,:),'k'); % EMG channel

% Save Channels of presses, photocells, EMG
% BottonPres=samples(54,:);
BottonPres=samples(TRIGGERind,:)';
Photocell=samples(ISOauxind,:)';
% EMG=samples(34,:);
EMG=samples(BIP01ind,:)';
EEG=samples(1:32,:)';
numChan=32;

% View EEG
% plot(time',samples(2:33,:)');
% plot(time',samples(1:32,:)');
plot(time,EEG);

%% load TMSi for EEGLAB
% open /home/zhibin/Documents/GitHub/TMSi_SAGA_matlab/SAGA_interface/testdatatryout2022.m

% cd /ssd/zhibin/1overf/20220324
cd /ssd/zhibin/1overf/20220331
cd /ssd/zhibin/1overf/
[file, pathname] = uigetfile({'*.Poly5';'*.S00';'*.TMS32'},'Pick your file');
[path,filename,extension] = fileparts(file);
% open Poly5toEEGlab

% Step 2: Plot a single channel.
plot((0:(d.num_samples - 1)) / d.sample_rate, d.samples(2, :));

% Step 3: Save dataset in the same directory as the *.Poly5 file.
% Open EEGlab
eeglab

% Load TMSi's channel location file
load('EEGChannels64TMSi.mat', 'ChanLocs');

% Transform data to eeglab. 
eegdataset = toEEGLab(d, ChanLocs);

% Save dataset in the same directory as the *.Poly5 file.
pop_saveset(eegdataset,'filename',filename,'filepath',pathname)
disp(['Data saved as EEGlab dataset (.set) in this folder: ',pathname])

%% load stimulus data
% cd /home/zhibin/Documents/GitHub/1overf/stimulus_data_storage % hnlb
% load('20220324.mat')
% load('20220331.mat')
cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/stimulus_data_storage
load('20220408.mat')

%% Covert Photocell AND  BottonPres signals to time points
% Photocell=samples(38,:);
% BottonPres=samples(54,:);
% unique(BottonPres) 

time;

% open /home/zhibin/Documents/GitHub/Motor_cordination/EEGanalysis/20211102/organize_photocells_step2.m

% for photocell
% view the time course of photocell signals
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(Photocell);xlabel('time');ylabel('photocell signal');
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
Halfhigh1=3/4*(max(Photocell)-min(Photocell)); % value for 'MinPeakProminence'
% Check if need to adjust the Halfhigh cutoff
close;figure;
findpeaks(Photocell,time,'MinPeakProminence',Halfhigh1,'Annotate','extents');
yline(Halfhigh1,'m','MinPeakProminence');

% locate the trial sessions % pks=value of the peak % locs=time of the peak
[pks,locs] = findpeaks(Photocell,time,'MinPeakProminence',Halfhigh1,'Annotate','extents');

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

% for botton presses
% botton1data=BottonPres(Startpoint:Endpoint);
botton1data=BottonPres;
% look for values other than 0 and 255
plot(time,botton1data,'bo'); % view time course of botton press values
PresInd=find(botton1data ~= 255 & botton1data ~= 0); % extract Index of real key presses in the values
plot(PresInd,ones(1,length(PresInd)),'ro'); % look at the above Index (one press produced several indices)
threshold = NumFramesInterval*ifi*sr/4; % determine a threshold of key press interval
BottonPresTimeInd=PresInd(find([0 diff(PresInd')]>threshold)); % exact index of key press onset in datatimes (reduce several indices into one)

% create a time series that assign botton presses as 1, all other as 0
BottonPresTime01=zeros(size(time));
BottonPresTime01(BottonPresTimeInd)=1;
plot(BottonPresTime01)

%% Examine photocell amd botton presses timing
PhotocellTime=locs; % locs are values in time
BottonPressTime=time(BottonPresTimeInd);

% view and compare
close;figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(PhotocellTime',ones(length(PhotocellTime),1),'bo');
hold on;
plot(BottonPressTime',ones(length(BottonPressTime),1),'r.');
xlabel('time');title('Photocell / Botton press indicator');
hold off;
legend('Photocell','Botton press');
% I ended up doing synchronization all the time


% examine key press interval
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(1:length(diff(BottonPresTimeInd)), diff(BottonPresTimeInd)/sr,'r.'); 
xlabel('key press');ylabel('Intervals (s)');
title('Differences of Time indices for botton presses');

% use a threshold to segment conditions
threshold2=max(diff(BottonPresTimeInd)/sr)/2; % automatic 
threshold2=8; % manual inspection
hold on; yline(threshold2,'m'); % show on top of previous plot
separationsInd=find((diff(BottonPresTimeInd)/sr)>threshold2); % indices of the first button press from the second condition in BottonPresTimeInd
separationsInd=[1;separationsInd+1]; %  indices from the first button press in each condition in BottonPresTimeInd
separationsTimeInd=BottonPresTimeInd(separationsInd);% indices in time 

% plot to examine
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(time,Photocell);xlabel('time');ylabel('photocell signal');
hold on; xline(time(separationsTimeInd(1)),'r');
xline(time(separationsTimeInd(2)),'r');
xline(time(separationsTimeInd(3)),'r');hold off;

%% Compute errors
open /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/1P_testing/Error_3plots.m

%% auto correlation and spectrum
open /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/1P_testing/Error_3plots.m
