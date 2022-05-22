%% load TMSi data
clear;close all;
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing;
cd /ssd/zhibin/1overf/20220515_2P
cd /ssd/zhibin/1overf/20220517_2P
cd /ssd/zhibin/1overf/20220518_2P

[EEGfileNameL]=uigetfile('*.Poly5');% select the Left player EEG
Path_filenameL=[pwd '/' EEGfileNameL];
[timeL,samplesL,TRIGGERindL,srL,channels_infoL] = LoadTMSi(Path_filenameL);

[EEGfileNameR]=uigetfile('*.Poly5');% select the Right player EEG
Path_filenameR=[pwd '/' EEGfileNameR];
[timeR,samplesR,TRIGGERindR,srR,channels_infoR] = LoadTMSi(Path_filenameR);

%% check bottom-right photocell light detector
labelsL=channels_infoL.labels;
% look for the second ISO aux channel for the photocell 
% ISOauxind1=find(labelsL=='AUX 1-1');% ISO aux
% ISOauxind2=find(labelsL=='AUX 1-2');
ISOauxind=find(labelsL=='ISO aux')
% plot(samples(38,:),'b'); %  ISO aux = analog
% examine
figure;subplot(2,1,1);
plot(samplesL(35,:)','r'); title('ISOauxind1')%  ISO aux = analog
subplot(2,1,2);
plot(samplesL(36,:)','b'); title('ISOauxind2')
% select a good one
Photocells=samplesL(35,:);

%% Extract Trigger channel info
TriggersL=samplesL(TRIGGERindL,:)';unique(TriggersL)
TriggersR=samplesR(TRIGGERindR,:)';unique(TriggersR)
% feedback from other = 251 (255-2^2); 
% self key presses = 239 (255-2^4); 
% 235 = 255-2^2-2^4 when feedback and self key press overlapped
% pacer = 127 (255-2^7); 
% 111 = 255 -2^7-2^4 % when pacer and self key press overlapped


% find the first photocell and align Left and Right in time
TriggersL(find(TriggersL == 255)) = 0; % replace 255 with 0
indL=unique([find(TriggersL == 127);find(TriggersL == 111)]); % pacer indices
figure;plot(TriggersL,'ro'); title('Recording L (Trigger Channel)'); % xline(indL(1)); xline(indL(1));
TriggersLcut=TriggersL(indL(1):end); % remove trigger data before the first photocell;
figure;plot(TriggersLcut,'ro');% Update TriggersL: chop off the head before first photocell;

TriggersR(find(TriggersR == 255)) = 0; % replace 255 with 0
indR=find(TriggersR == 127); % photocell indices
figure;plot(TriggersR,'bo');title('Recording R (Trigger Channel)');% xline(indR(1)); xline(indR(3));
TriggersRcut=TriggersR(indR(1):end); % remove trigger data before the first photocell;
figure;plot(TriggersRcut,'bo'); % Update TriggersL: chop off the head before first photocell;

% Align with the first phtocell and examine in two subplots
figure('units','normalized','outerposition',[0 0 1 0.6]);
subplot(2,1,1);
plot(TriggersLcut,'ro');
title('Left Player, 127=stimulus, 239=self, 251=feedback');
xlim([0 max([length(TriggersLcut) length(TriggersRcut)])]); 
subplot(2,1,2);
plot(TriggersRcut,'b.');
title('Right recorder, 127=stimulus, 239=self, 251=feedback');
xlim([0 max([length(TriggersLcut) length(TriggersRcut)])]); 

% Align with the first phtocell and examine in one plot (color code with recording device)
figure('units','normalized','outerposition',[0 0 1 0.6]);
plot(TriggersLcut,'ro');
xlim([0 max([length(TriggersLcut) length(TriggersRcut)])]); 
hold on;
plot(TriggersRcut,'b.');
title('Left recorder (red), Right Player (blue), 127=stimulus, 239=self, 251=feedback');

%% Align with the first phtocell and examine in one plot (color code with player)
figure('units','normalized','outerposition',[0 0 1 0.6]);
hold on;
% Recording device on the Left
plot((TriggersL==0)*0,'wo');plot((TriggersL==127)*1,'go');plot((TriggersL==239)*2,'ro');plot((TriggersL==251)*3,'b.'); 
% Recording device on the Right
plot((TriggersR==0)*0,'wo');plot((TriggersR==127)*1,'g.');plot((TriggersR==239)*2,'b.');plot((TriggersR==251)*3,'ro');
hold off;
% Combine into one plot
legend({'Recording device L','Pacer L','Presses L','feedback from R','Recording device R','Pacer R','Presses R','feedback from L'});hold off;
title('Left Player (red), Right Player (blue), Stimulus (green); 1=Pacer, 2=buttom press, 3=feedback to the other side');
% yline(127);yline(239);yline(251);
ylim([0 60]);
xlabel('time (0.5 ms / sample)');
% xlim([0 max([length(TriggersL) length(TriggersR)])]);

% zoom in and check delay
[x, ~] = ginput(2); % read two mouse clicks on the plot % x were index, y were real values
% get the proximate index
string(x)
ind1=round(x(1))
ind2=round(x(2))
(ind1-ind2)/srL % delay in ms

%% Align and Only plot the presses (color code with player)
figure('units','normalized','outerposition',[0 0 1 0.6]);
hold on;
% Recording device on the Left
plot((TriggersL==127)*1,'ro');plot((TriggersL==239)*2,'r.');%plot((TriggersL==251)*2,'bo'); 
% Recording device on the Right
plot((TriggersR==127)*1,'b.');plot((TriggersR==239)*2,'b.');%plot((TriggersR==251)*2,'ro');
hold off;
% Combine into one plot
legend({'Pacer L','Presses L','Pacer R','Presses R'});hold off;
title('Left Player (red), Right Player (blue), Stimulus (green); 1=Pacer, 2=buttom press');
% yline(127);yline(239);yline(251);
ylim([0 60]);
xlabel('time (0.5 ms / sample)');
% xlim([0 max([length(TriggersL) length(TriggersR)])]);

% zoom in and check delay
[x, ~] = ginput(2); % read two mouse clicks on the plot % x were index, y were real values
% get the proximate index
string(x)
ind1=round(x(1))
ind2=round(x(2))
(ind1-ind2)/srL % delay in ms

%% Align and Only plot the diff (color code with player)

figure('units','normalized','outerposition',[0 0 1 0.6]);
hold on;
% Recording device on the Left
plot((TriggersL==127)*1,'ro');plot(diff((TriggersL==239)*2),'r');%plot((TriggersL==251)*2,'bo'); 
% Recording device on the Right
plot((TriggersR==127)*1,'b.');plot(diff((TriggersR==239)*2),'b');%plot((TriggersR==251)*2,'ro');
hold off;
% Combine into one plot
legend({'Pacer L','Press interval L','Pacer R','Press interval R'});hold off;
title('Left Player (red), Right Player (blue), Stimulus (green); 1=Pacer, 2=buttom press');
% yline(127);yline(239);yline(251);
ylim([0 ]);
xlabel('time (0.5 ms / sample)');
% xlim([0 max([length(TriggersL) length(TriggersR)])]);

% zoom in and check delay
[x, ~] = ginput(2); % read two mouse clicks on the plot % x were index, y were real values
% get the proximate index
string(x)
ind1=round(x(1))
ind2=round(x(2))
(ind1-ind2)/srL % delay in ms


%% plot the buttom press intervals
plot(TriggersL==239);ylim([0 2]);
BoL=(TriggersL==239);%boolean array: true values when there are button presses
pressIndL=find([0;diff(BoL)]==1); % take each button press onset
pressIntervalsL=[pressIndL(1);diff(pressIndL)]; 
figure;plot(pressIntervalsL,'r');
ylim([0 3000]);

plot(TriggersR==239);ylim([0 2]);
BoR=(TriggersR==239);%boolean array: true values when there are button presses
pressIndR=find([0;diff(BoR)]==1); % take each button press onset
pressIntervalsR=[pressIndR(1);diff(pressIndR)]; 
figure;plot(pressIntervalsR,'b');
ylim([0 3000]);

figure;
plot(pressIntervalsL,'r');
hold on;
plot(pressIntervalsR,'b');
hold off;ylim([0 3000]);

figure;
subplot(2,1,1);
plot((TriggersL==239),'r');
subplot(2,1,2);
plot((TriggersR==239),'b');
ylim([0 2]);


