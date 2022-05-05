% %% load TMSi data
clear;close all;

%  Step 1: load Poly5 file.
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing;
open timeing_2P_check.m
cd /ssd/zhibin/1overf/20220504_2P
% filenames=dir; 
[EEGfileNameL]=uigetfile('*.Poly5');% select the Left player EEG
Path_filenameL=[pwd '/' EEGfileNameL];
[timeL,samplesL,TRIGGERindL,srL,channels_infoL] = LoadTMSi(Path_filenameL);
[EEGfileNameR]=uigetfile('*.Poly5');% select the Right player EEG
Path_filenameR=[pwd '/' EEGfileNameR];
[timeR,samplesR,TRIGGERindR,srR,channels_infoR] = LoadTMSi(Path_filenameR);

% Extract Trigger info
TriggersL=samplesL(TRIGGERindL,:)';unique(TriggersL)
TriggersR=samplesR(TRIGGERindR,:)';unique(TriggersR)
% feedback from other = 251 (255-2^2); 
% self key presses = 239 (255-2^4); 
% stimulus photocell = 127 (255-2^7); 
figure;plot(TriggersL,'ro');hold on;plot(TriggersR,'bo');hold off; % examine (not aligned in time);

% find the first photocell and align Left and Right in time
TriggersL(find(TriggersL== 255)) = 0; % replace 255 with 0
indL=find(TriggersL == 127); % photocell indices
plot(TriggersL,'ro');xline(indL(1));
TriggersL=TriggersL(indL(1):end);plot(TriggersL,'ro');% Update TriggersL: chop off the head before first photocell;
TriggersR(find(TriggersR== 255)) = 0; % replace 255 with 0
indR=find(TriggersR == 127); % photocell indices
plot(TriggersR,'bo');xline(indR(1));
TriggersR=TriggersR(indR(1):end);plot(TriggersR,'bo'); % Update TriggersL: chop off the head before first photocell;

% align with the first phtocell and examine
figure('units','normalized','outerposition',[0 0 1 0.6]);
subplot(2,1,1);
plot(TriggersL,'ro');
title('Left Player, 127=stimulus, 239=self, 251=feedback');
xlim([0 3e5]); 
subplot(2,1,2);
plot(TriggersR,'bo');
title('Right Player, 127=stimulus, 239=self, 251=feedback');
xlim([0 3e5]);

% align with the first phtocell and examine
figure('units','normalized','outerposition',[0 0 1 0.6]);
plot(TriggersL,'ro');
xlim([0 3e5]); 
hold on;
plot(TriggersR,'b.');
title('Left Player (red), Right Player (blue), 127=stimulus, 239=self, 251=feedback');

% align with the first phtocell and examine
figure('units','normalized','outerposition',[0 0 1 0.6]);
hold on; plot((TriggersL==255)*255,'ko'); plot((TriggersL==127)*127,'go');plot((TriggersL==239)*239,'ro');plot((TriggersL==251)*251,'bo');hold off;
xlim([0 3e5]); 
hold on; plot((TriggersR==255)*255,'ko'); plot((TriggersR==127)*127,'go');plot((TriggersR==239)*239,'bo');plot((TriggersR==251)*251,'ro');hold off;
title('Left Player (red), Right Player (blue), Stimulus (green); 127=stimulus, 239=buttom press, 251=feedback to the other side');
yline(127);yline(239);yline(251);
xlabel('time (0.5 ms / sample)');

% zoom in and check delay
[x, ~] = ginput(2); % read two mouse clicks on the plot % x were index, y were real values
% get the proximate index
string(x)
ind1=round(x(1))
ind2=round(x(2))
(ind1-ind2)/srL % delay in second


%% Local function
% https://www.mathworks.com/help/matlab/matlab_prog/local-functions-in-scripts.html
function [time,samples,TRIGGERind,sr,channels_info] = LoadTMSi(Path_filename);
    % d = TMSiSAGA.Poly5.read([pwd '/Right20220504-20220504T143808.DATA.Poly5']);
    d = TMSiSAGA.Poly5.read(Path_filename);

    samples=d.samples;
    sr=d.sample_rate;

    channels=d.channels;
    numbers=num2str([1:length(channels)]');
    labels=strings(length(channels),1);
    units=strings(length(channels),1);
    for i=1:length(channels)
        labels(i)=channels{i}.alternative_name;
        units(i)=channels{i}.unit_name;
    end
    channels_info=table(numbers,labels,units)

    % Create time stamps
    % num2str(d.time)
    time=[1/sr:1/sr:d.time]';

    % Plot channels of Key presses, photocells, EMG
    % look for TRIGGERS channel;
    TRIGGERind=find(labels=='TRIGGERS');
    % plot(samples(TRIGGERind,:),'ro');
    % unique(samples(TRIGGERind,:))
    % feedback from other = 251 (255-2^2); 
    % self key presses = 239 (255-2^4); 
    % stimulus photocell = 127 (255-2^7); 

    % BottonPres=samples(TRIGGERind,:)';
    % EEG=samples(1:32,:)';
end





