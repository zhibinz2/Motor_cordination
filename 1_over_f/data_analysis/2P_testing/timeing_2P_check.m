%% load TMSi data
clear;close all;

% Step 1: Open Poly5 file.
addpath(genpath('/home/zhibin/Documents/GitHub/TMSi_SAGA_matlab')); % hnlb
cd /ssd/zhibin/1overf/20220504_2P
% dL = TMSiSAGA.Poly5.read([pwd '/Left20220504-20220504T143810.DATA.Poly5']);
% samplesL=dL.samples;
% srL=dL.sample_rate;
% channelsL=dL.channels;

d = TMSiSAGA.Poly5.read([pwd '/Right20220504-20220504T143808.DATA.Poly5']);
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
num2str(d.time)
time=[1/sr:1/sr:d.time]';

% Plot channels of Key presses, photocells, EMG
% look for TRIGGERS channel;
TRIGGERind=find(labels=='TRIGGERS');
plot(samples(TRIGGERind,:),'ro');
unique(samples(TRIGGERind,:))
% feedback from other = 251 (255-2^2); 
% self key presses = 239 (255-2^4); 
% stimulus photocell = 127 (255-2^7); 





