%% load LSL labrecorder file

% addpath  C:\Users\zhibi\Documents\GitHub\xdf-Matlab
% addpath /home/zhibin/Documents/GitHub/xdf-Matlab
% addpath D:\360MoveData\Users\alienware\Documents\GitHub\xdf-Matlab
addpath(genpath('C:\Users\zhibi\Documents\GitHub\xdf-Matlab')); 

% addpath C:\Users\NIRS\Documents\GitHub\xdf-Matlab

% cd /home/zhibin/Documents/Artinis NIRS/zhibin/20220213fNIRS_LSL
% 
% cd C:\Users\NIRS\Documents\GitHub\xdf-Matlab
% addpath C:\Users\NIRS\Documents\GitHub\xdf-Matlab
% 
% 
% cd C:\Users\NIRS\Desktop\Acquisition\zhibin
% cd F:\UCI_dataset\fNIR\20220215\2022021508
% cd F:\UCI_dataset\fNIR\20220215\
cd C:\Users\zhibi\Desktop\Artinis_NIRS\zhibin\20220223
%% load data

% Load with toolbox format
oxysoft2matlab

% load LSL streams
clear;close all;
fileList = dir('*.xdf');
fileList.name
[streams,fileheader] = load_xdf(fileList.name);

% [streams,fileheader] = load_xdf('sub-P001_ses-S001_task-Default_run-001_20220313.xdf');
% [streams,fileheader] = load_xdf('sub-P001_ses-S001_task-Default_run-001_eeg2022021503.xdf');
% [streams,fileheader] = load_xdf('202202150*.xdf');


% cfg = [];
% cfg.dataset = '2022022301.oxy5';
% [data]=ft_preprocessing(cfg);

%% load behaviral data
fileList = dir('*.mat');
fileList.name
load(fileList.name);
% load 2022021503.mat

%%
info1=streams{1,1}.info;
info1.name % OxySoft
time_series1=streams{1,1}.time_series;
time_stamps1=streams{1,1}.time_stamps;
segments1=streams{1,1}.segments;
srate=info1.effective_srate;
clock_offsets=info1.clock_offsets;

figure('units','normalized','outerposition',[0 0 0.8 0.4]);
plot(time_stamps1,time_series1');
hold on;
% plotx(time_stamps1,time_series1');
%%
info2=streams{1,2}.info;
info2.name % OxySoft Event Marker / MyMarkerStream
time_series2=streams{1,2}.time_series; % noting here
time_stamps2=streams{1,2}.time_stamps; % noting here

% plot(time_stamps2(1),2, 'bo');
% plot(time_stamps2(2),2.5, 'ro'); baseline_duration=num2str(time_stamps2(2)-time_stamps2(1));
% plot(time_stamps2(3:2:25),2, 'bo'); 
% plot(time_stamps2(4:2:26),2.5, 'ro');trials_durations=num2str(time_stamps2(4:2:26)-time_stamps2(3:2:25));

index_Baseline_Start = find(strcmp(horzcat({time_series2{:}}),'Baseline-Start')==1); % index_Baseline_Start(1)=[];
index_Baseline_End = find(strcmp(horzcat({time_series2{:}}),'Baseline-End')==1); % index_Baseline_End(1)=[];
index_Trial_Start = find(strcmp(horzcat({time_series2{:}}),'Trial-Start')==1); % index_Trial_Start(1)=[];
index_Trial_End = find(strcmp(horzcat({time_series2{:}}),'Trial-End')==1);
index_Screen_Flip = find(strcmp(horzcat({time_series2{:}}),'Screen_Flip')==1);

plot(time_stamps2(index_Baseline_Start),2, 'bo');
plot(time_stamps2(index_Baseline_End),2.5, 'ro'); 
baseline_duration=num2str(time_stamps2(index_Baseline_End(1))-time_stamps2(index_Baseline_Start(1)))
% baseline_duration=num2str(time_stamps2(index_Baseline_End(2))-time_stamps2(index_Baseline_Start(2)))
plot(time_stamps2(index_Trial_Start),2, 'bo'); 
plot(time_stamps2(index_Trial_End),2.5, 'ro');
trials_durations=num2str(time_stamps2(index_Trial_End(1:end))-time_stamps2(index_Trial_Start(1:end)))
trials_intervals=num2str(time_stamps2(index_Trial_End(1:end-1))-time_stamps2(index_Trial_Start(2:end)))


% plot(time_stamps2(index_Screen_Flip),30, 'go');
inter_trials_intervals=num2str(time_stamps2(index_Trial_Start(2:end))-time_stamps2(index_Trial_End(1:end-1)))

% plot(time_stamps2(index_Screen_Flip),3, 'go');
% 
% figure; plot(diff(time_stamps2(index_Screen_Flip)),'ro');
% diff(time_stamps2(index_Screen_Flip))

% cd D:\360MoveData\Users\alienware\Documents\GitHub\Motor_cordination\fNIR_test
cd C:\Users\zhibi\Documents\GitHub\Motor_cordination\fNIR_test
% map index in time_stamps2 to time_stamps1
[stamp1index_Baseline_Start,stamp1index_Baseline_End,stamp1index_Trial_Start,stamp1index_Trial_End] = getStamp1indices(time_stamps1,time_stamps2,time_series2);

%%
% load behaviral
open xdfchannels.m
%%
% repmat_stamp1=repmat(time_stamps1',[1,length(time_stamps2)]);
% [minValue,closestIndex_stamp1] = min(abs(repmat_stamp1-time_stamps2));
% closestValue = time_stamps1(closestIndex_stamp1) 
% 
% time_stamps2(index_Baseline_Start)
% time_stamps2(index_Baseline_End)
% time_stamps2(index_Trial_Start)
% time_stamps2(index_Trial_End)

% stamp1index_Baseline_Start
% stamp1index_Baseline_End
% stamp1index_Trial_Start
% stamp1index_Trial_End

% find(ismember(time_stamps1,time_stamps2(index_Baseline_Start)))
% find(ismember(time_stamps1,time_stamps2(index_Baseline_End)))
% find(ismember(time_stamps1,time_stamps2(index_Trial_Start)))
% find(ismember(time_stamps1,time_stamps2(index_Trial_End)))

% time_stamps2(index_Baseline_Start)
% time_stamps2(index_Baseline_End)
% time_stamps2(index_Trial_Start)
% time_stamps2(index_Trial_End)



%%
% info3=streams{1,3}.info;
% info3.name %BioSemi
% time_series3=streams{1,3}.time_series;
% time_stamps3=streams{1,3}.time_stamps;  num2str(time_stamps3(1))
% segments3=streams{1,3}.segments;
% 
% plot(time_stamps3, time_series3, 'ro')
% 
% plot(time_stamps3(1),2, 'bo');
% plot(time_stamps3(2),2.5, 'ro');
% plot(time_stamps3(3:2:25),2, 'bo');
% plot(time_stamps3(4:2:26),2.5, 'ro');
%%
% info4=streams{1,4}.info;
% info4.name % MyMarkerStream
% time_series4=streams{1,4}.time_series;
% time_stamps4=streams{1,4}.time_stamps; num2str(time_stamps4(1))
% 
% plot(time_stamps4, 2, 'bo');
% 
% figure; bar(time_series4);
% figure; scatter(time_series4);
% 
% n=1;
% num2str(time_stamps1(n))
% num2str(time_stamps3(n))
% num2str(time_stamps4(n))
% 
% num2str(time_stamps1(end))
% num2str(time_stamps3(end))
% num2str(time_stamps4(end))

