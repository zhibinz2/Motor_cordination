%% load LSL labrecorder file

% addpath  C:\Users\zhibi\Documents\GitHub\xdf-Matlab
addpath /home/zhibin/Documents/GitHub/xdf-Matlab
cd /home/zhibin/Documents/Artinis NIRS/zhibin/20220213fNIRS_LSL

cd C:\Users\NIRS\Documents\GitHub\xdf-Matlab
addpath C:\Users\NIRS\Documents\GitHub\xdf-Matlab


cd C:\Users\NIRS\Desktop\Acquisition\zhibin

%%

[streams,fileheader] = load_xdf('sub-P001_ses-S001_task-Default_run-001_20220313.xdf');
[streams,fileheader] = load_xdf('sub-P001_ses-S001_task-Default_run-001_eeg2022021503.xdf');
[streams,fileheader] = load_xdf('2022021507.xdf');
%%
info1=streams{1,1}.info;
info1.name % OxySoft
time_series1=streams{1,1}.time_series;
time_stamps1=streams{1,1}.time_stamps;
segments1=streams{1,1}.segments;
srate=info1.effective_srate;
clock_offsets=info1.clock_offsets;


plot(time_stamps1,time_series1(1,:));hold on;
%%
info2=streams{1,2}.info;
info2.name % OxySoft Event Marker
time_series2=streams{1,2}.time_series; % noting here
time_stamps2=streams{1,2}.time_stamps; % noting here

% plot(time_stamps2(1),2, 'bo');
% plot(time_stamps2(2),2.5, 'ro'); baseline_duration=num2str(time_stamps2(2)-time_stamps2(1));
% plot(time_stamps2(3:2:25),2, 'bo'); 
% plot(time_stamps2(4:2:26),2.5, 'ro');trials_durations=num2str(time_stamps2(4:2:26)-time_stamps2(3:2:25));

index_Baseline_Start = find(strcmp(horzcat({time_series2{:}}),'Baseline-Start')==1);
index_Baseline_End = find(strcmp(horzcat({time_series2{:}}),'Baseline-End')==1);
index_Trial_Start = find(strcmp(horzcat({time_series2{:}}),'Trial-Start')==1);
index_Trial_End = find(strcmp(horzcat({time_series2{:}}),'Trial-End')==1);
index_Screen_Flip = find(strcmp(horzcat({time_series2{:}}),'Screen_Flip')==1);


plot(time_stamps2(index_Baseline_Start),2, 'bo');
plot(time_stamps2(index_Baseline_End),2.5, 'ro'); 
baseline_duration=num2str(time_stamps2(index_Baseline_End(1))-time_stamps2(index_Baseline_Start(1)))
plot(time_stamps2(index_Trial_Start),2, 'bo'); 
plot(time_stamps2(index_Trial_End),2.5, 'ro');
trials_durations=num2str(time_stamps2(index_Trial_End(1:end))-time_stamps2(index_Trial_Start(1:end)))
plot(time_stamps2(index_Screen_Flip),3, 'go');

figure; plot(diff(time_stamps2(index_Screen_Flip)),'ro');
diff(time_stamps2(index_Screen_Flip))
%%
info3=streams{1,3}.info;
info3.name %BioSemi
time_series3=streams{1,3}.time_series;
time_stamps3=streams{1,3}.time_stamps;  num2str(time_stamps3(1))
segments3=streams{1,3}.segments;

plot(time_stamps3, time_series3, 'ro')

plot(time_stamps3(1),2, 'bo');
plot(time_stamps3(2),2.5, 'ro');
plot(time_stamps3(3:2:25),2, 'bo');
plot(time_stamps3(4:2:26),2.5, 'ro');
%%
info4=streams{1,4}.info;
info4.name % MyMarkerStream
time_series4=streams{1,4}.time_series;
time_stamps4=streams{1,4}.time_stamps; num2str(time_stamps4(1))

plot(time_stamps4, 2, 'bo');

figure; bar(time_series4);
figure; scatter(time_series4);

n=1;
num2str(time_stamps1(n))
num2str(time_stamps3(n))
num2str(time_stamps4(n))

num2str(time_stamps1(end))
num2str(time_stamps3(end))
num2str(time_stamps4(end))