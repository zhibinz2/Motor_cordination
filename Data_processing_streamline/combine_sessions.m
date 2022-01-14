%% Stream line the steps
clear

% Load behaviral data
cd /home/zhibin/Documents/GitHub/Motor_cordination/behaviraldata
open All_subjects.m
% load matlab20211031.mat;
% load 2021111802.mat
cd /home/zhibin/Documents/GitHub/Motor_cordination/EEGanalysis/20211102
run organize_Behavioral_step1.m

% Load EEG data and apply filter to get filtered_data
cd /home/zhibin/Documents/Acquisition/
% cd /home/zhibin/Documents/Acquisition/bimanual_Reach_Jack_2021111802
% EEGfileName='bimanual_Reach_Hiro_20211119'; EEG=loadcurry([pwd '/' EEGfileName '.cdt']);
fileList = dir('*.cdt');EEG=loadcurry([pwd '/' fileList.name]);
cd /home/zhibin/Documents/GitHub/Motor_cordination/EEGanalysis/20211102
run organize_photocells_step2.m
run EEG_filter_step2b.m
% Just to examine
filtered_data;

% To get data_trials;
run integrate_EEG_into_data_trials_step3.m
% Just to examine
data_trials;
% "500ms padding + 500ms green dot baseline + 500ms plan phase + 1000ms movement + 1000ms show bonus + 500ms padding"

% step4 resample the movements
run Integrate_Behavioral_into_data_trials_step4.m
% Plot movements (just to check)
run draw_movements_all_conditions_step5.m
% Speed profile
cd /home/zhibin/Documents/GitHub/Motor_cordination/PerformanceAnaylisis
open speed_derivative.m

% Get goodchans, and a boolean goodepochs array
% run Remove_EEG_artifact_Rereference_step6.m
run artifact_removal_step6b.m
% Just to examine
reRef_data;
goodchans;
Bgoodepochs;% in boolean
sum(Bgoodepochs) % number of bad epochs

% Run ICA on reRef_data
% this should be step 6b
open RunICA_step6b.m
% input reRef_data and get mixedsig
mixedsig;

% Just to examine the final clean data
figure('units','normalized','outerposition',[0 0 1 0.5]);
plot(mixedsig(goodchans,:)');


%% Use Goodchans to do Laplacian on mixedsig (don't do laplacian)
run Use_Goodchans_to_Laplacian_step7.m
% now filtered_data is after laplacian
laplacian_data;%Just to examine

%% Cut the data and reorganzied into data_trials format (skip this when do without laplacian)
% try using EEGLAB
% eeglab

% try using my own method
% Initialize the data_trials matrix
laplacian_trials=zeros(500+NumTrialtimepoints,NumEEGChannels,length(goodepochs));% add the 500 ms before trial as baseline
for ntr=1:length(goodepochs) 
    laplacian_trials((1:(500+NumTrialtimepoints)),1:NumEEGChannels,ntr)=laplacian_data((IndEvents(ntr)):(IndEnds(ntr)),1:NumEEGChannels); 
end
% Just to examine
laplacian_trials;

%% Cut the data and reorganzied into data_trials format 
% try using EEGLAB
% eeglab

% try using my own method
mixedsig=mixedsig';
% Initialize the data_trials matrix
afterICA_trials=zeros(NumTrialtimepoints,NumEEGChannels,length(goodepochs));% add the 500 ms before trial as baseline
for ntr=1:length(goodepochs) % ntr=length(goodepochs)
    afterICA_trials(1:NumTrialtimepoints,1:NumEEGChannels,ntr)=mixedsig((IndEnds(ntr)-1999):(IndEnds(ntr)),1:NumEEGChannels); 
end
% Just to examine
afterICA_trials;
plotx(afterICA_trials(:,goodchans,2));

%% Need to step4 again to integrate behaviral data?
open Integrate_Behavioral_into_data_trials_step4.m

%% Step 8 
open Powerspetra.m

%% Step 9
open plotSpectrogram.m

%% Save data file for sharing
cd /home/zhibin/Documents/Acquisition/bimanual_Reach_Jack_2021111802/Cleaner_Data
run Cleaner_Data_explained.m




















%% combine data_trials from 3 sessions
% run section 1 again after each line of code

data_trials_day1=data_trials;
allPerm_day1=allPerm;
TrialScores_day1=TrialScores;
goodepochs_day1=Bgoodepochs;

% run section 1 again with the data in day2
data_trials_day2=data_trials;
allPerm_day2=allPerm;
TrialScores_day2=TrialScores;
goodepochs_day2=Bgoodepochs;

% run section 1 again with the data in day3
data_trials_day3=data_trials;
allPerm_day3=allPerm;
TrialScores_day3=TrialScores;
goodepochs_day3=Bgoodepochs;


data_trials_alldays=cat(3,data_trials_day1,data_trials_day2,data_trials_day3);
allPerm_alldays=[allPerm_day1 allPerm_day2 allPerm_day3];
TrialScores_alldays=[TrialScores_day1 TrialScores_day2 TrialScores_day3];
goodepochs_alldays=[goodepochs_day1 goodepochs_day2 goodepochs_day3];

%% *100 and Baseline correct
% (So the laplacian is in microvolts/mm^2 right now.    more typical in published literature would be microvolts/cm^2.   in which case you should multiply by 100 the raw Laplacian values before FFT)
data_trials_alldays=data_trials_alldays.*100;
baselinesamps = 501:1000; % use the first 500ms as baseline
baselinecorrected_data_trial=zeros(size(data_trials_alldays)); 
for i=1:size(data_trials_alldays,3) % loop through trials
    erpdata=data_trials_alldays(:,:,i);
    newerp = baselinecorrect(erpdata,baselinesamps);
    baselinecorrected_data_trial(:,:,i)=newerp;
end
% plotx(mean(baselinecorrected_data_trial(:,1:128,:),3));
% xticks(linspace(0,4000,9));xticklabels({'-1000','-500','0','500','1000','1500','2000','2500','3000'});
baselinecorrected_data_trial;

%% Combine subjects




