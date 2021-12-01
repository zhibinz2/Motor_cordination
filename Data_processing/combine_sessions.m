%% Stream line the steps

% Load behaviral data
cd /home/zhibin/Documents/GitHub/Motor_cordination/behaviraldata
% load matlab20211031.mat;
cd /home/zhibin/Documents/GitHub/Motor_cordination/EEGanalysis/20211102
run organize_Behavioral_step1.m

% Load EEG data and apply filter to get filtered_data
cd /home/zhibin/Documents/Acquisition/
% EEGfileName='bimanual_Reach_Jack_20211118'; EEG=loadcurry([pwd '/' EEGfileName '.cdt']);
cd /home/zhibin/Documents/GitHub/Motor_cordination/EEGanalysis/20211102
run organize_EEG_filter_step2.m
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

% Get goodchans
run Remove_EEG_artifact_Rereference_step6.m
% Just to examine
goodchans;

% Use Goodchans to do Laplacian
run Use_Goodchans_to_Laplacian_step7.m
% now filtered_data is after laplacian
filtered_data;%Just to examine

% Need to run step3 and step4 again to get new data_trials
cd /home/zhibin/Documents/GitHub/Motor_cordination/EEGanalysis/20211102
run integrate_EEG_into_data_trials_step3.m
run Integrate_Behavioral_into_data_trials_step4.m

%% combine data_trials from 3 sessions
% run section 1 again after each line of code

data_trials_day1=data_trials;
allPerm_day1=allPerm;
TrialScores_day1=TrialScores;

% run section 1 again with the data in day2
data_trials_day2=data_trials;
allPerm_day2=allPerm;
TrialScores_day2=TrialScores;

% run section 1 again with the data in day3
data_trials_day3=data_trials;
allPerm_day3=allPerm;
TrialScores_day3=TrialScores;


data_trials_alldays=cat(3,data_trials_day1,data_trials_day2,data_trials_day3);
allPerm_alldays=[allPerm_day1 allPerm_day2 allPerm_day3];
TrialScores_alldays=[TrialScores_day1 TrialScores_day2 TrialScores_day3];

%% Combine subjects
