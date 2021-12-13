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
tic
run organize_EEG_filter_step2.m
toc
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

% Get goodchans, and a boolean goodepochs array
run Remove_EEG_artifact_Rereference_step6.m
% Just to examine
reRef_data
goodchans;
Bgoodepochs;% in boolean

% Run ICA on reRef_data
% this should be step 6b
reRef_data
% to be continue...


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




