% cd D:\360MoveData\Users\alienware\Documents\GitHub\Motor_cordination

% conditions
allPerm;
trialsL=find(allPerm==1);
trialsR=find(allPerm==2);

% channels
O2HbChans=1:2:43;
HHbChans=2:2:44;

O2HbChansL=1:2:21;
HHbChansL=2:2:22;

O2HbChansR=23:2:43;
HHbChansR=24:2:44;

time_series1(O2HbChansL);
time_series1(O2HbChansR);

time_series1(HHbChansL);
time_series1(HHbChansR);

% organized time_series1 into trials
time_series1_baseline=time_series1(stamp1index_Baseline_Start:stamp1index_Baseline_End);

trial_length=max(stamp1index_Trial_End-stamp1index_Trial_Start);
NumTrials=length(stamp1index_Trial_Start);
NumChans=size(time_series1,1);

time_series1_trials=zeros(4*trial_length,NumChans,NumTrials);
for i=1:NumTrials
    time_series1_trials(:,1:NumChans,i)=time_series1(:,stamp1index_Trial_Start(i)-trial_length:stamp1index_Trial_Start(i)+3*trial_length-1)';
end
time_series1_trials;

% baseline correction
baseline_averagesL=mean(mean(time_series1_trials(1:trial_length,:,trialsL),3),1);
baseline_averagesR=mean(mean(time_series1_trials(1:trial_length,:,trialsR),3),1);

time_series1_trialsL=time_series1_trials(:,:,trialsL)-baseline_averagesL.*ones(size(time_series1_trials,1),1).*ones(1,1,length(trialsL));
time_series1_trialsR=time_series1_trials(:,:,trialsR)-baseline_averagesR.*ones(size(time_series1_trials,1),1).*ones(1,1,length(trialsR));

