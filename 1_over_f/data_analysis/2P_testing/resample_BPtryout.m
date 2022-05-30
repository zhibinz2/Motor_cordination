%% resample key presses
cd /ssd/zhibin/1overf/20220515_2P/Segmented_data
cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/exp_uni_2p/Syncopation_2P

BPCondi1L;BPCondi1R;
% select the time series
BP01=BPCondi1L;
% convert to intervals
cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing
BPintervals=Calinterval(BP01');
AddshowframesCum=cumsum(BPintervals);

figure;
% subplot(3,1,1);
plot(BPintervals,'b.');hold on;
plot(AddshowframesCum,'mo');hold off;
xlim([0 60]);title('BP intervals & cumsum');
legend('BP intervals','Cumsum');

% resample
sr=unique([srL srR]) % TMSi SAGA recording sampling rate in Hz

% srShowframe=sr*ifi*waitframes % samples per screen flip rate interval
% sr/srShowframe % how many Psychtoolbox flips per second
% DownSampleRate=sr/(sr/srShowframe)

open /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing/timeing_2P_check.m
open /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/exp_uni_2p/5_conditions/TestFliping_time.m

FlipsFreq = (1/ifi)/waitframes % how many pyschotoolbox flip per second (feedback/photocell freq)
DownSampleRate = sr/FlipsFreq 

% Downsample
AddshowframesDown=round(AddshowframesCum/DownSampleRate);
% Addshowframes=round(Addshowframes/(sr/(1/ifi))/waitframes); % Downsampling rate at 40

% subplot(3,1,2);
figure;
plot(AddshowframesCum,'mo'); % before downsampling 
hold on;
plot(AddshowframesDown,'k+'); % downsampled frames for psychtoolbox
hold off;
xlim([0 60]);title('before downsampling & after downsampling');
legend('before downsampling', 'after downsampling');

% Aligned to the 30th tap
AddshowframesDown(1) 
NumFramesInterval2Hz*29+1 % num of farmes at the 30th tap
shift=AddshowframesDown(1)-(NumFramesInterval2Hz*29+1)
AddshowframesDown=AddshowframesDown-shift;
AddshowframesDown(1) == NumFramesInterval2Hz*29+1 % test if they were aligned
SyncoShowframes=[]; % create Showframes for Syncopated expt
SyncoShowframes=[1:NumFramesInterval2Hz:NumFramesInterval2Hz*30 AddshowframesDown(2:end)]; 


cd /ssd/zhibin/1overf/20220515_2P
load('20220515.mat')
figure;
plot(Showframes1,'go');hold on;
plot(SyncoShowframes,'k*');
xlim([0 60]);
title('syn showframes1 & synco showframe');
legend('syn showframes1', 'synco showframe');

figure;
plot(SyncoShowframes,'go');hold on;
SyncoSlowShowframes=SyncoShowframes*1.5;
plot(SyncoSlowShowframes,'k*');
xlim([0 60]);
title('before slowing & after slowing');
legend('before slowing', 'after slowing');

%% call the above in a function resamBP
sr=unique([srL srR]);sr=2000;
cd /ssd/zhibin/1overf/20220515_2P
load('20220515.mat')
ifi;waitframes;
cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/exp_uni_2p/Syncopation_2P
[SyncoSlowShowframes] = resamBP(BP01,ifi,sr,waitframes);


%% Randomization of the pacer stimuli in each condition ****************************************
%     % Mean stimulus interval for 2Hz pacing
%     MeanTapInterval2Hz=0.5; % second
%     NumFramesInterval2Hz=round(MeanTapInterval2Hz/(ifi*waitframes));  
%     
%     % condition 1-4 (paced the frist 30 taps, set the rest of the frames with value zeros)
%     Showframes1=[1:NumFramesInterval2Hz:NumFramesInterval2Hz*30 zeros(1,570)]; 
%     Showframes2=Showframes1;
%     Showframes3=Showframes1;
%     Showframes4=Showframes1;
%     
%     % Mean stimulus interval for 3Hz pacing
%     MeanTapInterval3Hz=1/3; % second
%     NumFramesInterval3Hz=round(MeanTapInterval3Hz/(ifi*waitframes));  % on average 72 frames per stimulus 
%     % condition 5
%     Showframes5=[1:NumFramesInterval3Hz:NumFramesInterval3Hz*30 zeros(1,570)]; 
%     Showframes6=Showframes5;
%     
%     % combine all 6 conditions
%     Showframes=[Showframes1;Showframes2;Showframes3;Showframes4;Showframes5;Showframes6];