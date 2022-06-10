function [SyncoSlowShowframes] = resamBP(BP01,ifi,sr,waitframes);
% This function resample the SAGA recording of 2000Hz button presses
% boolean timeseries into new showframes time series for psychtoolbox in
% the new Syncoption 2 Persion experiment.

% addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing
BPintervals=Calinterval(BP01');
AddshowframesCum=cumsum(BPintervals);

% calculate downsample rate
FlipsFreq = (1/ifi)/waitframes; % how many pyschotoolbox flip per second (feedback/photocell freq)
DownSampleRate = sr/FlipsFreq;

% Downsample
AddshowframesDown=round(AddshowframesCum/DownSampleRate);

% Aligned to the 30th tap
MeanTapInterval2Hz=0.5; % second % Mean stimulus interval for 2Hz pacing
NumFramesInterval2Hz=round(MeanTapInterval2Hz/(ifi*waitframes));  
shift=AddshowframesDown(1)-(NumFramesInterval2Hz*29+1);
AddshowframesDown=AddshowframesDown-shift;
SyncoShowframes=[1:NumFramesInterval2Hz:NumFramesInterval2Hz*30 AddshowframesDown(2:end)]; 
SyncoSlowShowframes=SyncoShowframes*1.5;% slow down by x1.5
end

