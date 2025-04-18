function [SyncoSlowShowframes] = resamBP(BP01,ifi,sr,waitframes);
% This function resample the SAGA recording of 2000Hz button presses
% boolean timeseries into new showframes time series for psychtoolbox in
% the new Syncoption 2 Persion experiment.

% addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing
BPintervals=Calinterval(BP01');
AddshowframesCum=cumsum([BPintervals]);

% calculate downsample rate
FlipsFreq = (1/ifi)/waitframes; % how many pyschotoolbox flip per second (feedback/photocell freq)
DownSampleRate = sr/FlipsFreq;

% Downsample
AddshowframesDown=round(AddshowframesCum/DownSampleRate);

% Add pacers and ligned to the 30th tap
MeanTapInterval13Hz=1/1.3; % second % Mean stimulus interval for 2Hz pacing
NumFramesInterval13Hz=round(MeanTapInterval13Hz/(ifi*waitframes));  
shift=AddshowframesDown(1)-(NumFramesInterval13Hz*29+1);
AddshowframesDown=AddshowframesDown-shift;
SyncoShowframes=[1:NumFramesInterval13Hz:NumFramesInterval13Hz*30 AddshowframesDown(2:end)]; 
SyncoSlowShowframes=SyncoShowframes*1;% slow down by x1
end

