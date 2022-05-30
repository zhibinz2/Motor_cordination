% check number of BP in previous recording
cd /ssd/zhibin/1overf/20220515_2P/Segmented_data
[BPCondi_dataName]=uigetfile('BP*.mat');% select the BP file
Path_filename=[pwd '/' BPCondi_dataName];
load(Path_filename);
[BPCondi_dataName]=uigetfile('FB*.mat');% select the FB file
Path_filename=[pwd '/' BPCondi_dataName];
load(Path_filename);

cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/exp_uni_2p/Syncopation_2P
[ShowframesL] = resamBP(BPCondi1L,ifi,sr,waitframes);
[ShowframesR] = resamBP(BPCondi1R,ifi,sr,waitframes);

[ShowframesL] = resamBP(FBCondi1R,ifi,sr,waitframes);
[ShowframesR] = resamBP(BPCondi1L,ifi,sr,waitframes);
