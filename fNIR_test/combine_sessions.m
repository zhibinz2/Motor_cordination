% This script combine the data from 2022021503 -2022021508

%% initialize the first session
time_series1_trialsL_sessions=time_series1_trialsL;
time_series1_trialsR_sessions=time_series1_trialsR;

cd F:\UCI_dataset\fNIR\20220215\
save('20220215sessions.mat','time_series1_trialsL_sessions','time_series1_trialsR_sessions');

%% combine consecutive sessions
cd F:\UCI_dataset\fNIR\20220215\
load 20220215sessions.mat

time_series1_trialsL_sessions=cat(3,time_series1_trialsL_sessions,time_series1_trialsL);
time_series1_trialsR_sessions=cat(3,time_series1_trialsR_sessions,time_series1_trialsR);
save('20220215sessions.mat','time_series1_trialsL_sessions','time_series1_trialsR_sessions');

