% This script combine the data from 2022021503 -2022021508

%% initialize the first session
time_series1_trialsL_sessions=time_series1_trialsL;
time_series1_trialsR_sessions=time_series1_trialsR;

% cd F:\UCI_dataset\fNIR\20220215\
cd C:\Users\zhibi\Desktop\Artinis_NIRS\zhibin\20220223
save('20220223sessions.mat','time_series1_trialsL_sessions','time_series1_trialsR_sessions');

%% combine consecutive sessions
% cd F:\UCI_dataset\fNIR\20220215\
cd C:\Users\zhibi\Desktop\Artinis_NIRS\zhibin\20220223
load 20220223sessions.mat

time_series1_trialsL_sessions=cat(3,time_series1_trialsL_sessions,time_series1_trialsL);
time_series1_trialsR_sessions=cat(3,time_series1_trialsR_sessions,time_series1_trialsR);
save('20220223sessions.mat','time_series1_trialsL_sessions','time_series1_trialsR_sessions');

