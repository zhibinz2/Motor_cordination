addpath('/home/zhibin/Documents/GitHub/labstreaminglayer/LSL/liblsl-Matlab')


addpath('/home/zhibin/Documents/GitHub/xdf-Matlab')

%%
streamsFNIR = load_xdf('sub-P001_ses-S001_task-Default_run-001_eeg.xdf')

stream1=streamsFNIR{1};% markers sent
stream2=streamsFNIR{2};% oxysoft markers
stream3=streamsFNIR{3};% oxysoft fnir time series 36*788

%%
streamsABM = load_xdf('sub-P001_ses-S001_task-Default_run-001_eeg.xdf')

stream1=streamsABM{1};% markers sent
stream2=streamsABM{2};% oxysoft markers