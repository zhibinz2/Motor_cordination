addpath('/home/zhibin/Documents/GitHub/labstreaminglayer/LSL/liblsl-Matlab')


addpath('/home/zhibin/Documents/GitHub/xdf-Matlab')


streams = load_xdf('sub-P001_ses-S001_task-Default_run-001_eeg.xdf')

stream1=streams{1};% markers sent
stream2=streams{2};% oxysoft markers
stream3=streams{3};% oxysoft fnir time series 36*788