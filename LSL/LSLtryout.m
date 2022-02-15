%% Sending Data
% instantiate the library
disp('Loading library...');

% addpath C:\Users\zhibi\Documents\GitHub\labstreaminglayer\LSL\liblsl-Matlab\bin\
% addpath /home/hnl/Documents/GitHub/labstreaminglayer/LSL/liblsl-Matlab
addpath(genpath('/home/hnl/Documents/GitHub/labstreaminglayer/LSL/liblsl-Matlab'));
cd /home/hnl/Documents/GitHub/labstreaminglayer/LSL/liblsl-Matlab
lib = lsl_loadlib();

% make a new stream outlet
disp('Creating a new streaminfo...');
info = lsl_streaminfo(lib,'BioSemi','EEG',8,100,'cf_float32','sdfwerr32432');

disp('Opening an outlet...');
outlet = lsl_outlet(info);

while true
    % send data into the outlet, sample by sample
    disp('Now transmitting data...');
    outlet.push_sample(1);
end


%% Sending Marker
% instantiate the library
disp('Loading library...');

addpath C:\Users\zhibi\Documents\GitHub\labstreaminglayer\LSL\liblsl-Matlab\bin\
lib = lsl_loadlib();

% make a new stream outlet
% the name (here MyMarkerStream) is visible to the experimenter and should be chosen so that 
% it is clearly recognizable as your MATLAB software's marker stream
% The content-type should be Markers by convention, and the next three arguments indicate the 
% data format (1 channel, irregular rate, string-formatted).
% The so-called source id is an optional string that allows for uniquely identifying your 
% marker stream across re-starts (or crashes) of your script (i.e., after a crash of your script 
% other programs could continue to record from the stream with only a minor interruption).
disp('Creating a new marker stream info...');
info = lsl_streaminfo(lib,'MyMarkerStream','Markers',1,0,'cf_string','myuniquesourceid23443');

disp('Opening an outlet...');
outlet = lsl_outlet(info);

% send markers into the outlet
disp('Now transmitting data...');
markers = {'Test-Marker', 'Trial-Start', 'Trial-End'};
mrk = markers{1};

while
    disp(['now sending ' mrk]);
    outlet.push_sample({mrk});   % note that the string is wrapped into a cell-array
end


%% load LSL labrecorder file

% addpath  C:\Users\zhibi\Documents\GitHub\xdf-Matlab
addpath /home/zhibin/Documents/GitHub/xdf-Matlab
cd /home/zhibin/Documents/Artinis NIRS/zhibin/20220213fNIRS_LSL

[streams,fileheader] = load_xdf('sub-P001_ses-S001_task-Default_run-001_20220313.xdf');

info1=streams{1,1}.info;
info1.name % OxySoft
time_series1=streams{1,1}.time_series;
time_stamps1=streams{1,1}.time_stamps;
segments1=streams{1,1}.segments;

plot(time_stamps1,time_series1);

info2=streams{1,2}.info;
info2.name % OxySoft Event Marker
time_series2=streams{1,2}.time_series; % noting here
time_stamps2=streams{1,2}.time_stamps; % noting here


info3=streams{1,3}.info;
info3.name %BioSemi
time_series3=streams{1,3}.time_series;
time_stamps3=streams{1,3}.time_stamps;
segments3=streams{1,3}.segments;

info4=streams{1,4}.info;
info4.name % MyMarkerStream
time_series4=streams{1,4}.time_series;
time_stamps4=streams{1,4}.time_stamps;

figure; bar(time_series4);
figure; scatter(time_series4);