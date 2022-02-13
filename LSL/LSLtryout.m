%% Sending Data
% instantiate the library
disp('Loading library...');

addpath C:\Users\zhibi\Documents\GitHub\labstreaminglayer\LSL\liblsl-Matlab\bin\
lib = lsl_loadlib();

% make a new stream outlet
disp('Creating a new streaminfo...');
info = lsl_streaminfo(lib,'BioSemi','EEG',8,100,'cf_float32','sdfwerr32432');

disp('Opening an outlet...');
outlet = lsl_outlet(info);

% send data into the outlet, sample by sample
disp('Now transmitting data...');


outlet.push_sample(1);


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
disp(['now sending ' mrk]);
outlet.push_sample({mrk});   % note that the string is wrapped into a cell-array
