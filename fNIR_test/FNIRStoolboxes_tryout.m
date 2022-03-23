%% load oxy5 with Oxysoft2matlab
clear
cd C:\Users\NIRS\Desktop\USB documents\Matlab\matlab file converter
addpath Matlab

cd C:\Users\NIRS\Documents\zhibin

clear;
oxysoft2matlab

% oxysoft2matlab ({'20211214.oxy5', 'Measurement Examples.oxyproj'})

DataFNIR=ans{1};

oxyvals=DataFNIR.oxyvals;
dxyvals=DataFNIR.dxyvals;
time=DataFNIR.time;
SampleNo=DataFNIR.sampleNo;
wavelengths=DataFNIR.wavelengths;
DPF=DataFNIR.DPF;
diatance=DataFNIR.distance;
Fs=DataFNIR.Fs;
label=DataFNIR.label;
RxLabel=DataFNIR.RxLabel;
TxLabel=DataFNIR.TxLabel;
Rx_TxId=DataFNIR.Rx_TxId;
chanPos=DataFNIR.chanPos;
transPos=DataFNIR.transPos;
receiPos=DataFNIR.receiPos;

ADlabel=DataFNIR.ADlabel;
ADvalues=DataFNIR.ADvalues; 
figure;plot(ADvalues(:,14));% PortAd_Buttons
figure;plot(ADvalues(:,15));% PortAd_Input

events=DataFNIR.events;


%% FieldTrip load oxy5
restoredefaultpath
addpath C:\Users\NIRS\Documents\fieldtrip\fieldtrip-20220304
ft_defaults

clear
cd C:\Users\NIRS\Documents\zhibin\20220215
cfg             = [];
cfg.dataset     = '20220215.oxy5';
data_raw        = ft_preprocessing(cfg);

trial=cell2mat(data_raw.trial);
time=cell2mat(data_raw.time);
label=data_raw.label;
hdrlabel=data_raw.hdr.label;
figure;plotx((trial(end-6:end,:))');
ADC001=trial(45,:);figure;plot(ADC001');
ADC007=trial(51,:);figure;plot(ADC007');
ADC008=trial(52,:);figure;plot(ADC008'); % look like steps
ADC009=trial(53,:);figure;plot(ADC009');
ADC010=trial(end-5,:);figure;plot(ADC010');
ADC011=trial(end-4,:);figure;plot(ADC011');
ADC012=trial(end-3,:);figure;plot(ADC012');hold on; % H bottom
ADC013=trial(end-2,:);plot(ADC013'); % H and L bottom
ADC014=trial(end-1,:);figure;plot(ADC014'); % noisy analog/TTL signal
ADC015=trial(end,:);figure;plot(ADC015'); % better analog/TTL signal

%% Load LSL data
cd C:\Users\NIRS\Documents\zhibin\20220215
addpath(genpath('C:\Users\zhibin\Documents\GitHub\xdf-Matlab')); 
fileList = dir('*.xdf');
fileList.name
[streams,fileheader] = load_xdf(fileList.name);
info1=streams{1,1}.info;
info1.name % OxySoft
time_series1=streams{1,1}.time_series;
time_stamps1=streams{1,1}.time_stamps;
segments1=streams{1,1}.segments;
srate=info1.effective_srate;
clock_offsets=info1.clock_offsets;

% Replace trial with LSL xdf data
data_raw.trial=time_series1;
data_raw.time=time_stamps1;
data_raw.hdr.nSamples=length(time_stamps1);
data_raw.sampleinfo=[1,length(time_stamps1)];
data_raw.label=label(1:size(time_series1,1));
data_raw.hdr.label=hdrlabel(1:size(time_series1,1));
data_raw.hdr.nChans=size(time_series1,1);
data_raw.hdr.chantype(size(time_series1,1)+1:end)=[];
data_raw.hdr.chanunit(size(time_series1,1)+1:end)=[];
data_raw.hdr.Fs=srate;
data_raw.fsample=srate;

%% field trip functions again
cfg = [];
cfg.ylim = 'maxabs';
ft_databrowser(cfg, data_raw);

cfg = [];
cfg.ylim = 'maxmin';
% cfg.channel = {'Rx4b-Tx5 [860nm]', 'Rx4b-Tx5 [764nm]'};  % you can also use wildcards like 'Rx4b-Tx5*'
cfg.channel = {'Rx1-Tx1*'}; 
ft_databrowser(cfg, data_raw);

% Exercise 1
cfg = [];
cfg.artfctdef.zvalue.channel = {'Rx1-Tx1 [760nm]', 'Rx1-Tx1 [842nm]'};
cfg.artfctdef.zvalue.cutoff = 5;
cfg.artfctdef.zvalue.hpfilter = 'yes';
cfg.artfctdef.zvalue.hpfreq = 0.1;
cfg.artfctdef.zvalue.rectify = 'yes';
cfg.artfctdef.zvalue.artpadding = 2;
% cfg.artfctdef.zvalue.interactive = 'yes'; % the interactive display makes more sense after segmentating data in trials
[cfg, artifact] = ft_artifact_zvalue(cfg, data_raw);

% Exercise 2
cfg = [];
cfg.dpf = 5.9;
cfg.channel = {'Rx1-Tx1 [760nm]', 'Rx1-Tx1 [842nm]'};
data_conc = ft_nirs_transform_ODs(cfg, data_raw);

% Exercise 3 (bandpass filter)
cfg = [];
cfg.ylim = 'maxmin';
% cfg.channel = {''Rx1-Tx1 [760nm]', 'Rx1-Tx1 [842nm]''};  % you can also use wildcards like 'Rx4b-Tx5*'
cfg.channel = {'Rx1-Tx1*'}; 
ft_databrowser(cfg, data_conc);

cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = [0.01 0.1];
data_filtered = ft_preprocessing(cfg, data_conc);

% Define epochs of interest
help ft_definetrial

cfg = [];
cfg.dataset = '2022031203.oxy3';
cfg.trialdef.eventtype = '?';

ft_definetrial(cfg);

cfg.trialdef.eventtype  = 'event';
cfg.trialdef.eventvalue = 'A';
cfg.trialdef.prestim    = 10;
cfg.trialdef.poststim   = 35;

cfg = ft_definetrial(cfg);

%% NIRS_SPM 
% refer to user manual in zotero 
% (Need Artinis_OXYMON_SampleData.nir in below directory)
% C:\Users\NIRS\Documents\NIRS_SPM_v4_r1\Sample_data\NIRS_data_file  
NIRS_SPM

%% NIRSTORM
% https://www.artinis.com/blogpost-all/2021/fnirs-analysis-toolbox-series-nirstorm
% https://neuroimage.usc.edu/brainstorm/Installation
cd C:\Users\NIRS\Documents\brainstorm3
brainstorm

cd C:\Users\NIRS\Documents\nirstorm-0.51
run nst_install.m
nst_install('copy');

% https://github.com/Nirstorm/nirstorm/wiki#tutorials

% This tutorial contain example data
% https://neuroimage.usc.edu/brainstorm/Tutorials/NIRSFingerTapping

% This tutorial contain example data
%https://github.com/Nirstorm/nirstorm/wiki/Workshop-PERFORM-Week-2018


%% NIRS-KIT 
% https://www.nitrc.org/frs/?group_id=1423&release_id=4540
% Zetero manual (Install directory contain the sample data)
% cd C:\Users\NIRS\Documents\NIRS_KIT_V2.0_202107\Sample_Data

%% Brain AnalyzIR (NIRS-TOOLBOX) (loading data not working)
clear
cd C:\Users\NIRS\Documents\GitHub\nirs-toolbox\+nirs\+io
file = 'C:\Users\NIRS\Documents\NIRS_SPM_v4_r1\Sample_data\NIRS_data_file\Artinis_OXYMON_SampleData.nir';
file = 'C:\Users\NIRS\Documents\zhibin\20220308LSL_analog\homer3_2022030801.snirf';
file = 'C:\Users\NIRS\Documents\zhibin\20211215left.snirf';
file = 'C:\Users\NIRS\Documents\zhibin\2022030806.nirs';
raw = nirs.io.loadDotNirs(file); % empty
raw.draw;
close;
raw.gui;
jobs = nirs.modules.TrimBaseline();
jobs.preBaseline = 5; % time before the first event  
jobs.postBaseline = 20; % time after the last event 
raw_exp = jobs.run(raw);
raw_exp = nirs.viz.StimUtil(raw_exp);

%% Homer 3 (lots of error, no idea)
% https://www.artinis.com/blogpost-all/2021/fnirs-analysis-toolbox-series-homer
% install Signal processing / Image Processing / Curve Fitting / Symbolic
% Math / Control System toolboxes and Simulink first
cd C:\Users\NIRS\Documents\Homer3
setpaths
Homer3
cd C:\Users\NIRS\Documents\AtlasViewer
setpaths
AtlasViewerGUI

clear
cd C:\Users\NIRS\Documents\zhibin\example_data
Homer3 % lots of errors


%% Nap-toolbox (No idea where to find it)
% https://www.artinis.com/blogpost-all/2021/fnirs-analysis-toolboxes-introduction
% https://openfnirs.org/software/other-fnirs-software/
