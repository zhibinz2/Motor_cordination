%% load TMSi data
clear;close all;
cd /ssd/zhibin/1overf/20220324

[file, pathname] = uigetfile({'*.Poly5';'*.S00';'*.TMS32'},'Pick your file');
[path,filename,extension] = fileparts(file);5

fn=file;
addpath C:\Users\zhibi\Desktop\TMSi\SAGA Interface for Matlab (94-2303-0200-0-7)\SAGA_interface
cd  C:\Users\zhibi\Documents\TMSi\SAGA Interface for Matlab (94-2303-0200-0-7)\SAGA_interface
Poly5toEEGlab(fn)

cd C:\Users\zhibi\Desktop\TMSi
cd 'SAGA Interface for Matlab (94-2303-0200-0-7)'\
cd SAGA_interface
cd C:\Users\zhibi\Desktop\TMSi\SAGA Interface for Matlab (94-2303-0200-0-7)\SAGA_interface
d = TMSiSAGA.Poly5.read('C:\Users\HNL\Documents\TEST.DATA.Poly5');
d = TMSiSAGA.Poly5.read('C:\Users\zhibi\Documents\test2022031520220315T230153.DATA.Poly5');
d = TMSiSAGA.Poly5.read('C:\Users\HNL\Documents\testCedrus2022031620220316T011004.DATA.Poly5');
d = TMSiSAGA.Poly5.read('C:\Users\HNL\Documents\testCedrusRB2022032520220316T001249.DATA.Poly5');
d = TMSiSAGA.Poly5.read('C:\Users\HNL\Documents\20220316test20220316T012740.DATA.poly5');
d = TMSiSAGA.Poly5.read('C:\Users\zhibi\Documents\EMG_trigger20220320T223459.DATA.Poly5');
d = TMSiSAGA.Poly5.read('C:\Users\HNL\Documents\button_press2022-20220324T161705.DATA.Poly5');

samples=d.samples;
plot(samples(53,:),'ro'); % trigger
plot(samples(37,:),'ro'); hold on; % analog
plot(samples(38,:),'b'); % analog
plot(samples(39,:));
unique(samples(37,:))
C:\Users\HNL\Documents


plot(samples(53,:),'ro');

triggerSample=samples(53,:);

%% load LSL labrecorder file
addpath  C:\Users\zhibi\Documents\GitHub\xdf-Matlab
[streams,fileheader] = load_xdf('sub-P001_ses-S001_task-Default_run-001_eeg.xdf');


%%
cd C:\Users\HNL\Desktop\TMSi\SAGA Interface for Matlab (94-2303-0200-0-7)\SAGA_interface
d = TMSiSAGA.Poly5.read('C:\Users\HNL\Documents\TEST.DATA.Poly5');
d = TMSiSAGA.Poly5.read('C:\Users\HNL\Documents\testCedrusRB2022032520220316T001249.DATA.Poly5');
d = TMSiSAGA.Poly5.read('C:\Users\HNL\Documents\TMSimodel-20220316T013831.DATA.Poly5');
d = TMSiSAGA.Poly5.read('C:\Users\zhibi\Documents\photocell4000RB-20220322T215131.DATA.Poly5');

samples=d.samples;
plot(samples(53,:));
plot(samples(53,:),'ro');
unique(samples(53,:))
255-ans
plot(samples(53,:),'ro');
triggerSample=samples(53,:);

plot(samples(54,:));
unique(samples(54,:))


%-- 1/27/2022 4:11 PM --%
d = TMSiSAGA.Poly5.read('C:\Users\HNL\Documents\TEST.DATA.Poly5');
samples=d.samples;
plot(samples(53,:));
[file, pathname] = uigetfile({'*.Poly5';'*.S00';'*.TMS32'},'Pick your file');
clear
[file, pathname] = uigetfile({'*.Poly5';'*.S00';'*.TMS32'},'Pick your file');
[path,filename,extension] = fileparts(file);
[file, pathname] = uigetfile({'*.Poly5';'*.S00';'*.TMS32'},'Pick your file');
[path,filename,extension] = fileparts(file);
fn='TEST.DATA.Poly5';
Poly5toEEGlab(fn)
fn='TEST.DATA';
Poly5toEEGlab(fn)
fn=file;
Poly5toEEGlab(fn)
d = TMSiSAGA.Poly5.read('C:\Users\HNL\Documents\TEST.DATA.Poly5');
samples=d.samples;
d = TMSiSAGA.Poly5.read('C:\Users\HNL\Documents\TEST.DATA.Poly5');
clear
d = TMSiSAGA.Poly5.read('C:\Users\HNL\Documents\TEST.DATA.Poly5');
cd C:\Users\zhibi\Desktop\TMSi\SAGA Interface for Matlab (94-2303-0200-0-7)\SAGA_interface
d = TMSiSAGA.Poly5.read('C:\Users\HNL\Documents\TEST.DATA.Poly5');
cd C:\Users\zhibi\Desktop\TMSi\SAGA Interface for Matlab (94-2303-0200-0-7)\SAGA_interface
cd C:\Users\HNL\Documents
cd C:\Users\zhibi\Desktop\TMSi\SAGA Interface for Matlab (94-2303-0200-0-7)\SAGA_interface
cd C:\Users\zhibi\Desktop\TMSi
cd SAGA Interface for Matlab (94-2303-0200-0-7)
cd 'SAGA Interface for Matlab (94-2303-0200-0-7)'\
cd SAGA_interface
d = TMSiSAGA.Poly5.read('C:\Users\HNL\Documents\TEST.DATA.Poly5');
samples=d.samples;


d = TMSiSAGA.Poly5.read('C:\Users\zhibi\Documents\default-20220305T214217.DATA.Poly5');
samples=d.samples;
plot(samples(37,:));

%% load behaviral data
cd /home/zhibin/Documents/GitHub/1overf/stimulus_data_storage
load('20220324.mat')

%% auto correlation and spectrum