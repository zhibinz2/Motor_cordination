% EEGLAB history file generated on the 17-Dec-2021
% ------------------------------------------------
EEG = eeg_checkset( EEG );
EEG = eeg_checkset( EEG );
figure;pop_topoplot(EEG,0, 1, 'EEG Data epochs',[1 1] ,0, 'electrodes', 'on', 'masksurf', 'on');
EEG = eeg_checkset( EEG );
EEG.chanlocs=pop_chanedit(EEG.chanlocs,  'plotrad',0.5);
EEG = pop_loadset( 'filename', 'eeglab_data_epochs_ica.set', 'filepath', '/data/common/matlab/eeglab/sample_data/');
EEG=pop_chanedit(EEG,  'lookup', '/data/common/matlab/eeglab/plugins/dipfit2.0/standard_BESA/standard-10-5-cap385.elp');
EEG = pop_saveset( EEG,  'savemode', 'resave');
EEG=pop_chanedit(EEG,  'eval', 'pop_writelocs( chans, ''/data/common/matlab/eeglab/sample_data/eeglab_chan32_2.locs'',  ''filetype'', ''loc'', ''format'',{ ''channum'', ''theta'', ''radius'', ''labels''}, ''header'', ''off'', ''customheader'','''');', 'load',{ '/data/common/matlab/eeglab/sample_data/eeglab_chan32_2.locs', 'filetype', 'autodetect'}, 'lookup', '/data/common/matlab/eeglab/plugins/dipfit2.0/standard_BESA/standard-10-5-cap385.elp');
EEG=pop_chanedit(EEG,  'settype',{ '1:32', 'EEG'}, 'changefield',{2, 'type', 'EOG'}, 'changefield',{6, 'type', 'EOG'});
EEG = pop_loadset( 'filename', 'eeglab_data_epochs_ica.set', 'filepath', '/data/common/matlab/eeglab/sample_data/');
EEG = pop_loadset('filename','eeglab_data_epochs_ica.set','filepath','/Users/arno/eeglab/sample_data/');
EEG = eeg_checkset( EEG );
EEG.etc.eeglabvers = '14.1.2'; % this tracks which version of EEGLAB is being used, you may ignore it
EEG = pop_loadset('filename','eeglab_data_epochs_ica.set','filepath','/home/zhibin/Downloads/EEG CN/eeglab14_1_2b/sample_data/');
EEG = eeg_checkset( EEG );
figure; pop_plottopo(EEG, [1:32] , 'EEG Data epochs', 0, 'ydir',1);
