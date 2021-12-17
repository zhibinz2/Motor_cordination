% EEGLAB history file generated on the 17-Dec-2021
% ------------------------------------------------
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
EEG = pop_loadset('filename','eeglab_data_epochs_ica.set','filepath','/home/zhibin/Downloads/EEG CN/eeglab14_1_2b/sample_data/');
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
EEG = eeg_checkset( EEG );
figure; pop_plottopo(EEG, [1:32] , 'EEG Data epochs', 0, 'ydir',1);
pop_saveh( EEG.history, 'eeglabhist_plotmap_erp.m', '/home/zhibin/Documents/GitHub/Motor_cordination/EEGanalysis/EEGLABtryout/');
eeglab redraw;
