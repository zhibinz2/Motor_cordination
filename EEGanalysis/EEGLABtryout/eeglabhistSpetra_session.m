% EEGLAB history file generated on the 18-Dec-2021
% ------------------------------------------------
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
EEG = pop_importdata('dataformat','array','nbchan',128,'data','loadEEG','setname','loadEEG','srate',1000,'pnts',2000,'xmin',0.5,'chanlocs','/home/zhibin/Documents/GitHub/Motor_cordination/EEGanalysis/test3.ced');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'gui','off'); 
EEG=pop_chanedit(EEG, 'nosedir','+Y');
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset( EEG );
figure; pop_plottopo(EEG, [1:128] , 'loadEEG', 0, 'ydir',1);
pop_saveh( EEG.history, 'eeglab128hist.m', '/home/zhibin/Documents/GitHub/Motor_cordination/EEGanalysis/EEGLABtryout/');
EEG = eeg_checkset( EEG );
figure; pop_spectopo(EEG, 1, [500  2499], 'EEG' , 'percent', 15, 'freq', [6 10 22], 'freqrange',[2 25],'electrodes','off');
EEG = eeg_checkset( EEG );
figure; pop_spectopo(EEG, 1, [500  2499], 'EEG' , 'percent', 15, 'freq', [6 10 22], 'freqrange',[2 25],'electrodes','off');
pop_saveh( EEG.history, 'eeglabhistPLOTspetra.m', '/home/zhibin/Documents/GitHub/Motor_cordination/EEGanalysis/EEGLABtryout/');
eeglab redraw;
