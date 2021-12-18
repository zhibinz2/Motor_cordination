% EEGLAB history file generated on the 18-Dec-2021
% ------------------------------------------------

EEG.etc.eeglabvers = '14.1.2'; % this tracks which version of EEGLAB is being used, you may ignore it
EEG = pop_importdata('dataformat','array','nbchan',128,'data','loadEEG','setname','loadEEG','srate',1000,'pnts',2000,'xmin',0.5,'chanlocs','/home/zhibin/Documents/GitHub/Motor_cordination/EEGanalysis/test3.ced');
EEG = eeg_checkset( EEG );
EEG=pop_chanedit(EEG, 'nosedir','+Y');
EEG = eeg_checkset( EEG );
figure; pop_plottopo(EEG, [1:128] , 'loadEEG', 0, 'ydir',1);
figure; pop_spectopo(EEG, 1, [500  2499], 'EEG' , 'percent', 15, 'freq', [6 10 22], 'freqrange',[2 25],'electrodes','off');
figure; pop_spectopo(EEG, 1, [500  2499], 'EEG' , 'percent', 15, 'freq', [6 10 22], 'freqrange',[2 25],'electrodes','off');
