EEG = pop_loadset('filename','eeglab_data_epochs_ica.set','filepath','/home/zhibin/Downloads/EEG CN/eeglab14_1_2b/sample_data/');
%%
laplacian_data;
loadEEG=laplacian_data';
%%
% organized as  channels X pnt X trials
EEG.data=permute(baselinecorrected_laplacian100_trial,[2 1 3]);
EEG.data=pow'
EEG.nbchan=128;
EEG.chanlocs=chanlocs128; 

X = [chanlocs128.X];
Y = [chanlocs128.Y];
Z = [chanlocs128.Z];

ZZ=(0.25-(Z+0.05))./1; %Make higher Z smaller fraction

XXPLOT=(-1)*(10*Y)*0.95.*ZZ*5;
YYPLOT=(10*X)*0.8.*ZZ*6;

plot(XXPLOT,YYPLOT,'ro');

labels = {'FP1','FPZ','FP2','AF3','AF4','F11','F7','F5','F3','F1','FZ','F2','F4','F6','F8','F12','FT11','FC5','FC3','FC1','FCZ','FC2','FC4','FC6','FT12','T7','C5','C3','C1','CZ','C2','C4','C6','T8','TP7','CP5','CP3','CP1','CPZ','CP2','CP4','CP6','TP8','M1','M2','P7','P5','P3','P1','PZ','P2','P4','P6','P8','PO7','PO3','POZ','PO4','PO8','O1','OZ','O2','CB1','CB2','AFP1','AFP2','AF7','AF5','AFZ','AF6','AF8','AFF5H','AFF3H','AFF1H','AFF2H','AFF4H','AFF6H','F9','F10','FFT7H','FFC5H','FFC3H','FFC1H','FFC2H','FFC4H','FFC6H','FFT8H','FT9','FT7','FT8','FT10','FTT7H','FCC5H','FCC3H','FCC1H','FCC2H','FCC4H','FCC6H','FTT8H','TTP7H','CCP5H','CCP3H','CCP1H','CCP2H','CCP4H','CCP6H','TTP8H','TPP7H','CPP5H','CPP3H','CPP1H','CPP2H','CPP4H','CPP6H','TPP8H','P9','P10','PPO3H','PPO1H','PPO2H','PPO4H','PO9','PO5','PO1','PO2','PO6','PO10','CBZ'};;
text(XXPLOT,YYPLOT,labels,'VerticalAlignment','bottom','HorizontalAlignment','right');

%% Plot all powerspetra on scalp map
XX=(10*X+1)*0.455;YY=(10*Y+1)*0.455;
for chan=1:128
    subplot('Position',[XXPLOT YYPLOT 0.04 0.02]);
    plot(freqs,pow(:,chan));ylim([0 0.5]);
    %xlabel('freq');ylabel('pow');
    title([AllchanNames{chan}]);
end
for chan=1:128
    subplot('Position',[XXPLOT YYPLOT 0.04 0.02]);
    ylim([0 0.5]);
end
suptitle('session2021111802')
%%

EEG.pnts=size(baselinecorrected_laplacian100_trial,1);
EEG.pnts=46;
EEG.xmin=-0.5;
EEG.xmax=1.5;

uilist    = { { 'style' 'text' 'string' 'Channels to plot' } ...
              { 'style' 'edit' 'string' [ '1:' num2str( EEG.nbchan ) ] 'tag' 'chan' } ...
              { 'style' 'text' 'string' 'Plot title' } ...
              { 'style' 'edit' 'string' fastif(isempty(EEG.setname), '',EEG.setname)  'tag' 'title' } ...                  
              { 'style' 'text' 'string' 'Plot single trials' } ...
              { 'style' 'checkbox' 'string' '(set=yes)' 'tag' 'cbst' } ...
              { 'style' 'text' 'string' 'Plot in rect. array' } ...
              { 'style' 'checkbox' 'string' '(set=yes)' 'tag' 'cbra' } ...
              { 'style' 'text' 'string' 'Other plot options (see help)' } ...
              { 'style' 'edit' 'string' '''ydir'', 1' 'tag' 'opt' } };
geometry = { [1 1] [1 1] [1 1] [1 1] [1 1] };

[result userdata tmphalt restag ] = inputgui( 'uilist', uilist, 'geometry', geometry, 'helpcom', 'pophelp(''pop_plottopo'')', 'title', 'Topographic ERP plot - pop_plottopo()');
 
channels     = eval( [ '[' restag.chan ']' ] );           
plottitle    = restag.title;
addoptions   = eval( [ '{' restag.opt '}' ] );

figure('name', ' plottopo()');
options ={ 'frames' EEG.pnts 'limits' [EEG.xmin EEG.xmax 0 0]*1000 ...
           'title' plottitle 'chans' channels addoptions{:} };
options = { options{:} 'chanlocs' EEG.chanlocs };
           
% plottopo( EEG.data, options{:} );

plottopo(mean(EEG.data,3), options{:} );

%%
EEG = pop_importdata('dataformat','array','nbchan',128,'data','loadEEG','setname','loadEEG','srate',1000,'pnts',2000,'xmin',0.5,'chanlocs','/home/zhibin/Documents/GitHub/Motor_cordination/EEGanalysis/test3.ced');
EEG = eeg_checkset( EEG );
EEG.data=permute(baselinecorrected_laplacian100_trial,[2 1 3]);
EEG=pop_chanedit(EEG, 'nosedir','+Y');
figure; pop_plottopo(EEG, [1:128] , 'loadEEG', 0, 'ydir',1);
figure; pop_plottopo(EEG, [1:128] , 'loadEEG', 0, 'ydir',1);

%%
plottopo(pow', options{:} );