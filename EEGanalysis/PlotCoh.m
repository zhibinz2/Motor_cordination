%% plot the coherence

UniCondi=unique(CondiData);

% First, plot the scalp coherence
rate=Fs;
maxfreq=50;
win=500:1500; % After 500 ms

%% All channels

for u=1:length(UniCondi)% u=2
    
    indtemp=find(CondiData==UniCondi(u));
    
    % pick a condition from randomized set allPerm
    conditionSelected = UniCondi(u);
    
    subplot(1,length(UniCondi),u);
    
    % First, plot the scalp coherence
    [pow,freqs,df,eppow,corr,cprod,fcoef] = allspectra(reRef_data,rate,maxfreq,indtemp,win);
    % calculate the coherence
    coh=abs(corr).^2;
    
    % Plot coh of specific frequency
    imagesc(squeeze(coh(46,:,:)));
    xlabel('channels');ylabel('channels');title({'45 Hz', 'coh - condition: ' num2str(conditions(u))});colorbar;
    
    % label all channels
    xticks([1:128]);
    xticklabels({'FP1','FPZ','FP2','AF3','AF4','F11','F7','F5','F3','F1','FZ','F2','F4','F6','F8','F12','FT11','FC5','FC3','FC1','FCZ','FC2','FC4','FC6','FT12','T7','C5','C3','C1','CZ','C2','C4','C6','T8','TP7','CP5','CP3','CP1','CPZ','CP2','CP4','CP6','TP8','M1','M2','P7','P5','P3','P1','PZ','P2','P4','P6','P8','PO7','PO3','POZ','PO4','PO8','O1','OZ','O2','CB1','CB2','AFP1','AFP2','AF7','AF5','AFZ','AF6','AF8','AFF5H','AFF3H','AFF1H','AFF2H','AFF4H','AFF6H','F9','F10','FFT7H','FFC5H','FFC3H','FFC1H','FFC2H','FFC4H','FFC6H','FFT8H','FT9','FT7','FT8','FT10','FTT7H','FCC5H','FCC3H','FCC1H','FCC2H','FCC4H','FCC6H','FTT8H','TTP7H','CCP5H','CCP3H','CCP1H','CCP2H','CCP4H','CCP6H','TTP8H','TPP7H','CPP5H','CPP3H','CPP1H','CPP2H','CPP4H','CPP6H','TPP8H','P9','P10','PPO3H','PPO1H','PPO2H','PPO4H','PO9','PO5','PO1','PO2','PO6','PO10','CBZ','VEOG','HEOG','EMG1','EMG2','HL 1','HL 2','EMG3','EMG4','EMG5','EMG6','TRIGGER'}); % all 139 channels (including the additional chans)
    xtickangle(90);
    yticks([1:128]);
    yticklabels({'FP1','FPZ','FP2','AF3','AF4','F11','F7','F5','F3','F1','FZ','F2','F4','F6','F8','F12','FT11','FC5','FC3','FC1','FCZ','FC2','FC4','FC6','FT12','T7','C5','C3','C1','CZ','C2','C4','C6','T8','TP7','CP5','CP3','CP1','CPZ','CP2','CP4','CP6','TP8','M1','M2','P7','P5','P3','P1','PZ','P2','P4','P6','P8','PO7','PO3','POZ','PO4','PO8','O1','OZ','O2','CB1','CB2','AFP1','AFP2','AF7','AF5','AFZ','AF6','AF8','AFF5H','AFF3H','AFF1H','AFF2H','AFF4H','AFF6H','F9','F10','FFT7H','FFC5H','FFC3H','FFC1H','FFC2H','FFC4H','FFC6H','FFT8H','FT9','FT7','FT8','FT10','FTT7H','FCC5H','FCC3H','FCC1H','FCC2H','FCC4H','FCC6H','FTT8H','TTP7H','CCP5H','CCP3H','CCP1H','CCP2H','CCP4H','CCP6H','TTP8H','TPP7H','CPP5H','CPP3H','CPP1H','CPP2H','CPP4H','CPP6H','TPP8H','P9','P10','PPO3H','PPO1H','PPO2H','PPO4H','PO9','PO5','PO1','PO2','PO6','PO10','CBZ','VEOG','HEOG','EMG1','EMG2','HL 1','HL 2','EMG3','EMG4','EMG5','EMG6','TRIGGER'}) % all 139 channels (including the additional chans)

end

%% Divide left and right hemisphere

% Left hemisphere
Leftticks=[1 4 6 7 8 9 10 17 18 19 20 26 27 28 29 35 36 37 38 44 46 47 48 49 55 56 60 63 65 67 68 72 73 74 78 80 81 82 83 88 89 92 93 94 95 100 101 102 103 108 109 110 111 116 118 119 122 123 124];
% Right hemisphere
Rightticks=[3 5 12 13 14 15 16 22 23 24 25 31 32 33 34 40 41 42 43 45 51 52 53 54 58 59 62 64 66 70 71 75 76 77 79 84 85 86 87 90 91 96 97 98 99 104 105 106 107 112 113 114 115 117 120 121 125 126 127];


for u=1:length(UniCondi)% u=2
    
    indtemp=find(CondiData==UniCondi(u));
    
    % pick a condition from randomized set allPerm
    conditionSelected = UniCondi(u);
    
    subplot(1,length(UniCondi),u);
    
    % First, plot the scalp coherence
    [pow,freqs,df,eppow,corr,cprod,fcoef] = allspectra(reRef_data,rate,maxfreq,indtemp,win);
    % calculate the coherence
    coh=abs(corr).^2;
    

    % Plot coh of specific frequency between the two hemispheres
    imagesc(squeeze(coh(46,Leftticks,Rightticks)));
    xlabel('right hemisphere channels');ylabel('left hemisphere channels');
    title({'45 Hz', 'coh - condition: ' num2str(conditions(u))});colorbar;

    % Left hemisphere
    yticks([1:59]);
    % Left hemisphere
    yticklabels({'FP1','AF3','F11','F7','F5','F3','F1','FT11','FC5','FC3','FC1','T7','C5','C3','C1','TP7','CP5','CP3','CP1','M1','P7','P5','P3','P1','PO7','PO3','O1','CB1','AFP1','AF7','AF5','AFF5H','AFF3H','AFF1H','F9','FFT7H','FFC5H','FFC3H','FFC1H','FT9','FT7','FTT7H','FCC5H','FCC3H','FCC1H','TTP7H','CCP5H','CCP3H','CCP1H','TPP7H','CPP5H','CPP3H','CPP1H','P9','PPO3H','PPO1H','PO9','PO5','PO1'});
    % Right hemisphere
    xticks([1:59]);
    % Right hemisphere
    xticklabels({'FP2','AF4','F2','F4','F6','F8','F12','FC2','FC4','FC6','FT12','C2','C4','C6','T8','CP2','CP4','CP6','TP8','M2','P2','P4','P6','P8','PO4','PO8','O2','CB2','AFP2','AF6','AF8','AFF2H','AFF4H','AFF6H','F10','FFC2H','FFC4H','FFC6H','FFT8H','FT8','FT10','FCC2H','FCC4H','FCC6H','FTT8H','CCP2H','CCP4H','CCP6H','TTP8H','CPP2H','CPP4H','CPP6H','TPP8H','P10','PPO2H','PPO4H','PO2','PO6','PO10'});
    xtickangle(90);
    
end


%% plot just one condtion
u=6; % plot just one condtion
    
indtemp=find(CondiData==UniCondi(u));

% find the high peformance trial index within one condition
HighInd=indtemp(find(TrialScores(indtemp)>0.6));
LowInd=indtemp(find(TrialScores(indtemp)<0.6));

% pick a condition from randomized set allPerm
conditionSelected = UniCondi(u);

% First, plot the scalp coherence
[pow,freqs,df,eppow,corr,cprod,fcoef] = allspectra(reRef_data,rate,maxfreq,HighInd,win);
% calculate the coherence
coh=abs(corr).^2;


% Plot coh of specific frequency between the two hemispheres
imagesc(squeeze(coh(46,Leftticks,Rightticks)));
xlabel('right hemisphere channels');ylabel('left hemisphere channels');
title({'45 Hz', 'coh - condition: ' num2str(conditions(u))});colorbar;

% Left hemisphere
yticks([1:59]);
% Left hemisphere
yticklabels({'FP1','AF3','F11','F7','F5','F3','F1','FT11','FC5','FC3','FC1','T7','C5','C3','C1','TP7','CP5','CP3','CP1','M1','P7','P5','P3','P1','PO7','PO3','O1','CB1','AFP1','AF7','AF5','AFF5H','AFF3H','AFF1H','F9','FFT7H','FFC5H','FFC3H','FFC1H','FT9','FT7','FTT7H','FCC5H','FCC3H','FCC1H','TTP7H','CCP5H','CCP3H','CCP1H','TPP7H','CPP5H','CPP3H','CPP1H','P9','PPO3H','PPO1H','PO9','PO5','PO1'});
% Right hemisphere
xticks([1:59]);
% Right hemisphere
xticklabels({'FP2','AF4','F2','F4','F6','F8','F12','FC2','FC4','FC6','FT12','C2','C4','C6','T8','CP2','CP4','CP6','TP8','M2','P2','P4','P6','P8','PO4','PO8','O2','CB2','AFP2','AF6','AF8','AFF2H','AFF4H','AFF6H','F10','FFC2H','FFC4H','FFC6H','FFT8H','FT8','FT10','FCC2H','FCC4H','FCC6H','FTT8H','CCP2H','CCP4H','CCP6H','TTP8H','CPP2H','CPP4H','CPP6H','TPP8H','P10','PPO2H','PPO4H','PO2','PO6','PO10'});
xtickangle(90);


%% plot just one condtion and only the prefrontal chans

% Left prefrontal
Leftticks=[1 4 6 7 8 9 10 72 73 74];
% Right prefrontal
Rightticks=[3 5 12 13 14 15 16 75 76 77];


u=4; % plot just one condtion
    
indtemp=find(CondiData==UniCondi(u));

% find the high peformance trial index within one condition
HighInd=indtemp(find(TrialScores(indtemp)>0.6));
LowInd=indtemp(find(TrialScores(indtemp)<0.4));

% pick a condition from randomized set allPerm
conditionSelected = UniCondi(u);

% First, plot the scalp coherence
[pow,freqs,df,eppow,corr,cprod,fcoef] = allspectra(reRef_data,rate,maxfreq,HighInd,win);
% calculate the coherence
coh=abs(corr).^2;


% Plot coh of specific frequency between the two hemispheres
imagesc(squeeze(coh(46,Leftticks,Rightticks)));
xlabel('right hemisphere channels');ylabel('left hemisphere channels');
title({'45 Hz', 'coh - condition: ' num2str(conditions(u))});
colorbar;set(gca,'clim',[0 1]);colormap('jet');

% Left hemisphere
yticks([1:10]);
% Left hemisphere
yticklabels({'FP1','AF3','F11','F7','F5','F3','F1','AFF5H','AFF3H','AFF1H'});
% Right hemisphere
xticks([1:10]);
% Right hemisphere
xticklabels({'FP2','AF4','F2','F4','F6','F8','F12','AFF2H','AFF4H','AFF6H'});
xtickangle(90);

%% plot just one condtion and only the motor chans

% Left hemisphere
Leftticks=[18 19 20 27 28 29 93 94 95];
% Right hemisphere
Rightticks=[22 23 24 31 32 33 96 97 98];

u=3; % u=7; % plot just one condtion
    
indtemp=find(CondiData==UniCondi(u));

% find the high peformance trial index within one condition
HighInd=indtemp(find(TrialScores(indtemp)>0.6));
LowInd=indtemp(find(TrialScores(indtemp)<0.4));

% pick a condition from randomized set allPerm
conditionSelected = UniCondi(u);

% First, plot the scalp coherence
[pow,freqs,df,eppow,corr,cprod,fcoef] = allspectra(reRef_data,rate,maxfreq,LowInd,win);
% calculate the coherence
coh=abs(corr).^2;


% Plot coh of specific frequency between the two hemispheres
imagesc(squeeze(coh(46,Leftticks,Rightticks)));
xlabel('right hemisphere channels');ylabel('left hemisphere channels');
title({'45 Hz', 'coh - condition: ' num2str(conditions(u))});
colorbar;set(gca,'clim',[0 1]);colormap('jet');

% Left hemisphere
yticks([1:10]);
% Left hemisphere
yticklabels({'FC5','FC3','FC1','C5','C3','C1','FCC5H','FCC3H','FCC1H'});
% Right hemisphere
xticks([1:10]);
% Right hemisphere
xticklabels({'FC2','FC4','FC6','C2','C4','C6','FCC2H','FCC4H','FCC6H'});
xtickangle(90);

%% plot just one condtion and only the parietal chans

% Left hemisphere
Leftticks=[36 37 38 46 47 48 49 109 110 111];
% Right hemisphere
Rightticks=[40 41 42 51 52 53 54 112 113 114];

u=3; % u=7; % plot just one condtion
    
indtemp=find(CondiData==UniCondi(u));

% find the high peformance trial index within one condition
HighInd=indtemp(find(TrialScores(indtemp)>0.6));
LowInd=indtemp(find(TrialScores(indtemp)<0.4));

% pick a condition from randomized set allPerm
conditionSelected = UniCondi(u);

% First, plot the scalp coherence
[pow,freqs,df,eppow,corr,cprod,fcoef] = allspectra(reRef_data,rate,maxfreq,LowInd,win);
% calculate the coherence
coh=abs(corr).^2;


% Plot coh of specific frequency between the two hemispheres
imagesc(squeeze(coh(46,Leftticks,Rightticks)));
xlabel('right hemisphere channels');ylabel('left hemisphere channels');
title({'45 Hz', 'coh - condition: ' num2str(conditions(u))});
colorbar;set(gca,'clim',[0 1]);colormap('jet');

% Left hemisphere
yticks([1:10]);
% Left hemisphere
yticklabels({'CP5','CP3','CP1','P7','P5','P3','P1','CPP5H','CPP3H','CPP1H'});
% Right hemisphere
xticks([1:10]);
% Right hemisphere
xticklabels({'CP2','CP4','CP6','P2','P4','P6','P8','CPP2H','CPP4H','CPP6H'});
xtickangle(90);