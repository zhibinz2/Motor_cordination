% So the laplacian is in microvolts/mm^2 right now.    
% more typical in published literature would be microvolts/cm^2.   
% in which case you should multiply by 100 the raw Laplacian values before FFT
% ok multiply by 100 and you should be fine.

% laplacian100_trials=laplacian_trials*100;
% Replace laplacian100_trials with afterICA_trials hereafter

%%  baselinecorrect
cd /home/zhibin/Documents/GitHub/Motor_cordination/EEGanalysis/20211102
baselinesamps = 1:500; % use the first 500ms as baseline
baselinecorrected_trial=zeros(size(afterICA_trials)); % initialize this matric
% Then loop through each trial for baseline correction
for i=1:size(baselinecorrected_trial,3) % loop through trials
    trialdata=afterICA_trials(:,:,i);
    newtrialdata = baselinecorrect(trialdata,baselinesamps);
    baselinecorrected_trial(:,:,i)=newtrialdata;
end

% Just to examine the effect of baseline correction 
% (doesn't seem to do much change since my data is centered on zero, since they are average referenced)
figure;
trialdata=afterICA_trials(:,:,1);
baseline = mean(trialdata(baselinesamps,:),1);
newtrialdata = trialdata - ones(size(trialdata,1),1)*baseline;
subplot(2,1,1);
plot(trialdata);title('before baseline correct');
subplot(2,1,2);
plot(newtrialdata);title('after baseline correct');

% Just to examine the values of laplacian and baseline correction
% (doesn't seem change much either)
figure;
subplot(3,1,1);
plot(laplacian_trials(:,:,1));title('laplacian');
subplot(3,1,2);
plot(afterICA_trials(:,:,1));title('laplacian*100');
subplot(3,1,3);
plot(baselinecorrected_trial(:,:,1));title('laplacian after baseline correct');

%% don't baseline correct before wavelet (whatever, it doesn't change much anyway)
% baselinecorrected_laplacian100_trial=afterICA_trials;

%% First, settings for plotting
rate=Fs;
maxfreq=30;
win=1001:2000; % After 500 ms
% goodepochsIncluded=1:length(goodepochs);% all epochs in the matric
% First, plot the scalp coherence
% [pow,freqs,df,eppow,corr,cprod,fcoef] = allspectra(baselinecorrected_laplacian100_trial(:,1:128,:),rate,maxfreq,goodepochs,win);
% [pow,freqs,df] = allspectra(baselinecorrected_trial,rate,maxfreq,goodepochsIncluded,win);

% plot out pow the 128 channels
AllchanNames={'FP1','FPZ','FP2','AF3','AF4','F11','F7','F5','F3','F1','FZ','F2','F4','F6','F8','F12','FT11','FC5','FC3','FC1','FCZ','FC2','FC4','FC6','FT12','T7','C5','C3','C1','CZ','C2','C4','C6','T8','TP7','CP5','CP3','CP1','CPZ','CP2','CP4','CP6','TP8','M1','M2','P7','P5','P3','P1','PZ','P2','P4','P6','P8','PO7','PO3','POZ','PO4','PO8','O1','OZ','O2','CB1','CB2','AFP1','AFP2','AF7','AF5','AFZ','AF6','AF8','AFF5H','AFF3H','AFF1H','AFF2H','AFF4H','AFF6H','F9','F10','FFT7H','FFC5H','FFC3H','FFC1H','FFC2H','FFC4H','FFC6H','FFT8H','FT9','FT7','FT8','FT10','FTT7H','FCC5H','FCC3H','FCC1H','FCC2H','FCC4H','FCC6H','FTT8H','TTP7H','CCP5H','CCP3H','CCP1H','CCP2H','CCP4H','CCP6H','TTP8H','TPP7H','CPP5H','CPP3H','CPP1H','CPP2H','CPP4H','CPP6H','TPP8H','P9','P10','PPO3H','PPO1H','PPO2H','PPO4H','PO9','PO5','PO1','PO2','PO6','PO10','CBZ','VEOG','HEOG','EMG1','EMG2','HL 1','HL 2','EMG3','EMG4','EMG5','EMG6','TRIGGER'};
% AllchanNames{1}

%% plot power spectra of all chans (Skip)
for chan=1:128
    subplot(8,16,chan);
    plot(freqs,pow(:,chan));
    xlabel('freq');ylabel('pow');
    ylim([0 0.5]);
    title([AllchanNames{chan}]);
end

suptitle('session2021111802')


%% Calculate subplot coordinates for scalp map
cd /home/zhibin/Documents/GitHub/Motor_cordination/EEGanalysis
load('coordinates.mat');
X = [chanlocs128.X];
Y = [chanlocs128.Y];
Z = [chanlocs128.Z];
labels = {'FP1','FPZ','FP2','AF3','AF4','F11','F7','F5','F3','F1','FZ','F2','F4','F6','F8','F12','FT11','FC5','FC3','FC1','FCZ','FC2','FC4','FC6','FT12','T7','C5','C3','C1','CZ','C2','C4','C6','T8','TP7','CP5','CP3','CP1','CPZ','CP2','CP4','CP6','TP8','M1','M2','P7','P5','P3','P1','PZ','P2','P4','P6','P8','PO7','PO3','POZ','PO4','PO8','O1','OZ','O2','CB1','CB2','AFP1','AFP2','AF7','AF5','AFZ','AF6','AF8','AFF5H','AFF3H','AFF1H','AFF2H','AFF4H','AFF6H','F9','F10','FFT7H','FFC5H','FFC3H','FFC1H','FFC2H','FFC4H','FFC6H','FFT8H','FT9','FT7','FT8','FT10','FTT7H','FCC5H','FCC3H','FCC1H','FCC2H','FCC4H','FCC6H','FTT8H','TTP7H','CCP5H','CCP3H','CCP1H','CCP2H','CCP4H','CCP6H','TTP8H','TPP7H','CPP5H','CPP3H','CPP1H','CPP2H','CPP4H','CPP6H','TPP8H','P9','P10','PPO3H','PPO1H','PPO2H','PPO4H','PO9','PO5','PO1','PO2','PO6','PO10','CBZ'};

sort(Z) %just to examine
ZZ=((max(Z)-min(Z))*1.5-(Z-min(Z)))./((max(Z)-min(Z))*1.5); %Make higher Z smaller fraction
sort(ZZ) %just to examine

XXPLOT=(10*Y)*0.95;
YYPLOT=(10*X)*0.8;

XXPLOT=XXPLOT.*ZZ;
YYPLOT=YYPLOT.*ZZ;

XXPLOT=-1*XXPLOT;
YYPLOT=YYPLOT;

XXPLOT=XXPLOT*0.6;
YYPLOT=YYPLOT*0.7;

XXPLOT=XXPLOT+0.5;
YYPLOT=YYPLOT+0.5;

plot(XXPLOT,YYPLOT,'ro');

labels = {'FP1','FPZ','FP2','AF3','AF4','F11','F7','F5','F3','F1','FZ','F2','F4','F6','F8','F12','FT11','FC5','FC3','FC1','FCZ','FC2','FC4','FC6','FT12','T7','C5','C3','C1','CZ','C2','C4','C6','T8','TP7','CP5','CP3','CP1','CPZ','CP2','CP4','CP6','TP8','M1','M2','P7','P5','P3','P1','PZ','P2','P4','P6','P8','PO7','PO3','POZ','PO4','PO8','O1','OZ','O2','CB1','CB2','AFP1','AFP2','AF7','AF5','AFZ','AF6','AF8','AFF5H','AFF3H','AFF1H','AFF2H','AFF4H','AFF6H','F9','F10','FFT7H','FFC5H','FFC3H','FFC1H','FFC2H','FFC4H','FFC6H','FFT8H','FT9','FT7','FT8','FT10','FTT7H','FCC5H','FCC3H','FCC1H','FCC2H','FCC4H','FCC6H','FTT8H','TTP7H','CCP5H','CCP3H','CCP1H','CCP2H','CCP4H','CCP6H','TTP8H','TPP7H','CPP5H','CPP3H','CPP1H','CPP2H','CPP4H','CPP6H','TPP8H','P9','P10','PPO3H','PPO1H','PPO2H','PPO4H','PO9','PO5','PO1','PO2','PO6','PO10','CBZ'};;
text(XXPLOT,YYPLOT,labels,'VerticalAlignment','bottom','HorizontalAlignment','right');

%% Select good trials and condition
UniCondi=unique(CondiData);

% Use Bgoodepochs to select good trials only
CondiDataGoodTrials=CondiData(goodepochs);
allPermGoodTrials=allPerm(goodepochs);
TrialScoresGoodTrials=TrialScores(goodepochs);

% win=501:1000; % planning
% win=1001:2000; % movement
% u=1;
% u=7;

%% Plot all powerspetra on scalp map (skip)
indtemp=find(CondiDataGoodTrials==UniCondi(u));
HighInd=indtemp(find(TrialScoresGoodTrials(indtemp)>median(TrialScoresGoodTrials(indtemp))));
[pow,freqs,df] = allspectra(baselinecorrected_trial(win,:,HighInd),rate,maxfreq);

figure('units','normalized','outerposition',[0 0 1 1]);
for chan=1:128
    subplot('Position',[XXPLOT(chan) YYPLOT(chan) 0.04 0.03]);
    plot(freqs,pow(:,chan));
    if ~isempty(find([1:127]==chan))
    set(gca,'XTick',[]); set(gca,'YTick',[]); 
    end
    if chan==128
        xlabel('freq');ylabel('pow');
    end
    %ylim([0 1]);
    title([AllchanNames{chan}]);
end
% suptitle('session2021111802')

%% Plot all powerspetra on scalp map for each condition 
conditionNames={'0:4' '1:4' '1:2' '1:1' '2:1' '4:1' '4:0'}; 
Ymax=2;

% all performance trials
figure('units','normalized','outerposition',[0 0 1 1]);
for u=[1 7]; % sp=1:2
    ConditionColor=colors(u,:);
    indtemp=find(CondiDataGoodTrials==UniCondi(u));
    
    [pow,freqs,df] = allspectra(baselinecorrected_trial(win,:,indtemp),rate,maxfreq);

    for chan=1:128
        subplot('Position',[XXPLOT(chan) YYPLOT(chan) 0.02 0.03]);
        
        plot(freqs,pow(:,chan),'color',ConditionColor);
        hold on
        if ~isempty(find([1:127]==chan))
        set(gca,'XTick',[]); set(gca,'YTick',[]); 
        end
        if chan==128
            xlabel('freq');ylabel('pow');
        end
        ylim([0 Ymax]);
        title([AllchanNames{chan}]);
    end
end
hold off;
legend((conditionNames([1 7])));




% High performance trials
figure('units','normalized','outerposition',[0 0 0.6 1]);
for u=[1 7]; % sp=1:2
    ConditionColor=colors(u,:);
    indtemp=find(CondiDataGoodTrials==UniCondi(u));
    HighInd=indtemp(find(TrialScoresGoodTrials(indtemp)>median(TrialScoresGoodTrials(indtemp))));
    LowInd=indtemp(find(TrialScoresGoodTrials(indtemp)<median(TrialScoresGoodTrials(indtemp))));
    
    [pow,freqs,df] = allspectra(baselinecorrected_trial(win,:,HighInd),rate,maxfreq);

    for chan=1:128
        subplot('Position',[XXPLOT(chan) YYPLOT(chan) 0.02 0.03]);
        
        plot(freqs,pow(:,chan),'color',ConditionColor);
        hold on
        if ~isempty(find([1:127]==chan))
        set(gca,'XTick',[]); set(gca,'YTick',[]); 
        end
        if chan==128
            xlabel('freq');ylabel('pow');
        end
        ylim([0 Ymax]);
        title([AllchanNames{chan}]);
    end
end
hold off;
legend((conditionNames([1 7])));
suptitle('high performance trials');

% Low performance trials
figure('units','normalized','outerposition',[0 0 0.6 1]);
for u=[1 7]; % sp=1:2
    ConditionColor=colors(u,:);
    indtemp=find(CondiDataGoodTrials==UniCondi(u));
    HighInd=indtemp(find(TrialScoresGoodTrials(indtemp)>median(TrialScoresGoodTrials(indtemp))));
    LowInd=indtemp(find(TrialScoresGoodTrials(indtemp)<median(TrialScoresGoodTrials(indtemp))));
    
    [pow,freqs,df] = allspectra(baselinecorrected_trial(win,:,LowInd),rate,maxfreq);

    for chan=1:128
        subplot('Position',[XXPLOT(chan) YYPLOT(chan) 0.02 0.03]);
        
        plot(freqs,pow(:,chan),'color',ConditionColor);
        hold on
        if ~isempty(find([1:127]==chan))
        set(gca,'XTick',[]); set(gca,'YTick',[]); 
        end
        if chan==128
            xlabel('freq');ylabel('pow');
        end
        ylim([0 Ymax]);
        title([AllchanNames{chan}]);
    end
end
hold off;
legend((conditionNames([1 7])));
suptitle('low performance trials');
%% Plot all erp of planning on scalp map for each condition (skip)
indtemp=find(CondiDataGoodTrials==UniCondi(u));
HighInd=indtemp(find(TrialScoresGoodTrials(indtemp)>median(TrialScoresGoodTrials(indtemp))));

figure('units','normalized','outerposition',[0 0 1 1]);
for chan=1:128
    subplot('Position',[XXPLOT(chan) YYPLOT(chan) 0.04 0.03]);
    plot(1:length(win),mean(baselinecorrected_trial(win,chan,HighInd),3));
    if ~isempty(find([1:127]==chan))
    set(gca,'XTick',[]); set(gca,'YTick',[]); 
    end
    if chan==128
        xlabel('time');ylabel('laplacian power');
    end
    ylim([-.5 .5]);
    title([AllchanNames{chan}]);
end
%suptitle('session2021111802')


%%  plot only C3/C4 28/32 F3/F4 9/13 FC3/FC4 19/23 (skip)
rate=Fs;
maxfreq=50;
win=1501:2500; % During the movement, after 500 ms of planning
win=1001:1500; % planning;
win=2501:3500; % Seeing Bonus;

% plot out pow the 128 channels
AllchanNames={'FP1','FPZ','FP2','AF3','AF4','F11','F7','F5','F3','F1','FZ','F2','F4','F6','F8','F12','FT11','FC5','FC3','FC1','FCZ','FC2','FC4','FC6','FT12','T7','C5','C3','C1','CZ','C2','C4','C6','T8','TP7','CP5','CP3','CP1','CPZ','CP2','CP4','CP6','TP8','M1','M2','P7','P5','P3','P1','PZ','P2','P4','P6','P8','PO7','PO3','POZ','PO4','PO8','O1','OZ','O2','CB1','CB2','AFP1','AFP2','AF7','AF5','AFZ','AF6','AF8','AFF5H','AFF3H','AFF1H','AFF2H','AFF4H','AFF6H','F9','F10','FFT7H','FFC5H','FFC3H','FFC1H','FFC2H','FFC4H','FFC6H','FFT8H','FT9','FT7','FT8','FT10','FTT7H','FCC5H','FCC3H','FCC1H','FCC2H','FCC4H','FCC6H','FTT8H','TTP7H','CCP5H','CCP3H','CCP1H','CCP2H','CCP4H','CCP6H','TTP8H','TPP7H','CPP5H','CPP3H','CPP1H','CPP2H','CPP4H','CPP6H','TPP8H','P9','P10','PPO3H','PPO1H','PPO2H','PPO4H','PO9','PO5','PO1','PO2','PO6','PO10','CBZ','VEOG','HEOG','EMG1','EMG2','HL 1','HL 2','EMG3','EMG4','EMG5','EMG6','TRIGGER'};
% AllchanNames{1}

conditionNames={'0:4' '1:4' '1:2' '1:1' '2:1' '4:1' '4:0'}; 

%% For all 3 days (skip)
CondiData=allPerm_alldays(logical(goodepochs_alldays));
TrialScores=TrialScores_alldays(logical(goodepochs_alldays));
% plot(TrialScores,'ro');hold on;yline(median(TrialScores),'m');

%% plot power spetra in each condtion (skip)
UniCondi=unique(CondiData);

% Use Bgoodepochs to select good trials only
CondiDataGoodTrials=CondiData(goodepochs);
allPermGoodTrials=allPerm(goodepochs);
TrialScoresGoodTrials=TrialScores(goodepochs);
colors=[0 0 1; 0 0.4 0.85; 0 0.8 0.7; 0 1 0; 0.7 0.8 0; 0.85 0.4 0; 1 0 0];


win=501:1000; % planning
win=1001:2000; % movement

% figure;tic;% plot 7 conditions

for u=1:length(UniCondi); % sp=1:2
%     u=uall(sp);
    indtemp=find(CondiDataGoodTrials==UniCondi(u));
    HighInd=indtemp(find(TrialScoresGoodTrials(indtemp)>median(TrialScoresGoodTrials(indtemp))));
    LowInd=indtemp(find(TrialScoresGoodTrials(indtemp)<median(TrialScoresGoodTrials(indtemp))));

%     [out,idx] = sort(TrialScores(indtemp)); % example: [out,idx] = sort([14 8 91 19])
%     
%     % Average every n elements
%     n=10;
%     % OutMean = arrayfun(@(i) mean(out(i:i+n-1)),1:n:length(out)-n+1);
%     OutMean = arrayfun(@(i) mean(out(i:i+n-1)),1:1:length(out)-n+1);
%     OutStd = arrayfun(@(i) std(out(i:i+n-1)),1:1:length(out)-n+1);
%     
%     
%     CohMean=zeros(1,length(out)-n+1);

%     for i=1:1:length(out)-n+1
%         % First, plot the scalp coherence
%         [pow,freqs,df,eppow,corr,cprod,fcoef] = allspectra(reRef_data,rate,maxfreq,idx(i:i+n-1),win);
%         % calculate the coherence
%         coh=abs(corr).^2;
%         % Plot coh of specific frequency between the two hemispheres
%         CohMean(i)=mean(squeeze(coh(FreqI+1,Leftticks,Rightticks)),'all');
%     end
    
    %*********************************************************************

    Ymax=2;
    figure(1);
    
    [pow,freqs,df,eppow,corr,cprod,fcoef] = allspectra(baselinecorrected_trial,rate,maxfreq,HighInd,win);
    
    subplot(1,length(UniCondi),length(UniCondi)+1-u); % C3/C4 28/32 
    
%     plotx(freqs,pow(:,28),'color',[0.9290 0.6940 0.1250]);
    plotx(freqs,pow(:,28),'color', colors(1,:));
    xlabel('freq');ylabel('pow (\muV^2/cm^2)'); xlim([0 25]);ylim([0 Ymax]);
%     title([conditionNames{UniCondi(u)} '  ' AllchanNames{28}]);
    title(['high performance trials ' conditionNames{UniCondi(u)}]);
    hold on
    
%     plotx(freqs,pow(:,32),'color',[0 0.4470 0.7410]);
    plotx(freqs,pow(:,32),'color', colors(7,:));
%     title([conditionNames{UniCondi(u)} '  ' AllchanNames{32}]);
    
%     plotx(freqs,pow(:,9),'color',[1 0 1]);% F3/F4 9/13 
    plotx(freqs,pow(:,9),'color', colors(2,:));
%     title([conditionNames{UniCondi(u)} '  ' AllchanNames{9}]);
    
%     plotx(freqs,pow(:,13),'color',[0 1 0]);
    plotx(freqs,pow(:,13),'color', colors(6,:));
%     title([conditionNames{UniCondi(u)} '  ' AllchanNames{13}]);
    
%     plotx(freqs,pow(:,19),'color',[1 0 0]);% red; FC3/FC4 19/23
    plotx(freqs,pow(:,19),'color', colors(3,:));
%     title([conditionNames{UniCondi(u)} '  ' AllchanNames{19}]);
    
%     plotx(freqs,pow(:,23),'color', [0 1 1]);% cyan; FC4 
    plotx(freqs,pow(:,23),'color', colors(5,:));
%     title([conditionNames{UniCondi(u)} '  ' AllchanNames{23}]);

    legend({'C3','C4','F3','F4','FC3','FC4'});
    
    %*********************************************************************
    
    figure(2);
    
    [pow,freqs,df,eppow,corr,cprod,fcoef] = allspectra(baselinecorrected_trial,rate,maxfreq,LowInd,win);
    
    subplot(1,length(UniCondi),length(UniCondi)+1-u); % C3/C4 28/32 
    
%     plotx(freqs,pow(:,28),'color',[0.9290 0.6940 0.1250]);
    plotx(freqs,pow(:,28),'color', colors(1,:));
    xlabel('freq');ylabel('pow (\muV^2/cm^2)'); xlim([0 25]);ylim([0 Ymax]);
%     title([conditionNames{UniCondi(u)} '  ' AllchanNames{28}]);
    title(['low performance trials ' conditionNames{UniCondi(u)}]);
    hold on
    
%     plotx(freqs,pow(:,32),'color',[0 0.4470 0.7410]);
    plotx(freqs,pow(:,32),'color', colors(7,:));
%     title([conditionNames{UniCondi(u)} '  ' AllchanNames{32}]);
    
%     plotx(freqs,pow(:,9),'color',[1 0 1]);% F3/F4 9/13 
    plotx(freqs,pow(:,9),'color', colors(2,:));
%     title([conditionNames{UniCondi(u)} '  ' AllchanNames{9}]);
    
%     plotx(freqs,pow(:,13),'color',[0 1 0]);
    plotx(freqs,pow(:,13),'color', colors(6,:));
%     title([conditionNames{UniCondi(u)} '  ' AllchanNames{13}]);
    
%     plotx(freqs,pow(:,19),'color',[1 0 0]);% red; FC3/FC4 19/23
    plotx(freqs,pow(:,19),'color', colors(3,:));
%     title([conditionNames{UniCondi(u)} '  ' AllchanNames{19}]);
    
%     plotx(freqs,pow(:,23),'color', [0 1 1]);% cyan; FC4 
    plotx(freqs,pow(:,23),'color', colors(5,:));
%     title([conditionNames{UniCondi(u)} '  ' AllchanNames{23}]);

    legend({'C3','C4','F3','F4','FC3','FC4'});
end
% toc;
% suptitle('Jack session 2')
% clf(figure(1));clf(figure(2)); %clear figure;
%% plot 6 channels over different conditions (skip)
win=1001:1500; % planning;
win=1501:2500; % During the movement, after 500 ms of planning
win=2501:3500; % Seeing Bonus;

UniCondi=unique(CondiData);

% figure;tic;
% for i=1:6
%     subplot(1,6,i);
% end

% clf(figure(1));clf(figure(2));

figure(1);
ax1=subplot(3,2,1);
ax2=subplot(3,2,2);
ax3=subplot(3,2,3);
ax4=subplot(3,2,4);
ax5=subplot(3,2,5);
ax6=subplot(3,2,6);

figure(2);
ax12=subplot(3,2,1);
ax22=subplot(3,2,2);
ax32=subplot(3,2,3);
ax42=subplot(3,2,4);
ax52=subplot(3,2,5);
ax62=subplot(3,2,6);

% colors=[1 0 1; 0 1 1; 0 0.4470 0.7410; 0 1 0; 0.9290 0.6940 0.1250; 1 0 1; 1 0 0];
colors=[0 0 1; 0 0.4 0.85; 0 0.8 0.7; 0 1 0; 0.7 0.8 0; 0.85 0.4 0; 1 0 0];

for u=1:length(UniCondi); % u=2;
    indtemp=find(CondiDataGoodTrials==UniCondi(u));
    HighInd=indtemp(find(TrialScoresGoodTrials(indtemp)>median(TrialScoresGoodTrials(indtemp))));
    LowInd=indtemp(find(TrialScoresGoodTrials(indtemp)<median(TrialScoresGoodTrials(indtemp))));

    figure(1);
    [pow,freqs,df,eppow,corr,cprod,fcoef] = allspectra(baselinecorrected_trial,rate,maxfreq,HighInd,win);
    
    ConditionColor=colors(u,:);

    axes(ax1); 
    hold on;
    plotx(freqs,pow(:,9),'color',ConditionColor);% F3/F4 9/13 
    hold off

    axes(ax2); 
    hold on;
    plotx(freqs,pow(:,13),'color',ConditionColor);
    hold off

    axes(ax3); 
    hold on;
    plotx(freqs,pow(:,19),'color',ConditionColor);% red; FC3/FC4 19/23
    hold off

    axes(ax4); 
    hold on;
    plotx(freqs,pow(:,23),'color',ConditionColor);% cyan; FC4 
    hold off
    
    axes(ax5); % C3/C4 28/32 
    hold on;
    plotx(freqs,pow(:,28),'color',ConditionColor);
%     legend{u}=(conditionNames{UniCondi(u)});
    hold off

    axes(ax6); 
    hold on;
    plotx(freqs,pow(:,32),'color',ConditionColor);
    hold off

    %*********************************************************************
    figure(2);
    [pow,freqs,df,eppow,corr,cprod,fcoef] = allspectra(baselinecorrected_trial,rate,maxfreq,LowInd,win);
    
    ConditionColor=colors(u,:);

    axes(ax12); 
    hold on;
    plotx(freqs,pow(:,9),'color',ConditionColor);% F3/F4 9/13 
    hold off

    axes(ax22); 
    hold on;
    plotx(freqs,pow(:,13),'color',ConditionColor);
    hold off

    axes(ax32); 
    hold on;
    plotx(freqs,pow(:,19),'color',ConditionColor);% red; FC3/FC4 19/23
    hold off

    axes(ax42); 
    hold on;
    plotx(freqs,pow(:,23),'color',ConditionColor);% cyan; FC4 
    hold off
    
    axes(ax52); % C3/C4 28/32 
    hold on;
    plotx(freqs,pow(:,28),'color',ConditionColor);
%     legend{u}=(conditionNames{UniCondi(u)});
    hold off

    axes(ax62); 
    hold on;
    plotx(freqs,pow(:,32),'color',ConditionColor);
    hold off

end
% toc;
% suptitle('Jack session 2')
Ymax=1;

figure(1);
axes(ax1);
xlabel('freq');ylabel('pow (\muV^2/cm^2)');ylim([0 Ymax]);xlim([0 25]);title(['high ' AllchanNames{9}]); %F3
legend(conditionNames);
axes(ax2);
xlabel('freq');ylabel('pow (\muV^2/cm^2)');ylim([0 Ymax]);xlim([0 25]);title(['high ' AllchanNames{13}]); %F4
legend(conditionNames);
axes(ax3);
xlabel('freq');ylabel('pow (\muV^2/cm^2)');ylim([0 Ymax]);xlim([0 25]);title(['high ' AllchanNames{19}]); %FC3
legend(conditionNames);
axes(ax4);
xlabel('freq');ylabel('pow (\muV^2/cm^2)');ylim([0 Ymax]);xlim([0 25]);title(['high ' AllchanNames{23}]); %FC4
legend(conditionNames);
axes(ax5);
xlabel('freq');ylabel('pow (\muV^2/cm^2)');ylim([0 Ymax]);xlim([0 25]);title(['high ' AllchanNames{28}]); %C3
legend(conditionNames);
axes(ax6);
xlabel('freq');ylabel('pow (\muV^2/cm^2)');ylim([0 Ymax]);xlim([0 25]);title(['high ' AllchanNames{32}]); %C4
legend(conditionNames);


figure(2);
axes(ax12);
xlabel('freq');ylabel('pow (\muV^2/cm^2)');ylim([0 Ymax]);xlim([0 25]);title(['low ' AllchanNames{9}]); %F3
legend(conditionNames);
axes(ax22);
xlabel('freq');ylabel('pow (\muV^2/cm^2)');ylim([0 Ymax]);xlim([0 25]);title(['low ' AllchanNames{13}]); %F4
legend(conditionNames);
axes(ax32);
xlabel('freq');ylabel('pow (\muV^2/cm^2)');ylim([0 Ymax]);xlim([0 25]);title(['low ' AllchanNames{19}]); %FC3
legend(conditionNames);
axes(ax42);
xlabel('freq');ylabel('pow (\muV^2/cm^2)');ylim([0 Ymax]);xlim([0 25]);title(['low ' AllchanNames{23}]); %FC4
legend(conditionNames);
axes(ax52);
xlabel('freq');ylabel('pow (\muV^2/cm^2)');ylim([0 Ymax]);xlim([0 25]);title(['low ' AllchanNames{28}]); %C3
legend(conditionNames);
axes(ax62);
xlabel('freq');ylabel('pow (\muV^2/cm^2)');ylim([0 Ymax]);xlim([0 25]);title(['low ' AllchanNames{32}]); %C4
legend(conditionNames);


%% examine 1:2 data and power of FC3 19 AND C3 28 (skip)
conditionNames;
u=3;% condition 1:2
indtemp=find(CondiDataGoodTrials==UniCondi(u));
HighInd=indtemp(find(TrialScoresGoodTrials(indtemp)>median(TrialScoresGoodTrials(indtemp))));
LowInd=indtemp(find(TrialScoresGoodTrials(indtemp)<median(TrialScoresGoodTrials(indtemp))));
ConnectedData=baselinecorrected_trial(:,:,LowInd(1));
% zeros(size(baselinecorrected_laplacian100_trial,1)*length(LowInd),128);
for i=2:length(LowInd)
    ConnectedData=cat(1,ConnectedData,baselinecorrected_trial(:,:,LowInd(i)));
end
figure;
subplot(2,1,1);
plot(1:size(ConnectedData,1),ConnectedData(:,:));%ylim([-20 20]);
xlines=2000:2000:size(ConnectedData,1);
hold on;
for i=1:length(xlines)
    xline(xlines(i));
end
hold off;
title('All channels in 1:2 low performance trials');
subplot(2,1,2);
plot(1:size(ConnectedData,1),ConnectedData(:,[19 28]));%ylim([-20 20]);
hold on;
for i=1:length(xlines)
    xline(xlines(i));
end
title('channel FC3 C3 in 1:2 low performance trials');
legend({'FC3','C3'});
hold off;

% So it is trial number 3 9, 21 77,
TrialScoresGoodTrials(LowInd(3))
TrialScoresGoodTrials(LowInd(9))
% Scores are not too bad

[pow,freqs,df,eppow,corr,cprod,fcoef] = allspectra(baselinecorrected_trial(win,:,LowInd),rate,maxfreq);
figure;
plot(freqs,pow(:,[19 28]));