% So the laplacian is in microvolts/mm^2 right now.    
% more typical in published literature would be microvolts/cm^2.   
% in which case you should multiply by 100 the raw Laplacian values before FFT
% ok multiply by 100 and you should be fine.
laplacian100_trials=laplacian_trials*100;


%%  baselinecorrect
cd /home/zhibin/Documents/GitHub/Motor_cordination/EEGanalysis/20211102
baselinesamps = 1:500; % use the first 500ms as baseline
baselinecorrected_laplacian100_trial=zeros(size(laplacian100_trials)); % initialize this matric
% Then loop through each trial for baseline correction
for i=1:size(baselinecorrected_laplacian100_trial,3) % loop through trials
    trialdata=laplacian100_trials(:,:,i);
    newtrialdata = baselinecorrect(trialdata,baselinesamps);
    baselinecorrected_laplacian100_trial(:,:,i)=newtrialdata;
end

% Just to examine the effect of baseline correction 
% (doesn't seem to do much change since my data is centered on zero, since they are average referenced)
figure;
trialdata=laplacian100_trials(:,:,1);
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
plot(laplacian100_trials(:,:,1));title('laplacian*100');
subplot(3,1,3);
plot(baselinecorrected_laplacian100_trial(:,:,1));title('laplacian after baseline correct');

%% don't baseline correct before wavelet (whatever, it doesn't change much anyway)
% baselinecorrected_laplacian100_trial=laplacian100_trials;

%% First, plot the scalp coherence
rate=Fs;
maxfreq=30;
win=501:2000; % After 500 ms
goodepochsIncluded=1:length(goodepochs);% all epochs in the matric
% First, plot the scalp coherence
% [pow,freqs,df,eppow,corr,cprod,fcoef] = allspectra(baselinecorrected_laplacian100_trial(:,1:128,:),rate,maxfreq,goodepochs,win);
[pow,freqs,df] = allspectra(baselinecorrected_laplacian100_trial,rate,maxfreq,goodepochsIncluded,win);

% plot out pow the 128 channels
AllchanNames={'FP1','FPZ','FP2','AF3','AF4','F11','F7','F5','F3','F1','FZ','F2','F4','F6','F8','F12','FT11','FC5','FC3','FC1','FCZ','FC2','FC4','FC6','FT12','T7','C5','C3','C1','CZ','C2','C4','C6','T8','TP7','CP5','CP3','CP1','CPZ','CP2','CP4','CP6','TP8','M1','M2','P7','P5','P3','P1','PZ','P2','P4','P6','P8','PO7','PO3','POZ','PO4','PO8','O1','OZ','O2','CB1','CB2','AFP1','AFP2','AF7','AF5','AFZ','AF6','AF8','AFF5H','AFF3H','AFF1H','AFF2H','AFF4H','AFF6H','F9','F10','FFT7H','FFC5H','FFC3H','FFC1H','FFC2H','FFC4H','FFC6H','FFT8H','FT9','FT7','FT8','FT10','FTT7H','FCC5H','FCC3H','FCC1H','FCC2H','FCC4H','FCC6H','FTT8H','TTP7H','CCP5H','CCP3H','CCP1H','CCP2H','CCP4H','CCP6H','TTP8H','TPP7H','CPP5H','CPP3H','CPP1H','CPP2H','CPP4H','CPP6H','TPP8H','P9','P10','PPO3H','PPO1H','PPO2H','PPO4H','PO9','PO5','PO1','PO2','PO6','PO10','CBZ','VEOG','HEOG','EMG1','EMG2','HL 1','HL 2','EMG3','EMG4','EMG5','EMG6','TRIGGER'};
% AllchanNames{1}

%% plot power spectra of all chans
for chan=1:128
    subplot(8,16,chan);
    plot(freqs,pow(:,chan));
    xlabel('freq');ylabel('pow');
    title([AllchanNames{chan}]);
end

suptitle('session2021111802')

%% Plot all powerspetra on topoplot
load('coordinates.mat');
X = [chanlocs128.X];
Y = [chanlocs128.Y];
labels = {'FP1','FPZ','FP2','AF3','AF4','F11','F7','F5','F3','F1','FZ','F2','F4','F6','F8','F12','FT11','FC5','FC3','FC1','FCZ','FC2','FC4','FC6','FT12','T7','C5','C3','C1','CZ','C2','C4','C6','T8','TP7','CP5','CP3','CP1','CPZ','CP2','CP4','CP6','TP8','M1','M2','P7','P5','P3','P1','PZ','P2','P4','P6','P8','PO7','PO3','POZ','PO4','PO8','O1','OZ','O2','CB1','CB2','AFP1','AFP2','AF7','AF5','AFZ','AF6','AF8','AFF5H','AFF3H','AFF1H','AFF2H','AFF4H','AFF6H','F9','F10','FFT7H','FFC5H','FFC3H','FFC1H','FFC2H','FFC4H','FFC6H','FFT8H','FT9','FT7','FT8','FT10','FTT7H','FCC5H','FCC3H','FCC1H','FCC2H','FCC4H','FCC6H','FTT8H','TTP7H','CCP5H','CCP3H','CCP1H','CCP2H','CCP4H','CCP6H','TTP8H','TPP7H','CPP5H','CPP3H','CPP1H','CPP2H','CPP4H','CPP6H','TPP8H','P9','P10','PPO3H','PPO1H','PPO2H','PPO4H','PO9','PO5','PO1','PO2','PO6','PO10','CBZ'};;
XX=(10*X+1)*0.455;YY=(10*Y+1)*0.455;
for chan=1:128
    subplot('Position',[1-YY(chan) XX(chan) 0.04 0.02]);
    plot(freqs,pow(:,chan));
    %xlabel('freq');ylabel('pow');
    title([AllchanNames{chan}]);
end
%suptitle('session2021111802')

%%  plot only C3/C4 28/32 F3/F4 9/13 FC3/FC4 19/23
rate=Fs;
maxfreq=50;
win=1501:2500; % During the movement, after 500 ms of planning
win=1001:1500; % planning;
win=2501:3500; % Seeing Bonus;

% plot out pow the 128 channels
AllchanNames={'FP1','FPZ','FP2','AF3','AF4','F11','F7','F5','F3','F1','FZ','F2','F4','F6','F8','F12','FT11','FC5','FC3','FC1','FCZ','FC2','FC4','FC6','FT12','T7','C5','C3','C1','CZ','C2','C4','C6','T8','TP7','CP5','CP3','CP1','CPZ','CP2','CP4','CP6','TP8','M1','M2','P7','P5','P3','P1','PZ','P2','P4','P6','P8','PO7','PO3','POZ','PO4','PO8','O1','OZ','O2','CB1','CB2','AFP1','AFP2','AF7','AF5','AFZ','AF6','AF8','AFF5H','AFF3H','AFF1H','AFF2H','AFF4H','AFF6H','F9','F10','FFT7H','FFC5H','FFC3H','FFC1H','FFC2H','FFC4H','FFC6H','FFT8H','FT9','FT7','FT8','FT10','FTT7H','FCC5H','FCC3H','FCC1H','FCC2H','FCC4H','FCC6H','FTT8H','TTP7H','CCP5H','CCP3H','CCP1H','CCP2H','CCP4H','CCP6H','TTP8H','TPP7H','CPP5H','CPP3H','CPP1H','CPP2H','CPP4H','CPP6H','TPP8H','P9','P10','PPO3H','PPO1H','PPO2H','PPO4H','PO9','PO5','PO1','PO2','PO6','PO10','CBZ','VEOG','HEOG','EMG1','EMG2','HL 1','HL 2','EMG3','EMG4','EMG5','EMG6','TRIGGER'};
% AllchanNames{1}


conditionNames={'0:4' '1:4' '1:2' '1:1' '2:1' '4:1' '4:0'}; 
%% For all 3 days
CondiData=allPerm_alldays(logical(goodepochs_alldays));
TrialScores=TrialScores_alldays(logical(goodepochs_alldays));
% plot(TrialScores,'ro');hold on;yline(median(TrialScores),'m');

%% plot power spetra in each condtion
UniCondi=unique(CondiData);

% figure;tic;% plot 7 conditions

for u=1:length(UniCondi); % sp=1:2
%     u=uall(sp);
    indtemp=find(CondiData==UniCondi(u));
    HighInd=indtemp(find(TrialScores(indtemp)>median(TrialScores(indtemp))));
    LowInd=indtemp(find(TrialScores(indtemp)<median(TrialScores(indtemp))));

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

    Ymax=0.1;
    figure(1);
    
    [pow,freqs,df,eppow,corr,cprod,fcoef] = allspectra(baselinecorrected_laplacian100_trial(:,1:128,:),rate,maxfreq,HighInd,win);
    
    subplot(1,length(UniCondi),length(UniCondi)+1-u); % C3/C4 28/32 
    
    plotx(freqs,pow(:,28),'color',[0.9290 0.6940 0.1250]);
    xlabel('freq');ylabel('pow (\muV^2/cm^2)'); xlim([0 25]);ylim([0 Ymax]);
%     title([conditionNames{UniCondi(u)} '  ' AllchanNames{28}]);
    title(['high performance trials ' conditionNames{UniCondi(u)}]);
    hold on
    
    plotx(freqs,pow(:,32),'color',[0 0.4470 0.7410]);
%     title([conditionNames{UniCondi(u)} '  ' AllchanNames{32}]);
    
    plotx(freqs,pow(:,9),'color',[1 0 1]);% F3/F4 9/13 
%     title([conditionNames{UniCondi(u)} '  ' AllchanNames{9}]);
    
    plotx(freqs,pow(:,13),'color',[0 1 0]);
%     title([conditionNames{UniCondi(u)} '  ' AllchanNames{13}]);
    
    plotx(freqs,pow(:,19),'color',[1 0 0]);% red; FC3/FC4 19/23
%     title([conditionNames{UniCondi(u)} '  ' AllchanNames{19}]);
    
    plotx(freqs,pow(:,23),'color', [0 1 1]);% cyan; FC4 
%     title([conditionNames{UniCondi(u)} '  ' AllchanNames{23}]);

    legend({'C3','C4','F3','F4','FC3','FC4'});
    
    %*********************************************************************
    
    figure(2);
    
    [pow,freqs,df,eppow,corr,cprod,fcoef] = allspectra(baselinecorrected_laplacian100_trial(:,1:128,:),rate,maxfreq,LowInd,win);
    
    subplot(1,length(UniCondi),length(UniCondi)+1-u); % C3/C4 28/32 
    
    plotx(freqs,pow(:,28),'color',[0.9290 0.6940 0.1250]);
    xlabel('freq');ylabel('pow (\muV^2/cm^2)'); xlim([0 25]);ylim([0 Ymax]);
%     title([conditionNames{UniCondi(u)} '  ' AllchanNames{28}]);
    title(['low performance trials ' conditionNames{UniCondi(u)}]);
    hold on
    
    plotx(freqs,pow(:,32),'color',[0 0.4470 0.7410]);
%     title([conditionNames{UniCondi(u)} '  ' AllchanNames{32}]);
    
    plotx(freqs,pow(:,9),'color',[1 0 1]);% F3/F4 9/13 
%     title([conditionNames{UniCondi(u)} '  ' AllchanNames{9}]);
    
    plotx(freqs,pow(:,13),'color',[0 1 0]);
%     title([conditionNames{UniCondi(u)} '  ' AllchanNames{13}]);
    
    plotx(freqs,pow(:,19),'color',[1 0 0]);% red; FC3/FC4 19/23
%     title([conditionNames{UniCondi(u)} '  ' AllchanNames{19}]);
    
    plotx(freqs,pow(:,23),'color', [0 1 1]);% cyan; FC4 
%     title([conditionNames{UniCondi(u)} '  ' AllchanNames{23}]);

    legend({'C3','C4','F3','F4','FC3','FC4'});
end
% toc;
% suptitle('Jack session 2')
% clf(figure(1));clf(figure(2)); %clear figure;
%% plot 6 channels over different conditions
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

colors=[1 1 0; 0 1 1; 0 0.4470 0.7410; 0 1 0; 0.9290 0.6940 0.1250; 1 0 1; 1 0 0];

for u=1:length(UniCondi);
    indtemp=find(CondiData==UniCondi(u));
    HighInd=indtemp(find(TrialScores(indtemp)>median(TrialScores(indtemp))));
    LowInd=indtemp(find(TrialScores(indtemp)<median(TrialScores(indtemp))));

    figure(1);
    [pow,freqs,df,eppow,corr,cprod,fcoef] = allspectra(baselinecorrected_laplacian100_trial(:,1:128,:),rate,maxfreq,HighInd,win);
    
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
    [pow,freqs,df,eppow,corr,cprod,fcoef] = allspectra(baselinecorrected_laplacian100_trial(:,1:128,:),rate,maxfreq,LowInd,win);
    
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
Ymax=0.1;

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

