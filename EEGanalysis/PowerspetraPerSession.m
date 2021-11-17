% First, plot the scalp coherence
rate=Fs;
maxfreq=50;
win=500:1500; % After 500 ms

% First, plot the scalp coherence
[pow,freqs,df,eppow,corr,cprod,fcoef] = allspectra(reRef_data,rate,maxfreq,goodepochs,win);

% plot out pow the 128 channels
AllchanNames={'FP1','FPZ','FP2','AF3','AF4','F11','F7','F5','F3','F1','FZ','F2','F4','F6','F8','F12','FT11','FC5','FC3','FC1','FCZ','FC2','FC4','FC6','FT12','T7','C5','C3','C1','CZ','C2','C4','C6','T8','TP7','CP5','CP3','CP1','CPZ','CP2','CP4','CP6','TP8','M1','M2','P7','P5','P3','P1','PZ','P2','P4','P6','P8','PO7','PO3','POZ','PO4','PO8','O1','OZ','O2','CB1','CB2','AFP1','AFP2','AF7','AF5','AFZ','AF6','AF8','AFF5H','AFF3H','AFF1H','AFF2H','AFF4H','AFF6H','F9','F10','FFT7H','FFC5H','FFC3H','FFC1H','FFC2H','FFC4H','FFC6H','FFT8H','FT9','FT7','FT8','FT10','FTT7H','FCC5H','FCC3H','FCC1H','FCC2H','FCC4H','FCC6H','FTT8H','TTP7H','CCP5H','CCP3H','CCP1H','CCP2H','CCP4H','CCP6H','TTP8H','TPP7H','CPP5H','CPP3H','CPP1H','CPP2H','CPP4H','CPP6H','TPP8H','P9','P10','PPO3H','PPO1H','PPO2H','PPO4H','PO9','PO5','PO1','PO2','PO6','PO10','CBZ','VEOG','HEOG','EMG1','EMG2','HL 1','HL 2','EMG3','EMG4','EMG5','EMG6','TRIGGER'};
% AllchanNames{1}

%% all chans
for chan=1:128
    subplot(8,16,chan);
    plot(freqs,pow(:,chan));
    xlabel('freq');ylabel('pow');
    title([AllchanNames{chan}]);
end

suptitle('zhibin session 3 on 11/06/2021')
%%  C3/C4 28/32 F3/F4 9/13 FC3/FC4 19/23

% plot power spetra in each condtion

UniCondi=unique(CondiData);

figure;tic;% plot 7 conditions
for u=1:length(UniCondi); % sp=1:2
%     u=uall(sp);
    indtemp=find(CondiData==UniCondi(u));
    HighInd=indtemp(find(TrialScores(indtemp)>0.8));

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
    
    [pow,freqs,df,eppow,corr,cprod,fcoef] = allspectra(reRef_data,rate,maxfreq,HighInd,win);
    
    subplot(1,length(UniCondi),length(UniCondi)+1-u); % C3/C4 28/32 
    
    plotx(freqs,pow(:,28),'color',[0.9290 0.6940 0.1250]);
    xlabel('freq');ylabel('pow');ylim([0 0.3*10e-5]);xlim([0 50]);
%     title([conditionNames{UniCondi(u)} '  ' AllchanNames{28}]);
    title([conditionNames{UniCondi(u)}]);
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
toc;
suptitle('zhibin session 3 on 11/06/2021')

%% plot 6 channels over different conditions
UniCondi=unique(CondiData);

% figure;tic;
% for i=1:6
%     subplot(1,6,i);
% end

ax1=subplot(1,6,1);
ax2=subplot(1,6,2);
ax3=subplot(1,6,3);
ax4=subplot(1,6,4);
ax5=subplot(1,6,5);
ax6=subplot(1,6,6);

colors=[1 1 0; 0 1 1; 0 0.4470 0.7410; 0 1 0; 0.9290 0.6940 0.1250; 1 0 1; 1 0 0];

for u=1:length(UniCondi);
    indtemp=find(CondiData==UniCondi(u));
    HighInd=indtemp(find(TrialScores(indtemp)>0.8));

    [pow,freqs,df,eppow,corr,cprod,fcoef] = allspectra(reRef_data,rate,maxfreq,HighInd,win);
    
    ConditionColor=colors(u,:);

    axes(ax1); % C3/C4 28/32 
    hold on;
    plotx(freqs,pow(:,28),'color',ConditionColor);
%     legend{u}=(conditionNames{UniCondi(u)});
    hold off

    axes(ax2); 
    hold on;
    plotx(freqs,pow(:,32),'color',ConditionColor);
    hold off

    axes(ax3); 
    hold on;
    plotx(freqs,pow(:,9),'color',ConditionColor);% F3/F4 9/13 
    hold off

    axes(ax4); 
    hold on;
    plotx(freqs,pow(:,13),'color',ConditionColor);
    hold off

    axes(ax5); 
    hold on;
    plotx(freqs,pow(:,19),'color',ConditionColor);% red; FC3/FC4 19/23
    hold off

    axes(ax6); 
    hold on;
    plotx(freqs,pow(:,23),'color',ConditionColor);% cyan; FC4 
    hold off

end
toc;
suptitle('zhibin session 3 on 11/06/2021')

axes(ax1);
xlabel('freq');ylabel('pow');ylim([0 0.3*10e-5]);xlim([0 50]);title([AllchanNames{28}]); %C3
legend(conditionNames);
axes(ax2);
xlabel('freq');ylabel('pow');ylim([0 0.3*10e-5]);xlim([0 50]);title([AllchanNames{32}]); %C4
legend(conditionNames);
axes(ax3);
xlabel('freq');ylabel('pow');ylim([0 0.3*10e-5]);xlim([0 50]);title([AllchanNames{9}]); %F3
legend(conditionNames);
axes(ax4);
xlabel('freq');ylabel('pow');ylim([0 0.3*10e-5]);xlim([0 50]);title([AllchanNames{13}]); %F4
legend(conditionNames);
axes(ax5);
xlabel('freq');ylabel('pow');ylim([0 0.3*10e-5]);xlim([0 50]);title([AllchanNames{19}]); %FC3
legend(conditionNames);
axes(ax6);
xlabel('freq');ylabel('pow');ylim([0 0.3*10e-5]);xlim([0 50]);title([AllchanNames{23}]); %FC4
legend(conditionNames);
    