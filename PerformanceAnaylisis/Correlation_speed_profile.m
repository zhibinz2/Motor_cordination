% This is after plotSpectrogram.m
%% Correlation of speed profile
%% XXX1
% time x trials 
ALLhL=ALLhL(goodepochs,:)';ALLhR=ALLhR(goodepochs,:)';% travel distance
ALLvL=ALLvL(goodepochs,:)';ALLvR=ALLvR(goodepochs,:)';% velocity
ALLaL=ALLaL(goodepochs,:)';ALLaR=ALLaR(goodepochs,:)';% acceleration

%% XXX2
% Structure of this matrix: 
% wfreq x time x chans x trials
%    10   1500    128     289
size(normPowcnorm_ALLchan_sgolay_ALLtrials)

% select a frequency to analysis
wfreq; % checking
FreqSelected=3; % theta =6
CondiDataGoodTrials;

% select a trial
i=1; % which trial
theta_trial_spectrogram=squeeze(normPowcnorm_ALLchan_sgolay_ALLtrials(FreqSelected,:,:,i)); % time x chans

%% Examine correlation in this one trial - left hand (skip)
clear corr
XXX1=ALLvL(:,1);XXX2=theta_trial_spectrogram(1:size(ALLvL,1),:);
[RHO,PVAL] = corr(XXX1,XXX2);

figure;
imagesc(RHO);colorbar;
xlabel('channels');ylabel('trials'):title('Correlation coefficent between velocity and spectrogram');

figure;
topoplot(mean(RHO,1)',test,'nosedir','+Y');ColorLim=1;
colormap jet; colorbar ; caxis([-1*ColorLim ColorLim]);
title('Correlation coefficent between velocity and spectrogram');
%% same as above - right hand (skip)
clear RHO PVAL corr
XXX1=ALLvR(:,1);XXX2=theta_trial_spectrogram(1:size(ALLvR,1),:);
[RHO,PVAL] = corr(XXX1,XXX2);

figure;
imagesc(RHO);colorbar;
xlabel('channels');ylabel('trials');title('Correlation coefficent between velocity and spectrogram');

figure;
topoplot(mean(RHO,1)',test,'nosedir','+Y');ColorLim=1;
colormap jet; colorbar ; caxis([-1*ColorLim ColorLim]);
title('Correlation coefficent between velocity and spectrogram');

%% compare speed profile of left and right hand in all good trials (skip)
subplot(1,2,1);imagesc(ALLvL(:,CondiDataGoodTrials));
subplot(1,2,2);imagesc(ALLvR(:,CondiDataGoodTrials));
%% compare speed profile of left and right hand in one trial (skip)
subplot(1,2,1);plot(ALLvL(:,CondiDataGoodTrials(i)));
subplot(1,2,2);plot(ALLvR(:,CondiDataGoodTrials(i)));
%% correlation of left hand in one trial (skip)
clear corr
XXX1=ALLaL(:,i);XXX2=theta_trial_spectrogram(1:size(ALLaL,1),:);
[RHO,PVAL] = corr(XXX1,XXX2);

figure;
imagesc(RHO);colorbar;
xlabel('channels');ylabel('trials'):title('Correlation coefficent between velocity and spectrogram');

figure;
topoplot(mean(RHO,1)',test,'nosedir','+Y');ColorLim=1;
colormap jet; colorbar ; caxis([-1*ColorLim ColorLim]);
title('Correlation coefficent between velocity and spectrogram');
%% one trial for right hand (skip)
clear RHO PVAL corr
XXX1=ALLaR(:,i);XXX2=theta_trial_spectrogram(1:size(ALLaR,1),:);
[RHO,PVAL] = corr(XXX1,XXX2);

figure;
imagesc(RHO);colorbar;
xlabel('channels');ylabel('trials');title('Correlation coefficent between velocity and spectrogram');

figure;
topoplot(mean(RHO,1)',test,'nosedir','+Y');ColorLim=1;
colormap jet; colorbar ; caxis([-1*ColorLim ColorLim]);
title('Correlation coefficent between velocity and spectrogram');

%% compare left and right acceleration, all good trials (skip)
subplot(1,2,1);imagesc(ALLaL(:,CondiDataGoodTrials));
subplot(1,2,2);imagesc(ALLaR(:,CondiDataGoodTrials));
%% compare left and right acceleration, in one trials (skip)
subplot(1,2,1);plot(ALLaL(:,CondiDataGoodTrials(i)));ylim([-7e-3 7e-3]);
subplot(1,2,2);plot(ALLaR(:,CondiDataGoodTrials(i)));ylim([-7e-3 7e-3]);

%% plot the average (Scatter plot)
clear corr;
figure('units','normalized','outerposition',[0 0 1 1]);
for conditionSelected=1:length(UniCondi);

%     indtemp=find(CondiData==UniCondi(conditionSelected)); % all trials
%     HighInd=indtemp(find(TrialScores(indtemp)>median(TrialScores(indtemp))));
%     LowInd=indtemp(find(TrialScores(indtemp)<median(TrialScores(indtemp))));

    indtemp=find(CondiDataGoodTrials==UniCondi(conditionSelected));
    HighInd=indtemp(find(TrialScoresGoodTrials(indtemp)>median(TrialScoresGoodTrials(indtemp))));
    LowInd=indtemp(find(TrialScoresGoodTrials(indtemp)<median(TrialScoresGoodTrials(indtemp))));
    
    % get the average velocity and acceleration
    averageALLvL=mean(ALLvL(:,indtemp),2);
    averageALLvR=mean(ALLvR(:,indtemp),2);
    % examine
%     figure;plot(averageALLvL,'r');hold on;plot(averageALLvR,'b');
%     plot(average_spectrogram(1:length(averageALLvL),1),'k');

    averageALLaL=mean(ALLaL(:,indtemp),2);
    averageALLaR=mean(ALLaR(:,indtemp),2);
    
    % select a frequency
    wfreq; 
    FreqSelected=3;
    % select trials and average the spectrogram across the trials
    trials_spectrogram=squeeze(normPowcnorm_ALLchan_sgolay_ALLtrials(FreqSelected,:,:,HighInd));
    average_spectrogram=mean(trials_spectrogram,3);
    
    % plot the velocity corr left
    subplot(4,7,7+1-conditionSelected);
    RHO = corr(averageALLvL,average_spectrogram(1:1499,goodchans));
%   figure;plot(averageALLvL,'r');hold on;plot(averageALLvR,'b');plotx(average_spectrogram(1:length(averageALLvL),goodchans),'k');
%   figure;subplot(1,2,1);topoplot(RHO,test(goodchans),'nosedir','+Y','conv','on');ColorLim=1;subplot(1,2,2);
    plot(1:length(goodchans),RHO,'ro');xlabel('channels');ylabel('correlation coefficient');ylim([-1 1]);
%     topoplot(RHO,test(goodchans),'nosedir','+Y');colorbar; caxis([-1 1]);
    title(['condition ' conditionNames{UniCondi(conditionSelected)} ' Lv']);

    % plot the velocity corr right
    subplot(4,7,14+1-conditionSelected);
    RHO = corr(averageALLvR,average_spectrogram(1:length(averageALLvL),goodchans));
    plot(1:length(goodchans),RHO,'ro');xlabel('channels');ylabel('correlation coefficient');ylim([-1 1]);
%     topoplot(RHO,test(goodchans),'nosedir','+Y');colorbar; caxis([-1 1]);
    title(['condition ' conditionNames{UniCondi(conditionSelected)} ' Rv']);
    
    % plot the accelartaion left
    subplot(4,7,21+1-conditionSelected);
    RHO = corr(averageALLaL,average_spectrogram(1:length(averageALLaL),goodchans));
    plot(1:length(goodchans),RHO,'ro');xlabel('channels');ylabel('correlation coefficient');ylim([-1 1]);
%     topoplot(RHO,test(goodchans),'nosedir','+Y');colorbar; caxis([-1 1]);
    title(['condition ' conditionNames{UniCondi(conditionSelected)} ' La']);

    % plot the accelartaion left
    subplot(4,7,28+1-conditionSelected);
    RHO = corr(averageALLaR,average_spectrogram(1:length(averageALLaL),goodchans));
    plot(1:length(goodchans),RHO,'ro');xlabel('channels');ylabel('correlation coefficient');ylim([-1 1]);
%     topoplot(RHO,test(goodchans),'nosedir','+Y');colorbar; caxis([-1 1]);
    title(['condition ' conditionNames{UniCondi(conditionSelected)} ' Ra']);
end
%% plot the average (Scatter topoplot)
clear corr;
figure('units','normalized','outerposition',[0 0 1 1]);
for conditionSelected=1:length(UniCondi);

%     indtemp=find(CondiData==UniCondi(conditionSelected)); % all trials
%     HighInd=indtemp(find(TrialScores(indtemp)>median(TrialScores(indtemp))));
%     LowInd=indtemp(find(TrialScores(indtemp)<median(TrialScores(indtemp))));

    indtemp=find(CondiDataGoodTrials==UniCondi(conditionSelected));
    HighInd=indtemp(find(TrialScoresGoodTrials(indtemp)>median(TrialScoresGoodTrials(indtemp))));
    LowInd=indtemp(find(TrialScoresGoodTrials(indtemp)<median(TrialScoresGoodTrials(indtemp))));
    
    % get the average velocity and acceleration
    averageALLvL=mean(ALLvL(:,indtemp),2);
    averageALLvR=mean(ALLvR(:,indtemp),2);
    % examine
%     figure;plot(averageALLvL,'r');hold on;plot(averageALLvR,'b');
%     plot(average_spectrogram(1:length(averageALLvL),1),'k');

    averageALLaL=mean(ALLaL(:,indtemp),2);
    averageALLaR=mean(ALLaR(:,indtemp),2);
    
    % select a frequency
    wfreq; 
    FreqSelected=3;
    % select trials and average the spectrogram across the trials
    trials_spectrogram=squeeze(normPowcnorm_ALLchan_sgolay_ALLtrials(FreqSelected,:,:,HighInd));
    average_spectrogram=mean(trials_spectrogram,3);
    
    % plot the velocity corr left
    subplot(4,7,7+1-conditionSelected);
    RHO = corr(averageALLvL,average_spectrogram(1:1499,goodchans));
%   figure;plot(averageALLvL,'r');hold on;plot(averageALLvR,'b');plotx(average_spectrogram(1:length(averageALLvL),goodchans),'k');
%   figure;subplot(1,2,1);topoplot(RHO,test(goodchans),'nosedir','+Y','conv','on');ColorLim=1;subplot(1,2,2);
%     plot(1:length(goodchans),RHO,'ro');xlabel('channels');ylabel('correlation coefficient');ylim([-1 1]);
    topoplot(RHO,test(goodchans),'nosedir','+Y');colorbar; caxis([-1 1]);
    title(['condition ' conditionNames{UniCondi(conditionSelected)} ' Lv']);

    % plot the velocity corr right
    subplot(4,7,14+1-conditionSelected);
    RHO = corr(averageALLvR,average_spectrogram(1:length(averageALLvL),goodchans));
%     plot(1:length(goodchans),RHO,'ro');xlabel('channels');ylabel('correlation coefficient');ylim([-1 1]);
    topoplot(RHO,test(goodchans),'nosedir','+Y');colorbar; caxis([-1 1]);
    title(['condition ' conditionNames{UniCondi(conditionSelected)} ' Rv']);
    
    % plot the accelartaion left
    subplot(4,7,21+1-conditionSelected);
    RHO = corr(averageALLaL,average_spectrogram(1:length(averageALLaL),goodchans));
%     plot(1:length(goodchans),RHO,'ro');xlabel('channels');ylabel('correlation coefficient');ylim([-1 1]);
    topoplot(RHO,test(goodchans),'nosedir','+Y');colorbar; caxis([-1 1]);
    title(['condition ' conditionNames{UniCondi(conditionSelected)} ' La']);

    % plot the accelartaion left
    subplot(4,7,28+1-conditionSelected);
    RHO = corr(averageALLaR,average_spectrogram(1:length(averageALLaL),goodchans));
%     plot(1:length(goodchans),RHO,'ro');xlabel('channels');ylabel('correlation coefficient');ylim([-1 1]);
    topoplot(RHO,test(goodchans),'nosedir','+Y');colorbar; caxis([-1 1]);
    title(['condition ' conditionNames{UniCondi(conditionSelected)} ' Ra']);
end
%% PLS Regression (skip)
% https://www.mathworks.com/help/stats/plsregress.html#mw_635e159a-1335-4e52-af7f-5fc82592798a
[XLOADINGS,YLOADINGS] = plsregress(average_spectrogram(1:length(averageALLvL),goodchans),averageALLvL);
%% PLS Regression tryout (matlab example)  (skip)
[XL,yl,XS,YS,beta,PCTVAR,MSE,stats] = plsregress(average_spectrogram(1:length(averageALLvL),goodchans),averageALLvL,119,'cv',5,'Options',statset('UseParallel',true));
[XL,yl,XS,YS,beta,PCTVAR,MSE,stats] = plsregress(average_spectrogram(1:length(averageALLvL),goodchans),averageALLvL,119,'cv',5);
[XL,yl,XS,YS,beta,PCTVAR,MSE,stats] = plsregress(average_spectrogram(1:length(averageALLvL),goodchans),averageALLvL,10,'cv',5);
W0 = stats.W ./ sqrt(sum(stats.W.^2,1));
p = size(XL,1);
sumSq = sum(XS.^2,1).*sum(yl.^2,1);
vipScore = sqrt(p* sum(sumSq.*(W0.^2),2) ./ sum(sumSq,2));
indVIP = find(vipScore >= 1);
scatter(1:length(vipScore),vipScore,'x')
hold on
scatter(indVIP,vipScore(indVIP),'rx')
plot([1 length(vipScore)],[1 1],'--k')
hold off
axis tight
xlabel('Channels')
ylabel('variable importance in projection (VIP) scores')
%% PLS Regression (5 fold)
% https://www.mathworks.com/help/stats/plsregress.html
AllPredictChans=[];
figure('units','normalized','outerposition',[0 0 1 1]);
for conditionSelected=1:length(UniCondi);

%     indtemp=find(CondiData==UniCondi(conditionSelected)); % all trials
%     HighInd=indtemp(find(TrialScores(indtemp)>median(TrialScores(indtemp))));
%     LowInd=indtemp(find(TrialScores(indtemp)<median(TrialScores(indtemp))));

    indtemp=find(CondiDataGoodTrials==UniCondi(conditionSelected));
    HighInd=indtemp(find(TrialScoresGoodTrials(indtemp)>median(TrialScoresGoodTrials(indtemp))));
    LowInd=indtemp(find(TrialScoresGoodTrials(indtemp)<median(TrialScoresGoodTrials(indtemp))));
    
    % get the average velocity and acceleration
    averageALLvL=mean(ALLvL(:,indtemp),2);
    averageALLvR=mean(ALLvR(:,indtemp),2);
    % examine
%     figure;plot(averageALLvL,'r');hold on;plot(averageALLvR,'b');
%     plot(average_spectrogram(1:length(averageALLvL),goodchans),'k');

    averageALLaL=mean(ALLaL(:,indtemp),2);
    averageALLaR=mean(ALLaR(:,indtemp),2);
    
    % select a frequency
    wfreq; 
    FreqSelected=3;
    % select trials and average the spectrogram across the trials
    trials_spectrogram=squeeze(normPowcnorm_ALLchan_sgolay_ALLtrials(FreqSelected,:,:,HighInd));
    average_spectrogram=mean(trials_spectrogram,3);
    
    % plot the velocity corr left
    subplot(4,7,7+1-conditionSelected);
    [XL,yl,XS,YS,beta,PCTVAR,MSE,stats] = plsregress(average_spectrogram(1:length(averageALLvL),goodchans),averageALLvL,10,'cv',5);
    W0 = stats.W ./ sqrt(sum(stats.W.^2,1));
    p = size(XL,1);
    sumSq = sum(XS.^2,1).*sum(yl.^2,1);
    vipScore = sqrt(p* sum(sumSq.*(W0.^2),2) ./ sum(sumSq,2));
    indVIP = find(vipScore >= 1);
    
    scatter(1:length(vipScore),vipScore,'x')
    hold on
    scatter(indVIP,vipScore(indVIP),'rx')
    plot([1 length(vipScore)],[1 1],'--k')
    hold off
    axis tight
    xlabel('Channels')
    ylabel('VIP scores')
    title(['condition ' conditionNames{UniCondi(conditionSelected)} ' Lv']);
    ylim([-1 5]);
    
    vipScorelabels =  strsplit(num2str((1:length(vipScore)')));
    text(1:length(vipScore),vipScore,vipScorelabels,'VerticalAlignment','top','HorizontalAlignment','right');
    MaxGoodChan=find(vipScore==max(vipScore));
    vipScorelabels = {num2str(MaxGoodChan)};
    text(MaxGoodChan,max(vipScore),vipScorelabels,'VerticalAlignment','top','HorizontalAlignment','right','color',[1 0 0]);
    
    [bigvalues, bigidx] = sort(vipScore, 'descend');
    PredictChans=bigidx(1:5);
    xlabel(['PredictChans ' num2str(PredictChans')]);
    AllPredictChans=[AllPredictChans;PredictChans];
    
    % plot the velocity corr right
    subplot(4,7,14+1-conditionSelected);
    [XL,yl,XS,YS,beta,PCTVAR,MSE,stats] = plsregress(average_spectrogram(1:length(averageALLvR),goodchans),averageALLvR,10,'cv',5);
    W0 = stats.W ./ sqrt(sum(stats.W.^2,1));
    p = size(XL,1);
    sumSq = sum(XS.^2,1).*sum(yl.^2,1);
    vipScore = sqrt(p* sum(sumSq.*(W0.^2),2) ./ sum(sumSq,2));
    indVIP = find(vipScore >= 1);
    scatter(1:length(vipScore),vipScore,'x')
    hold on
    scatter(indVIP,vipScore(indVIP),'rx')
    plot([1 length(vipScore)],[1 1],'--k')
    hold off
    axis tight
    xlabel('Channels')
    ylabel('VIP scores')
    title(['condition ' conditionNames{UniCondi(conditionSelected)} ' Rv']);
    ylim([-1 5]);
    
    vipScorelabels =  strsplit(num2str((1:length(vipScore)')));
    text(1:length(vipScore),vipScore,vipScorelabels,'VerticalAlignment','top','HorizontalAlignment','right');
    MaxGoodChan=find(vipScore==max(vipScore));
    vipScorelabels = {num2str(MaxGoodChan)};
    text(MaxGoodChan,max(vipScore),vipScorelabels,'VerticalAlignment','top','HorizontalAlignment','right','color',[1 0 0]);
    
    [bigvalues, bigidx] = sort(vipScore, 'descend');
    PredictChans=bigidx(1:5);
    xlabel(['PredictChans ' num2str(PredictChans')]);
    AllPredictChans=[AllPredictChans;PredictChans];
    
    % plot the accelartaion left
    subplot(4,7,21+1-conditionSelected);
    [XL,yl,XS,YS,beta,PCTVAR,MSE,stats] = plsregress(average_spectrogram(1:length(averageALLaL),goodchans),averageALLaL,10,'cv',5);
    W0 = stats.W ./ sqrt(sum(stats.W.^2,1));
    p = size(XL,1);
    sumSq = sum(XS.^2,1).*sum(yl.^2,1);
    vipScore = sqrt(p* sum(sumSq.*(W0.^2),2) ./ sum(sumSq,2));
    indVIP = find(vipScore >= 1);
    scatter(1:length(vipScore),vipScore,'x')
    hold on
    scatter(indVIP,vipScore(indVIP),'rx')
    plot([1 length(vipScore)],[1 1],'--k')
    hold off
    axis tight
    xlabel('Channels')
    ylabel('VIP scores')
    title(['condition ' conditionNames{UniCondi(conditionSelected)} ' La']);
    ylim([-1 5]);
    
    vipScorelabels =  strsplit(num2str((1:length(vipScore)')));
    text(1:length(vipScore),vipScore,vipScorelabels,'VerticalAlignment','top','HorizontalAlignment','right');
    MaxGoodChan=find(vipScore==max(vipScore));
    vipScorelabels = {num2str(MaxGoodChan)};
    text(MaxGoodChan,max(vipScore),vipScorelabels,'VerticalAlignment','top','HorizontalAlignment','right','color',[1 0 0]);
    
    [bigvalues, bigidx] = sort(vipScore, 'descend');
    PredictChans=bigidx(1:5);
    xlabel(['PredictChans ' num2str(PredictChans')]);
    AllPredictChans=[AllPredictChans;PredictChans];

    % plot the accelartaion left
    subplot(4,7,28+1-conditionSelected);
    [XL,yl,XS,YS,beta,PCTVAR,MSE,stats] = plsregress(average_spectrogram(1:length(averageALLaR),goodchans),averageALLaR,10,'cv',5);
    W0 = stats.W ./ sqrt(sum(stats.W.^2,1));
    p = size(XL,1);
    sumSq = sum(XS.^2,1).*sum(yl.^2,1);
    vipScore = sqrt(p* sum(sumSq.*(W0.^2),2) ./ sum(sumSq,2));
    indVIP = find(vipScore >= 1);
    scatter(1:length(vipScore),vipScore,'x')
    hold on
    scatter(indVIP,vipScore(indVIP),'rx')
    plot([1 length(vipScore)],[1 1],'--k')
    hold off
    axis tight
    xlabel('Channels')
    ylabel('VIP scores')
    title(['condition ' conditionNames{UniCondi(conditionSelected)} ' Ra']);
    ylim([-1 5]);
    
    vipScorelabels =  strsplit(num2str((1:length(vipScore)')));
    text(1:length(vipScore),vipScore,vipScorelabels,'VerticalAlignment','top','HorizontalAlignment','right');
    MaxGoodChan=find(vipScore==max(vipScore));
    vipScorelabels = {num2str(MaxGoodChan)};
    text(MaxGoodChan,max(vipScore),vipScorelabels,'VerticalAlignment','top','HorizontalAlignment','right','color',[1 0 0]);
    
    [bigvalues, bigidx] = sort(vipScore, 'descend');
    PredictChans=bigidx(1:5);
    xlabel(['PredictChans ' num2str(PredictChans')]);
    AllPredictChans=[AllPredictChans;PredictChans];
end

AllPredictChans;
AllPredictChans=unique(AllPredictChans);
%% Look for what channels they are
labels;
PredictChansToMark=goodchans([67]);
labels(goodchans([61]))
labels(goodchans([61 97 27]))  % {'75-AFF2H'}    {'113-CPP4H'}    {'32-C4'}
labels(goodchans(AllPredictChans'))
labels(goodchans([13 49 67 63 75 89]))
labels(PredictChansToMark)


%% Calculate subplot coordinates for scalp map
cd /home/zhibin/Documents/GitHub/Motor_cordination/EEGanalysis
load('coordinates.mat');
X = [chanlocs128.X];
Y = [chanlocs128.Y];
Z = [chanlocs128.Z];
% labels = {'FP1','FPZ','FP2','AF3','AF4','F11','F7','F5','F3','F1','FZ','F2','F4','F6','F8','F12','FT11','FC5','FC3','FC1','FCZ','FC2','FC4','FC6','FT12','T7','C5','C3','C1','CZ','C2','C4','C6','T8','TP7','CP5','CP3','CP1','CPZ','CP2','CP4','CP6','TP8','M1','M2','P7','P5','P3','P1','PZ','P2','P4','P6','P8','PO7','PO3','POZ','PO4','PO8','O1','OZ','O2','CB1','CB2','AFP1','AFP2','AF7','AF5','AFZ','AF6','AF8','AFF5H','AFF3H','AFF1H','AFF2H','AFF4H','AFF6H','F9','F10','FFT7H','FFC5H','FFC3H','FFC1H','FFC2H','FFC4H','FFC6H','FFT8H','FT9','FT7','FT8','FT10','FTT7H','FCC5H','FCC3H','FCC1H','FCC2H','FCC4H','FCC6H','FTT8H','TTP7H','CCP5H','CCP3H','CCP1H','CCP2H','CCP4H','CCP6H','TTP8H','TPP7H','CPP5H','CPP3H','CPP1H','CPP2H','CPP4H','CPP6H','TPP8H','P9','P10','PPO3H','PPO1H','PPO2H','PPO4H','PO9','PO5','PO1','PO2','PO6','PO10','CBZ'};
labels = {'1-FP1', '2-FPZ', '3-FP2', '4-AF3', '5-AF4', '6-F11', '7-F7', '8-F5', '9-F3', '10-F1', '11-FZ', '12-F2', '13-F4', '14-F6', '15-F8', '16-F12', '17-FT11', '18-FC5', '19-FC3', '20-FC1', '21-FCZ', '22-FC2', '23-FC4', '24-FC6', '25-FT12', '26-T7', '27-C5', '28-C3', '29-C1', '30-CZ', '31-C2', '32-C4', '33-C6', '34-T8', '35-TP7', '36-CP5', '37-CP3', '38-CP1', '39-CPZ', '40-CP2', '41-CP4', '42-CP6', '43-TP8', '44-M1', '45-M2', '46-P7', '47-P5', '48-P3', '49-P1', '50-PZ', '51-P2', '52-P4', '53-P6', '54-P8', '55-PO7', '56-PO3', '57-POZ', '58-PO4', '59-PO8', '60-O1', '61-OZ', '62-O2', '63-CB1', '64-CB2', '65-AFP1', '66-AFP2', '67-AF7', '68-AF5', '69-AFZ', '70-AF6', '71-AF8', '72-AFF5H', '73-AFF3H', '74-AFF1H', '75-AFF2H', '76-AFF4H', '77-AFF6H', '78-F9', '79-F10', '80-FFT7H', '81-FFC5H', '82-FFC3H', '83-FFC1H', '84-FFC2H', '85-FFC4H', '86-FFC6H', '87-FFT8H', '88-FT9', '89-FT7', '90-FT8', '91-FT10', '92-FTT7H', '93-FCC5H', '94-FCC3H', '95-FCC1H', '96-FCC2H', '97-FCC4H', '98-FCC6H', '99-FTT8H', '100-TTP7H', '101-CCP5H', '102-CCP3H', '103-CCP1H', '104-CCP2H', '105-CCP4H', '106-CCP6H', '107-TTP8H', '108-TPP7H', '109-CPP5H', '110-CPP3H', '111-CPP1H', '112-CPP2H', '113-CPP4H', '114-CPP6H', '115-TPP8H', '116-P9', '117-P10', '118-PPO3H', '119-PPO1H', '120-PPO2H', '121-PPO4H', '122-PO9', '123-PO5', '124-PO1', '125-PO2', '126-PO6', '127-PO10', '128-CBZ'};

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
%% plot scalp map and color predictive chans red on the scalp map

figure('units','normalized','outerposition',[0 0 1 1]);
plot(XXPLOT,YYPLOT,'r.');

% labels = {'FP1','FPZ','FP2','AF3','AF4','F11','F7','F5','F3','F1','FZ','F2','F4','F6','F8','F12','FT11','FC5','FC3','FC1','FCZ','FC2','FC4','FC6','FT12','T7','C5','C3','C1','CZ','C2','C4','C6','T8','TP7','CP5','CP3','CP1','CPZ','CP2','CP4','CP6','TP8','M1','M2','P7','P5','P3','P1','PZ','P2','P4','P6','P8','PO7','PO3','POZ','PO4','PO8','O1','OZ','O2','CB1','CB2','AFP1','AFP2','AF7','AF5','AFZ','AF6','AF8','AFF5H','AFF3H','AFF1H','AFF2H','AFF4H','AFF6H','F9','F10','FFT7H','FFC5H','FFC3H','FFC1H','FFC2H','FFC4H','FFC6H','FFT8H','FT9','FT7','FT8','FT10','FTT7H','FCC5H','FCC3H','FCC1H','FCC2H','FCC4H','FCC6H','FTT8H','TTP7H','CCP5H','CCP3H','CCP1H','CCP2H','CCP4H','CCP6H','TTP8H','TPP7H','CPP5H','CPP3H','CPP1H','CPP2H','CPP4H','CPP6H','TPP8H','P9','P10','PPO3H','PPO1H','PPO2H','PPO4H','PO9','PO5','PO1','PO2','PO6','PO10','CBZ'};;
text(XXPLOT,YYPLOT,labels,'VerticalAlignment','bottom','HorizontalAlignment','right');

% color it red on the scalp map
PredictChanlabels = labels(PredictChansToMark);
text(XXPLOT(PredictChansToMark),YYPLOT(PredictChansToMark),PredictChanlabels,'VerticalAlignment','bottom','HorizontalAlignment','right','color',[1 0 0]);
    