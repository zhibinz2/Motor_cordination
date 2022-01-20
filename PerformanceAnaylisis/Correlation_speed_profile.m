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

% select a frequency
wfreq; 
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

%% compare speed profile of left and right hand in all good trials
subplot(1,2,1);imagesc(ALLvL(:,CondiDataGoodTrials));
subplot(1,2,2);imagesc(ALLvR(:,CondiDataGoodTrials));
%% compare speed profile of left and right hand in one trial
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

%% compare left and right acceleration, all good trials
subplot(1,2,1);imagesc(ALLaL(:,CondiDataGoodTrials));
subplot(1,2,2);imagesc(ALLaR(:,CondiDataGoodTrials));
%% compare left and right acceleration, in one trials
subplot(1,2,1);plot(ALLaL(:,CondiDataGoodTrials(i)));ylim([-7e-3 7e-3]);
subplot(1,2,2);plot(ALLaR(:,CondiDataGoodTrials(i)));ylim([-7e-3 7e-3]);

%% plot the average
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
    plot(1:112,RHO,'ro');xlabel('channels');ylabel('correlation coefficient');ylim([-1 1]);
%     topoplot(RHO,test(goodchans),'nosedir','+Y');colorbar; caxis([-1 1]);
    title(['condition ' conditionNames{UniCondi(conditionSelected)} ' Lv']);

    % plot the velocity corr right
    subplot(4,7,14+1-conditionSelected);
    RHO = corr(averageALLvR,average_spectrogram(1:length(averageALLvL),goodchans));
    plot(1:112,RHO,'ro');xlabel('channels');ylabel('correlation coefficient');ylim([-1 1]);
%     topoplot(RHO,test(goodchans),'nosedir','+Y');colorbar; caxis([-1 1]);
    title(['condition ' conditionNames{UniCondi(conditionSelected)} ' Rv']);
    
    % plot the accelartaion left
    subplot(4,7,21+1-conditionSelected);
    RHO = corr(averageALLaL,average_spectrogram(1:length(averageALLaL),goodchans));
    plot(1:112,RHO,'ro');xlabel('channels');ylabel('correlation coefficient');ylim([-1 1]);
%     topoplot(RHO,test(goodchans),'nosedir','+Y');colorbar; caxis([-1 1]);
    title(['condition ' conditionNames{UniCondi(conditionSelected)} ' La']);

    % plot the accelartaion left
    subplot(4,7,28+1-conditionSelected);
    RHO = corr(averageALLaR,average_spectrogram(1:length(averageALLaL),goodchans));
    plot(1:112,RHO,'ro');xlabel('channels');ylabel('correlation coefficient');ylim([-1 1]);
%     topoplot(RHO,test(goodchans),'nosedir','+Y');colorbar; caxis([-1 1]);
    title(['condition ' conditionNames{UniCondi(conditionSelected)} ' Ra']);
end


%% PLS Regression
% https://www.mathworks.com/help/stats/plsregress.html#mw_635e159a-1335-4e52-af7f-5fc82592798a
[XLOADINGS,YLOADINGS] = plsregress(average_spectrogram(1:length(averageALLvL),goodchans),averageALLvL);
%% PLS Regression tryout (matlab example)
[XL,yl,XS,YS,beta,PCTVAR,MSE,stats] = plsregress(average_spectrogram(1:length(averageALLvL),goodchans),averageALLvL);
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
    [XL,yl,XS,YS,beta,PCTVAR,MSE,stats] = plsregress(average_spectrogram(1:length(averageALLvL),goodchans),averageALLvL);
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

    % plot the velocity corr right
    subplot(4,7,14+1-conditionSelected);
    [XL,yl,XS,YS,beta,PCTVAR,MSE,stats] = plsregress(average_spectrogram(1:length(averageALLvR),goodchans),averageALLvR);
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
    
    % plot the accelartaion left
    subplot(4,7,21+1-conditionSelected);
    [XL,yl,XS,YS,beta,PCTVAR,MSE,stats] = plsregress(average_spectrogram(1:length(averageALLaL),goodchans),averageALLaL);
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

    % plot the accelartaion left
    subplot(4,7,28+1-conditionSelected);
    [XL,yl,XS,YS,beta,PCTVAR,MSE,stats] = plsregress(average_spectrogram(1:length(averageALLaR),goodchans),averageALLaR);
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
end

%% Look for what channels they are
labels(goodchans([61]))
labels(goodchans([61 97 27]))  % {'75-AFF2H'}    {'113-CPP4H'}    {'32-C4'}