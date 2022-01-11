%% (skip)
cd /home/zhibin/Documents/GitHub/Motor_cordination/EEGanalysis/20211102
clear data_trials
run integrate_EEG_into_data_trials_step3.m % now included the baseline and padding
data_trials % check whether included the baseline
%% wavelet
% wfreq = 1:50; % frequency(s) of interest
wfreq = [2 4 6 8 10 14 18 24 30 40];
wfc = 1.5;
sr=1000;

UniCondi=unique(CondiData);

% plot out pow the 128 channels
AllchanNames={'FP1','FPZ','FP2','AF3','AF4','F11','F7','F5','F3','F1','FZ','F2','F4','F6','F8','F12','FT11','FC5','FC3','FC1','FCZ','FC2','FC4','FC6','FT12','T7','C5','C3','C1','CZ','C2','C4','C6','T8','TP7','CP5','CP3','CP1','CPZ','CP2','CP4','CP6','TP8','M1','M2','P7','P5','P3','P1','PZ','P2','P4','P6','P8','PO7','PO3','POZ','PO4','PO8','O1','OZ','O2','CB1','CB2','AFP1','AFP2','AF7','AF5','AFZ','AF6','AF8','AFF5H','AFF3H','AFF1H','AFF2H','AFF4H','AFF6H','F9','F10','FFT7H','FFC5H','FFC3H','FFC1H','FFC2H','FFC4H','FFC6H','FFT8H','FT9','FT7','FT8','FT10','FTT7H','FCC5H','FCC3H','FCC1H','FCC2H','FCC4H','FCC6H','FTT8H','TTP7H','CCP5H','CCP3H','CCP1H','CCP2H','CCP4H','CCP6H','TTP8H','TPP7H','CPP5H','CPP3H','CPP1H','CPP2H','CPP4H','CPP6H','TPP8H','P9','P10','PPO3H','PPO1H','PPO2H','PPO4H','PO9','PO5','PO1','PO2','PO6','PO10','CBZ','VEOG','HEOG','EMG1','EMG2','HL 1','HL 2','EMG3','EMG4','EMG5','EMG6','TRIGGER'};

% AllchanNames{1}

%% For all 3 days (skip)
CondiData=allPerm_alldays(logical(goodepochs_alldays));
TrialScores=TrialScores_alldays(logical(goodepochs_alldays));
%% all conditions
for u=1:length(UniCondi);

    figure('units','normalized','outerposition',[0 0 0.6 0.6]);

%     indtemp=find(CondiData==UniCondi(u));
%     HighInd=indtemp(find(TrialScores(indtemp)>median(TrialScores(indtemp))));
%     LowInd=indtemp(find(TrialScores(indtemp)<median(TrialScores(indtemp))));
    indtemp=find(CondiDataGoodTrials==UniCondi(u));
    HighInd=indtemp(find(TrialScoresGoodTrials(indtemp)>median(TrialScoresGoodTrials(indtemp))));
    LowInd=indtemp(find(TrialScoresGoodTrials(indtemp)<median(TrialScoresGoodTrials(indtemp))));
    
    % Compute ERP in all channels
    basedlinecorrected_ERP=mean(baselinecorrected_trial(:,1:128,HighInd),3);
    % plot(basedlinecorrected_ERP);

    cnorm = wavelet(basedlinecorrected_ERP,sr,wfc,wfreq);
    Powcnorm = abs(cnorm).^2;
    
    % Power normalization
    % one way to it:
    logPowcorm=log10(Powcnorm);
    baselineMean=log10(mean(Powcnorm(:,1:500,:),2));
    normPowcnorm = logPowcorm-(ones(1,size(logPowcorm,2),1).*baselineMean);
    % second way to do it:
%     baselineMean = mean(Powcnorm(:,501:1000,:),2);
%     baselineMeanArray = ones(1,4000,1).*baselineMean;
%     normPowcnorm = (Powcnorm-baselineMeanArray)./(baselineMeanArray);

    
    % open sgolay; open sgolayfilt
    ColorLim=2;

    subplot(3,2,1);
    imagesc((sgolayfilt(squeeze(normPowcnorm(:,:,9))',1,31))'); % F3/F4 9/13  
    colormap jet; colorbar ; caxis([-1*ColorLim ColorLim]);
    xlabel('time(ms)'); ylabel('frequencies(Hz)'); set(gca,'ydir','normal');
    yticks([1:length(wfreq)]);yticklabels({'2','4','6','8','10','14','18','24','30','40'});
    xticks(linspace(0,4000,9));xticklabels({'-1000','-500','0','500','1000','1500','2000','2500','3000'});
    title([AllchanNames{9}]);
    hold on;xline(1000,'g','linewidth',5);xline(2500,'r','linewidth',5);hold off;
    
    subplot(3,2,2);
    imagesc((sgolayfilt(squeeze(normPowcnorm(:,:,13))',1,31))'); % F3/F4 9/13  
    colormap jet; colorbar ; caxis([-1*ColorLim ColorLim]);
    xlabel('time(ms)'); ylabel('frequencies(Hz)'); set(gca,'ydir','normal');
    yticks([1:length(wfreq)]);yticklabels({'2','4','6','8','10','14','18','24','30','40'});
    xticks(linspace(0,4000,9));xticklabels({'-1000','-500','0','500','1000','1500','2000','2500','3000'});
    title([AllchanNames{13}]);
    hold on;xline(1000,'g','linewidth',5);xline(2500,'r','linewidth',5);hold off;

    subplot(3,2,3);
    imagesc((sgolayfilt(squeeze(normPowcnorm(:,:,19))',1,31))'); % FC3/FC4 19/23 
    colormap jet; colorbar ; caxis([-1*ColorLim ColorLim]);
    xlabel('time(ms)'); ylabel('frequencies(Hz)'); set(gca,'ydir','normal');
    yticks([1:length(wfreq)]);yticklabels({'2','4','6','8','10','14','18','24','30','40'});
    xticks(linspace(0,4000,9));xticklabels({'-1000','-500','0','500','1000','1500','2000','2500','3000'});
    title([AllchanNames{19}]);
    hold on;xline(1000,'g','linewidth',5);xline(2500,'r','linewidth',5);hold off;

    subplot(3,2,4);
    imagesc((sgolayfilt(squeeze(normPowcnorm(:,:,23))',1,31))'); % FC3/FC4 19/23 
    colormap jet; colorbar ; caxis([-1*ColorLim ColorLim]);
    xlabel('time(ms)'); ylabel('frequencies(Hz)'); set(gca,'ydir','normal');
    yticks([1:length(wfreq)]);yticklabels({'2','4','6','8','10','14','18','24','30','40'});
    xticks(linspace(0,4000,9));xticklabels({'-1000','-500','0','500','1000','1500','2000','2500','3000'});
    title([AllchanNames{23}]);
    hold on;xline(1000,'g','linewidth',5);xline(2500,'r','linewidth',5);hold off;

    subplot(3,2,5);
    imagesc((sgolayfilt(squeeze(normPowcnorm(:,:,28))',1,31))'); % C3 
    colormap jet; colorbar ; caxis([-1*ColorLim ColorLim]);
    xlabel('time(ms)'); ylabel('frequencies(Hz)'); set(gca,'ydir','normal');
    yticks([1:length(wfreq)]);yticklabels({'2','4','6','8','10','14','18','24','30','40'});
    xticks(linspace(0,4000,9));xticklabels({'-1000','-500','0','500','1000','1500','2000','2500','3000'});
    title([AllchanNames{28}]);
    hold on;xline(1000,'g','linewidth',5);xline(2500,'r','linewidth',5);hold off;

    subplot(3,2,6);
    imagesc((sgolayfilt(squeeze(normPowcnorm(:,:,32))',1,31))'); % C4 
    colormap jet; colorbar ; caxis([-1*ColorLim ColorLim]);
    xlabel('time(ms)'); ylabel('frequencies(Hz)'); set(gca,'ydir','normal');
    yticks([1:length(wfreq)]);yticklabels({'2','4','6','8','10','14','18','24','30','40'});
    xticks(linspace(0,4000,9));xticklabels({'-1000','-500','0','500','1000','1500','2000','2500','3000'});
    title([AllchanNames{32}]);
    hold on;xline(1000,'g','linewidth',5);xline(2500,'r','linewidth',5);hold off;
    
    suptitle([conditionNames{UniCondi(u)}]);

end

%% Only condition 4:0 and 0:4 but all channels on scalp map
conditionNames={'0:4' '1:4' '1:2' '1:1' '2:1' '4:1' '4:0'}; 
win=1001:2000; 

% all performance trials
for u=[1 7]; % sp=1:2
    
    figure('units','normalized','outerposition',[0 0 1 1]);
    indtemp=find(CondiDataGoodTrials==UniCondi(u));
    
    % Compute ERP in all channels
    basedlinecorrected_ERP=mean(baselinecorrected_trial(:,1:128,indtemp),3);
    % plot(basedlinecorrected_ERP);

    cnorm = wavelet(basedlinecorrected_ERP,sr,wfc,wfreq);
    Powcnorm = abs(cnorm).^2;
    
    % Power normalization
    % one way to it:
    logPowcorm=log10(Powcnorm);
    baselineMean=log10(mean(Powcnorm(:,1:500,:),2));
    normPowcnorm = logPowcorm-(ones(1,size(logPowcorm,2),1).*baselineMean);
    
    % open sgolay; open sgolayfilt
    ColorLim=2;
    
    for chan=goodchans
        subplot('Position',[XXPLOT(chan) YYPLOT(chan) 0.02 0.03]);
        imagesc((sgolayfilt(squeeze(normPowcnorm(:,:,chan))',1,31))');
        colormap jet; caxis([-1*ColorLim ColorLim]);
        set(gca,'ydir','normal');
        title([AllchanNames{chan}]);
        hold on;xline(500,'k','linewidth',1);xline(1000,'k','linewidth',1);hold off;
        if ~isempty(find([1:127]==chan))
        set(gca,'XTick',[]); set(gca,'YTick',[]); 
        end
        if chan==128
            subplot('Position',[XXPLOT(chan)-0.04 YYPLOT(chan)-0.02 0.1 0.08]);
            imagesc((sgolayfilt(squeeze(normPowcnorm(:,:,chan))',1,31))');
            colormap jet; caxis([-1*ColorLim ColorLim]);
            set(gca,'ydir','normal');
            title([AllchanNames{chan}]);
            hold on;xline(500,'k','linewidth',1);xline(1000,'k','linewidth',1);hold off;
            xlabel('time(ms)'); ylabel('frequencies(Hz)'); 
            yticks([1:length(wfreq)]);yticklabels({'2','4','6','8','10','14','18','24','30','40'});
            xticks(linspace(0,2000,5));xticklabels({'-1000','-500','0','500','1000'});
            colorbar;
        end
    end
    suptitle(conditionNames(u));
end



%% https://github.com/rameshsrinivasanuci/matlab/blob/master/jenny/WaveletTransform.m (skip)
% now let's recover the time course
wfreq = 1:50; % frequency(s) of interest
wfc = 3;
cnorm = wavelet(erp,sr,wfc,wfreq);

timecourse = complex(zeros(length(wfreq),4000));
for i = 1:4000
    timepoint = squeeze(cnorm(:,i,:));
    timesignal = timepoint * weights;
    timecourse(:,i) = timesignal;
end

% smooth it, fit a line to smooth the spectrum 
sgftimecourse = sgolayfilt(abs(timecourse.'),1,31);


% now let's loop through individual trial
trialtimecourse = complex(zeros(4000,length(wfreq),360));
for i = 1:length(finalgoodtrials)
        trialdata = data(:,:,finalgoodtrials(i));
        trialwavelet = wavelet(trialdata,sr,wfc,wfreq);
    for j = 1:length(wfreq)
        trialfreq = squeeze(trialwavelet(j,:,:));
        trialsignal = trialfreq * weights;
        trialtimecourse(:,j,i) = trialsignal;
    end
end

% extract the power at the target frequency


% making plots
figure;imagesc(abs(sgftimecourse)')
colormap jet
colorbar 
caxis ([0 2])
xlabel('time(s)');
ylabel('frequencies(Hz)');

figure;imagesc(abs(timecourse))
colormap jet
colorbar 
caxis 
xlabel('time(s)');
ylabel('frequencies(Hz)');


% smooth it, fit a line to smooth the spectrum 
sgftimecourse = sgolayfilt(abs(timecourse.'),1,31);
