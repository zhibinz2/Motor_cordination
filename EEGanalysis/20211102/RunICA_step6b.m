% Run ICA on reRef_data
% this should be step 6b

%% Recombine all trials into one big long time series
% % Slow Method
% EEGdataShort=reRef_data(:,:,1);
% for i=2:size(reRef_data,3)
%     %EEGdataShortTemp=reRef_data(:,:,i);
%     EEGdataShort=cat(1,EEGdataShort,reRef_data(:,:,i));
% end

% Fast Method
split_reRef_data = num2cell(reRef_data(:,:,goodepochs), [1 2]); %split A keeping dimension 1 and 2 intact
EEGdataShort=vertcat(split_reRef_data{:});

%% Recalculate the event array (for display and for EEGLAB)
% So the 4000 timepoints are structured as 
% "500ms padding + 500ms green dot baseline + 500ms plan phase + 1000ms movement + 1000ms show bonus + 500ms padding"

NumGoodEpochs=length(goodepochs); % NumTrials=350; or numtotal or NumTrialsRecorded

% Index for start of trials
% IndEvents=round(linspace(1000,4000*NumTrials,NumTrials));
IndStart=[1:2000:(2000*(NumGoodEpochs-1)+1)];
IndEnds=IndStart+2000-1;

%% Index for EEGLAB
% IndEventsEEGLAB=zeros(1,size(reRef_data,1));
% IndEventsEEGLAB(IndStart)=1;
% IndEventsEEGLAB(IndEnds)=2;
% plot(IndEventsEEGLAB,'ro');

%% Examine data before ICA
TimedataShort=1:size(EEGdataShort,1);
plot(TimedataShort, EEGdataShort(:,goodchans));
hold on;
for i=1:length(IndStart)
    xline(IndStart(i),'g-',{num2str(i)});
%     xline(IndEnds(i),'r');
end
hold off;

% trial 41 105 107 are very bad even only good channels were display
figure
DurationPlot=IndStart(40):IndStart(42);
plotx(TimedataShort(DurationPlot),EEGdataShort(DurationPlot,goodchans));hold on;
for i=40:42
    xline(IndStart(i),'g-',{num2str(i)});
end
hold off;

figure;
DurationPlot=IndStart(103):IndStart(108);
plotx(TimedataShort(DurationPlot),EEGdataShort(DurationPlot,goodchans));hold on;
for i=103:108
    xline(IndStart(i),'g-',{num2str(i)});
end
hold off;

%% Run ICA
tic
[icasig, A, W] = fastica(EEGdataShort');
toc

%% Plot ICA component
for i=1:size(A,2)
    SqA(i)=sumsqr(A(:,i));
end
figure;
plot(1:size(A,2),SqA,'ro');ylabel('sum of square of column in A');xlabel('ICs');
[B,I]=sort(SqA,'descend');


% Component 56 68 69 have the highest weights in the mixing matrics A
% topoplot to examine them
cd /home/zhibin/Documents/GitHub/Motor_cordination/EEGanalysis
load('Neuroscan_spherical__topoplot_chans.mat')
ComponentsExam=[50 4 36 94];
figure;
for i=1:length(ComponentsExam)
    subplot(length(ComponentsExam),1,i);
    topoplot(A(:,ComponentsExam(i)),test,'nosedir','+Y');title(['component' num2str(ComponentsExam(i))]);colorbar;
end
    
% subplot(3,1,1);topoplot(A(:,69),test,'nosedir','+Y');title('component 69');colorbar;
% subplot(3,1,2);topoplot(A(:,56),test,'nosedir','+Y');title('component 56');colorbar;
% subplot(3,1,3);topoplot(A(:,68),test,'nosedir','+Y');title('component 68'); colorbar;
% 
% figure;
% subplot(3,1,1);topoplot(A(:,91),test,'nosedir','+Y');title('component 91');colorbar;
% subplot(3,1,2);topoplot(A(:,76),test,'nosedir','+Y');title('component 76');colorbar;
% subplot(3,1,3);topoplot(A(:,19),test,'nosedir','+Y');title('component 19'); colorbar;
% 
% figure;
% subplot(3,1,1);topoplot(A(:,11),test,'nosedir','+Y');title('component 11');colorbar;
% subplot(3,1,2);topoplot(A(:,61),test,'nosedir','+Y');title('component 61');colorbar;
% subplot(3,1,3);topoplot(A(:,111),test,'nosedir','+Y');title('component 111'); colorbar;
% 
% figure;
% subplot(3,1,1);topoplot(A(:,3),test,'nosedir','+Y');title('component 3');colorbar;
% subplot(3,1,2);topoplot(A(:,53),test,'nosedir','+Y');title('component 53');colorbar;
% subplot(3,1,3);topoplot(A(:,101),test,'nosedir','+Y');title('component 101'); colorbar;

%% Display signal before ICA
% addpath /home/zhibin/Documents/GitHub/
% addpath(genpath('/home/zhibin/Documents/GitHub/'));

figure;
PlotStart=20;PlotEnd=25;
Duration=(IndStart(PlotStart)-1000):(IndEnds(PlotEnd)+1000);
TimesdataShort=1:size(EEGdataShort,1);
plotx(TimesdataShort(Duration),EEGdataShort(Duration,goodchans));
hold on;
for i=PlotStart:PlotEnd
    xline(IndStart(i),'r',{'Trial Start'});
    xline(IndEnds(i),'g',{'Trial End'});
end
hold off;
title('EEG Singal Before ICA');

%%  Plot ICA signal
figure;
subplot(2,1,2);
plotx(TimesdataShort(Duration),icasig(ComponentsExam,Duration));
hold on;
for i=PlotStart:PlotEnd
    xline(IndStart(i),'r',{'Trial Start'});
    xline(IndEnds(i),'g',{'Trial End'});
end
hold off;
legend(strsplit(num2str(ComponentsExam)));
title('IC time series');xlabel('time samples');ylabel('uV');
subplot(2,1,1);
plotx(TimesdataShort(Duration),icasig(:,Duration));
hold on;
for i=PlotStart:PlotEnd
    xline(IndStart(i),'r',{'Trial Start'});
    xline(IndEnds(i),'g',{'Trial End'});
end
hold off;
title('IC time series');ylabel('uV');

%% Calculate Correlation
% FP1 and FP2 are channel 1 and 3;
% compute correlation between FP1 and FP2;

[RHO1,PVAL1] = corr(EEGdataShort(:,1),icasig');
figure
subplot(1,2,1);
% yyaxis left
plot(1:length(RHO1),RHO1,'bo');hold on;ylabel('correlation coefficient');
% yyaxis right
plot(1:length(RHO1),PVAL1,'ro');ylabel('p value');xlabel('ICs');
legend({'correlation coefficient' 'p-values '});title('correlation with FP1');hold off;
subplot(1,2,2);
[RHO3,PVAL3] = corr(EEGdataShort(:,3),icasig');
% yyaxis left
plot(1:length(RHO1),RHO3,'bo');hold on;ylabel('correlation coefficient');
% yyaxis right
plot(1:length(RHO1),PVAL3,'ro');ylabel('p value');xlabel('ICs');
legend({'correlation coefficient' 'p-values '});title('correlation with FP2');hold off;

% component 56,3 seem to be highly correlated with FP1 and FP2
ComponentsExam=[50 63 36 62 108];
figure;
for i=1:length(ComponentsExam)
    subplot(length(ComponentsExam),1,i);
    topoplot(A(:,ComponentsExam(i)),test,'nosedir','+Y');title(['component' num2str(ComponentsExam(i))]);colorbar;
end
%% Plot ICs Topoplot Spetrogram ERP AND Power spectra,  (Similar to EEGLAB)
ComponentsExam=[50 36 62];

for i=1:length(ComponentsExam)
    ComponentPick=ComponentsExam(i);
    
    % Put together trial matrics
    icasig_trials=[];
    for i=1:length(goodepochs)
        icasig_trials(:,1,i)=icasig(ComponentPick,IndStart(i):IndEnds(i))';
    end

    % Just to check
    goodchans;

    figure('units','normalized','outerposition',[0 0 0.38 0.85]);
    % Topoplot
    subplot('Position',[0.05 0.55 0.4 0.4]);
    topoplot(A(:,ComponentPick),test,'nosedir','+Y');
    title(['component' num2str(ComponentPick)]);
    colorbar;

    % Plot Spectrogram
    subplot('Position',[0.55 0.65 0.4 0.3]);
    ColorLim=5;
    imagesc(squeeze(icasig_trials)');colormap jet; colorbar ; ylabel('Shorted Trials'); 
    caxis([-1*ColorLim ColorLim]);
    xline(500,'k');xline(1000,'k');
    title(['Component ' num2str(ComponentPick)]);

    % Plot IC ERP
    subplot('Position',[0.55 0.55 0.3 0.07]);
    plot(1:size(icasig_trials,1),mean(icasig_trials,3)); 
    xline(500,'k');xline(1000,'k');xlabel('time');

    % Plot Powerspetra
    subplot('Position',[0.1 0.05 0.8 0.4]);
    rate=Fs;maxfreq=50;
    [pow,freqs,df,eppow,corr,cprod,fcoef] = allspectra(icasig_trials,rate,maxfreq);
    plot(freqs,pow);xlabel('frequecy (Hz)');ylabel('Magnitude');
end

%% Deside which components to remove and mix back the signal and display
ComponentRemove=[50 36 62];
A(:,ComponentRemove)=0; icasig(ComponentRemove,:)=0;

mixedsig=A*icasig;

% Plot before and after
PlotStart=20;PlotEnd=25;
Duration=(IndStart(PlotStart)-1000):(IndEnds(PlotEnd)+1000);
TimesdataShort=1:size(EEGdataShort,1);

figure;
subplot(2,1,1);
% before ICA
plotx(TimesdataShort(Duration),EEGdataShort(Duration,:));
hold on;
for i=PlotStart:PlotEnd
    xline(IndStart(i)+500,'g-',{'Target'});
    xline(IndEnds(i)-1000,'g.-',{'Go'});
end
hold off;
title('EEG Singal Before ICA');

% after ICA
subplot(2,1,2);
% PlotStart=95;PlotEnd=100;
% Duration=(IndEvents(PlotStart)-1000):(IndEnds(PlotEnd)+1000);
% TimesdataShort=1:size(EEGdataShort,1);
plotx(TimesdataShort(Duration),mixedsig(:,Duration));
hold on;
for i=PlotStart:PlotEnd
    xline(IndStart(i)+500,'g-',{'Target'});
    xline(IndEnds(i)-1000,'g.-',{'Go'});
end
hold off;
title('Mixed Signal with ICs removed');



