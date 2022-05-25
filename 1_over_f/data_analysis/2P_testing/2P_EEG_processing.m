%% plot raw EEG
timeL;samplesL;TRIGGERindL;srL;channels_infoL;
timeR;samplesR;TRIGGERindR;srR;channels_infoR;

figure;
subplot(2,1,1);plot(samplesL(1:32,:)');
subplot(2,1,2);plot(samplesR(1:32,:)');

figure; % look at the bad chan P8
subplot(2,1,1);plot(samplesL(28,:)');
subplot(2,1,2);plot(samplesR(28,:)');

%% Extract EEG
EEGL=samplesL(1:32,:)';
EEGR=samplesR(1:32,:)';

%% detrend the EEG data (no padding needed)
detrend_dataL=detrend(EEGL,2);
detrend_dataR=detrend(EEGR,2);

figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(detrend_dataL); ylim([-1000 1000]);

figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(detrend_dataR); ylim([-1000 1000]);

%% filtfilthd method (hnl) high pass first then low pass
cd /usr/local/MATLAB/R2019a/toolbox/signal/signal
which filtfilt
which filtfilt -all
cd D:\Program Files\MATLAB\R2019b\toolbox\signal\signal\filtfilt.m   
% high pass (no paddings needed)
Hd = makefilter(srL,0.2,0.15,6,20,0); 
filtered_dataL1=filtfilthd(Hd,detrend_dataL);
filtered_dataR1=filtfilthd(Hd,detrend_dataR);

figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(filtered_dataL1); ylim([-300 300]);
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(filtered_dataR1); ylim([-300 300]);

% add padding (this step takes several minutes)
paddingL=zeros(round(size(filtered_dataL1,1)/10), size(filtered_dataL1,2));
filtered_dataL2=cat(1,paddingL,filtered_dataL1,paddingL);
paddingR=zeros(round(size(filtered_dataR1,1)/10), size(filtered_dataR1,2));
filtered_dataR2=cat(1,paddingR,filtered_dataR1,paddingR);
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(filtered_dataL2);ylim([-300 300]);title('filtered-dataL2');
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(filtered_dataR2);ylim([-300 300]);title('filtered-dataR2');

% low pass 
% (this will create short edge artifact)
% (if added zero paddings, edge artifact disappear)
% (remove existing filtered_data variable from workspace might fasten)
% filtered_dataL3=[];
Hd = makefilter(srL,50,51,6,20,0);  
filtered_dataL3=filtfilthd(Hd,filtered_dataL2);
filtered_dataR3=filtfilthd(Hd,filtered_dataR2);
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(filtered_dataL3);ylim([-300 300]);
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(filtered_dataR3);ylim([-300 300]);
% plotx(filtered_dataR3(ind1:ind2,:));ylim([-300 300]); % channel 30 bad

% remove padding
filtered_dataL4=filtered_dataL3((size(paddingL,1)+1):(size(paddingL,1)+size(detrend_dataL,1)),:);
filtered_dataR4=filtered_dataR3((size(paddingR,1)+1):(size(paddingR,1)+size(detrend_dataR,1)),:);
figure('units','normalized','outerposition',[0 0 1 0.3]);
plotx(filtered_dataL4);ylim([-300 300]);
figure('units','normalized','outerposition',[0 0 1 0.3]);
plotx(filtered_dataR4);ylim([-300 300]);
% ylim([-100 100]);

clearvars filtered_dataL1 filtered_dataL2 filtered_dataL3
clearvars filtered_dataR1 filtered_dataR2 filtered_dataR3

%% ICA
filtered_dataL4; filtered_dataR4;
% Left player
% run ICA
[icasigL, AL, WL] = fastica(filtered_dataL4');
% Plot ICA component
for i=1:size(AL,2)
    SqAL(i)=sumsqr(AL(:,i));
end
figure;
plot(1:size(AL,2),SqAL,'ro');ylabel('sum of square of column in A');xlabel('ICs');
[B,I]=sort(SqAL,'descend');
ComponentsExam=I(1:10);
% topoplot to examine them
cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/channels_info
load('chaninfo.mat')
figure;
for i=1:length(ComponentsExam)
    subplot(5,2,i);
    topoplot(AL(:,ComponentsExam(i)),chaninfo,'nosedir','+X');title(['component' num2str(ComponentsExam(i))]);colorbar;
end
% Calculate Correlation
% FP1 and FP2 are channel 1 and 3;
% compute correlation between FP1 and FP2;
[RHO1,PVAL1] = corr(filtered_dataL4(:,1),icasigL');
[RHO3,PVAL3] = corr(filtered_dataL4(:,3),icasigL');
[B1,I1]=sort(abs(RHO1),'descend');[B3,I3]=sort(abs(RHO3),'descend');
ComponentsExamL=unique([I1(1:2) I3(1:2)]);
figure;
for i=1:length(ComponentsExamL)
    subplot(length(ComponentsExamL),1,i);
    topoplot(AL(:,ComponentsExamL(i)),chaninfo,'nosedir','+X');title(['component' num2str(ComponentsExamL(i))]);colorbar;
end

% Right player
% run ICA
[icasigR, AR, WR] = fastica(filtered_dataR4');
% Plot ICA component
for i=1:size(AR,2)
    SqAR(i)=sumsqr(AR(:,i));
end
figure;
plot(1:size(AR,2),SqAR,'ro');ylabel('sum of square of column in A');xlabel('ICs');
[B,I]=sort(SqAR,'descend');
ComponentsExam=I(1:10);
figure;
for i=1:length(ComponentsExam)
    subplot(5,2,i);
    topoplot(AR(:,ComponentsExam(i)),chaninfo,'nosedir','+X');title(['component' num2str(ComponentsExam(i))]);colorbar;
end
% Calculate Correlation
% FP1 and FP2 are channel 1 and 3;
% compute correlation between FP1 and FP2;
[RHO1,PVAL1] = corr(filtered_dataR4(:,1),icasigR');
[RHO3,PVAL3] = corr(filtered_dataR4(:,3),icasigR');
[B1,I1]=sort(abs(RHO1),'descend');[B3,I3]=sort(abs(RHO3),'descend');
ComponentsExamR=unique([I1(1:2) I3(1:2)]);
% topoplot to examine them
figure;
for i=1:length(ComponentsExamR)
    subplot(length(ComponentsExamR),1,i);
    topoplot(AR(:,ComponentsExamR(i)),chaninfo,'nosedir','+X');title(['component' num2str(ComponentsExamR(i))]);colorbar;
end

% Remove components and mix back the signal and display the comparison
% (before vs after)
% Left player
ComponentRemoveL=ComponentsExamL;
AL(:,ComponentRemoveL)=0; icasigL(ComponentRemoveL,:)=0;
mixedsigL=AL*icasigL;
mixedsigL=mixedsigL';
% Plot before and after
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,1,1);
% before ICA
plotx(timeL,filtered_dataL4);
hold on;hold off; ylim([-200 200]);
title('EEG Singal Before ICA');
% after ICA
subplot(2,1,2);
plotx(timeL,mixedsigL);
hold on;hold off; ylim([-200 200]);
title('Mixed Signal with ICs removed');
% selet two point on x axis to zoom in
ind1=round(length(timeL)*0.55);
ind2=round(length(timeL)*0.6);
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,1,1);
% before ICA
plotx(timeL(ind1:ind2),filtered_dataL4(ind1:ind2,:));
hold on;hold off; ylim([-200 200]);
title('EEG Singal Before ICA');
% after ICA
subplot(2,1,2);
plotx(timeL(ind1:ind2),mixedsigL(ind1:ind2,:));
hold on;hold off; ylim([-200 200]);
title('Mixed Signal with ICs removed');

% Right player
ComponentRemoveR=ComponentsExamR;
AR(:,ComponentRemoveR)=0; icasigR(ComponentRemoveR,:)=0;
mixedsigR=AR*icasigR;
mixedsigR=mixedsigR';
% Plot before and after
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,1,1);
% before ICA
plotx(timeR,filtered_dataR4);
hold on;hold off;ylim([-100 100]);
title('EEG Singal Before ICA');
% after ICA
subplot(2,1,2);
plotx(timeR,mixedsigR);
hold on;hold off;ylim([-100 100]);
title('Mixed Signal with ICs removed');
% selet two point on x axis to zoom in
ind1=round(length(timeL)*0.55);
ind2=round(length(timeL)*0.6);
% zoom in and plot again
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,1,1);
% before ICA
plotx(timeR(ind1:ind2),filtered_dataR4(ind1:ind2,:));
% plotx(timeR(ind1:ind2),filtered_dataR4(ind1:ind2,1:29)); 
hold on;hold off;ylim([-100 100]);
title('EEG Singal Before ICA');
% after ICA
subplot(2,1,2);
plotx(timeR(ind1:ind2),mixedsigR(ind1:ind2,:));
% plotx(timeR(ind1:ind2),mixedsigR(ind1:ind2,1:29));

hold on;hold off;ylim([-100 100]);
title('Mixed Signal with ICs removed');

%% examine spectra on scalp map
Timeselected=timeL(ind1:ind2);% Timeselected=timeR(ind1:ind2); 
EEGselected=filtered_dataL4(ind1:ind2,:);% EEGselected=filtered_dataR4(ind1:ind2,:);
EEGselected=mixedsigL(ind1:ind2,:);% EEGselected=mixedsigR(ind1:ind2,:);
% fft
fcoef=fft(EEGselected);
N=size(EEGselected,1); % N=length(Timeselected);
fcoef=fcoef/N;
halfN=floor(N/2);
% df=1/T; 
% fV=[0:df:(halfN-1)*df]; 
fV=linspace(0,srR/2,halfN+1); % same df
fcoef=2*fcoef(1:halfN,:);
amplitude = abs(fcoef);
% Plot on scalp map for the spectra (for my own examing)
figure('units','normalized','outerposition',[0 0 1 1]);
for chan=1:32
    subplot('Position',[XXPLOT(chan) YYPLOT(chan) 0.05 0.05]); % not showing, why
    % plot([1 1],[1 1],'ro');
    plot(fV(1:size(amplitude,1)),amplitude(:,chan));
    if ~isempty(find([1:30 32]==chan))
    set(gca,'XTick',[]); % set(gca,'YTick',[]); 
    end
    if chan==31
        xlabel('frequency');
        ylabel('amplitude (uV)');
    end
    xlim([0 25]);title([labels{chan}]);%ylim([0 1]);
end
suptitle('spectra of all channels on scalp map')


%% Segment EEG
% mixedsigL=mixedsigL';
% mixedsigR=mixedsigR';
% PacersL=PacerTimeIndL([1 2 3 4 5 34 35 36 37 66 67 68 69 98 99 100 101 130 131 132 133 162 163 164 165 194]);
% PacersR=PacerTimeIndR([1 2 3 4 5 34 35 36 37 66 67 68 69 98 99 100 101 130 131 132 133 162 163 164 165 194]);
EEGOpenEyeRestingL=mixedsigL(PacersL(1):PacersL(2),:);
EEGCloseEyeRestingL=mixedsigL(PacersL(3):PacersL(4),:);
EEGCondi1L=mixedsigL(PacersL(6):PacersL(7),:);
EEGResting2L=mixedsigL(PacersL(7):PacersL(8),:);
EEGCondi2L=mixedsigL(PacersL(10):PacersL(11),:);
EEGResting3L=mixedsigL(PacersL(11):PacersL(12),:);
EEGCondi3L=mixedsigL(PacersL(14):PacersL(15),:);
EEGResting4L=mixedsigL(PacersL(15):PacersL(16),:);
EEGCondi4L=mixedsigL(PacersL(18):PacersL(19),:);
EEGResting5L=mixedsigL(PacersL(19):PacersL(20),:);
EEGCondi5L=mixedsigL(PacersL(22):PacersL(23),:);
EEGResting6L=mixedsigL(PacersL(23):PacersL(24),:);
EEGCondi6L=mixedsigL(PacersL(26):end,:);

EEGOpenEyeRestingR=mixedsigR(PacersR(1):PacersR(2),:);
EEGCloseEyeRestingR=mixedsigR(PacersR(3):PacersR(4),:);
EEGCondi1R=mixedsigR(PacersR(6):PacersR(7),:);
EEGResting2R=mixedsigR(PacersR(7):PacersR(8),:);
EEGCondi2R=mixedsigR(PacersR(10):PacersR(11),:);
EEGResting3R=mixedsigR(PacersR(11):PacersR(12),:);
EEGCondi3R=mixedsigR(PacersR(14):PacersR(15),:);
EEGResting4R=mixedsigR(PacersR(15):PacersR(16),:);
EEGCondi4R=mixedsigR(PacersR(18):PacersR(19),:);
EEGResting5R=mixedsigR(PacersR(19):PacersR(20),:);
EEGCondi5R=mixedsigR(PacersR(22):PacersR(23),:);
EEGResting6R=mixedsigR(PacersR(23):PacersR(24),:);
EEGCondi6R=mixedsigR(PacersR(26):end,:);


%% segment button presses
BottonPresTimeL01; BottonPresTimeR01;

BPOpenEyeRestingL=BottonPresTimeL01(PacersL(1):PacersL(2),:);
BPCloseEyeRestingL=BottonPresTimeL01(PacersL(3):PacersL(4),:);
BPCondi1L=BottonPresTimeL01(PacersL(6):PacersL(7),:);
BPResting2L=BottonPresTimeL01(PacersL(7):PacersL(8),:);
BPCondi2L=BottonPresTimeL01(PacersL(10):PacersL(11),:);
BPResting3L=BottonPresTimeL01(PacersL(11):PacersL(12),:);
BPCondi3L=BottonPresTimeL01(PacersL(14):PacersL(15),:);
BPResting4L=BottonPresTimeL01(PacersL(15):PacersL(16),:);
BPCondi4L=BottonPresTimeL01(PacersL(18):PacersL(19),:);
BPResting5L=BottonPresTimeL01(PacersL(19):PacersL(20),:);
BPCondi5L=BottonPresTimeL01(PacersL(22):PacersL(23),:);
BPResting6L=BottonPresTimeL01(PacersL(23):PacersL(24),:);
BPCondi6L=BottonPresTimeL01(PacersL(26):end,:);

BPOpenEyeRestingR=BottonPresTimeR01(PacersR(1):PacersR(2),:);
BPCloseEyeRestingR=BottonPresTimeR01(PacersR(3):PacersR(4),:);
BPCondi1R=BottonPresTimeR01(PacersR(6):PacersR(7),:);
BPResting2R=BottonPresTimeR01(PacersR(7):PacersR(8),:);
BPCondi2R=BottonPresTimeR01(PacersR(10):PacersR(11),:);
BPResting3R=BottonPresTimeR01(PacersR(11):PacersR(12),:);
BPCondi3R=BottonPresTimeR01(PacersR(14):PacersR(15),:);
BPResting4R=BottonPresTimeR01(PacersR(15):PacersR(16),:);
BPCondi4R=BottonPresTimeR01(PacersR(18):PacersR(19),:);
BPResting5R=BottonPresTimeR01(PacersR(19):PacersR(20),:);
BPCondi5R=BottonPresTimeR01(PacersR(22):PacersR(23),:);
BPResting6R=BottonPresTimeR01(PacersR(23):PacersR(24),:);
BPCondi6R=BottonPresTimeR01(PacersR(26):end,:);

%% segment feedbacks from the other 
FeedbTimeL01; FeedbTimeR01;

FBOpenEyeRestingL=FeedbTimeL01(PacersL(1):PacersL(2),:);
FBCloseEyeRestingL=FeedbTimeL01(PacersL(3):PacersL(4),:);
FBCondi1L=FeedbTimeL01(PacersL(6):PacersL(7),:);
FBResting2L=FeedbTimeL01(PacersL(7):PacersL(8),:);
FBCondi2L=FeedbTimeL01(PacersL(10):PacersL(11),:);
FBResting3L=FeedbTimeL01(PacersL(11):PacersL(12),:);
FBCondi3L=FeedbTimeL01(PacersL(14):PacersL(15),:);
FBResting4L=FeedbTimeL01(PacersL(15):PacersL(16),:);
FBCondi4L=FeedbTimeL01(PacersL(18):PacersL(19),:);
FBResting5L=FeedbTimeL01(PacersL(19):PacersL(20),:);
FBCondi5L=FeedbTimeL01(PacersL(22):PacersL(23),:);
FBResting6L=FeedbTimeL01(PacersL(23):PacersL(24),:);
FBCondi6L=FeedbTimeL01(PacersL(26):end,:);

FBOpenEyeRestingR=FeedbTimeR01(PacersR(1):PacersR(2),:);
FBCloseEyeRestingR=FeedbTimeR01(PacersR(3):PacersR(4),:);
FBCondi1R=FeedbTimeR01(PacersR(6):PacersR(7),:);
FBResting2R=FeedbTimeR01(PacersR(7):PacersR(8),:);
FBCondi2R=FeedbTimeR01(PacersR(10):PacersR(11),:);
FBResting3R=FeedbTimeR01(PacersR(11):PacersR(12),:);
FBCondi3R=FeedbTimeR01(PacersR(14):PacersR(15),:);
FBResting4R=FeedbTimeR01(PacersR(15):PacersR(16),:);
FBCondi4R=FeedbTimeR01(PacersR(18):PacersR(19),:);
FBResting5R=FeedbTimeR01(PacersR(19):PacersR(20),:);
FBCondi5R=FeedbTimeR01(PacersR(22):PacersR(23),:);
FBResting6R=FeedbTimeR01(PacersR(23):PacersR(24),:);
FBCondi6R=FeedbTimeR01(PacersR(26):end,:);

