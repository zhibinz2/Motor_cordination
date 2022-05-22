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
plot(detrend_dataL);

figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(detrend_dataR);

%% filtfilthd method (hnl) high pass first then low pass
cd /usr/local/MATLAB/R2019a/toolbox/signal/signal
which filtfilt
% high pass (no paddings needed)
Hd = makefilter(srL,0.2,0.15,6,20,0); 
filtered_dataL1=filtfilthd(Hd,detrend_dataL);
filtered_dataR1=filtfilthd(Hd,detrend_dataR);
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(filtered_dataL1);
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(filtered_dataR1);

% add padding
paddingL=zeros(round(size(filtered_dataL1,1)/10), size(filtered_dataL1,2));
filtered_dataL2=cat(1,paddingL,filtered_dataL1,paddingL);
paddingR=zeros(round(size(filtered_dataR1,1)/10), size(filtered_dataR1,2));
filtered_dataR2=cat(1,paddingR,filtered_dataR1,paddingR);
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(filtered_dataL2);
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(filtered_dataR2);

% low pass 
% (this will create short edge artifact)
% (if added zero paddings, edge artifact disappear)
% (remove existing filtered_data variable from workspace might fasten)
% filtered_dataL3=[];
Hd = makefilter(srL,50,51,6,20,0);  
filtered_dataL3=filtfilthd(Hd,filtered_dataL2);
filtered_dataR3=filtfilthd(Hd,filtered_dataR2);
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(filtered_dataL3);
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(filtered_dataR3);
% ylim([-100 100]);

% remove padding
filtered_dataL4=filtered_dataL3((size(paddingL,1)+1):(size(paddingL,1)+size(detrend_dataL,1)),:);
filtered_dataR4=filtered_dataR3((size(paddingR,1)+1):(size(paddingR,1)+size(detrend_dataR,1)),:);
figure('units','normalized','outerposition',[0 0 1 0.3]);
plotx(filtered_dataL4);
figure('units','normalized','outerposition',[0 0 1 0.3]);
plotx(filtered_dataR4);
% ylim([-100 100]);

clearvars samplesL filtered_dataL1 filtered_dataL2 filtered_dataL3
clearvars samplesR filtered_dataR1 filtered_dataR2 filtered_dataR3
%% ICA
filtered_dataL4; filtered_dataR4
% Left player
% run ICA
[icasigL, AL, WL] = fastica(filtered_dataL4');
% Calculate Correlation
% FP1 and FP2 are channel 1 and 3;
% compute correlation between FP1 and FP2;
[RHO1,PVAL1] = corr(filtered_dataL4(:,1),icasigL');
[RHO3,PVAL3] = corr(filtered_dataL4(:,3),icasigL');
[B1,I1]=sort(abs(RHO1),'descend');[B3,I3]=sort(abs(RHO3),'descend');
ComponentsExamL=unique([I1(1) I3(1)]);
% topoplot to examine them
cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/channels_info
load('chaninfo.mat')
figure;
for i=1:length(ComponentsExamL)
    subplot(length(ComponentsExamL),1,i);
    topoplot(AL(:,ComponentsExamL(i)),chaninfo,'nosedir','+X');title(['component' num2str(ComponentsExamL(i))]);colorbar;
end

% Right player
% run ICA
[icasigR, AR, WR] = fastica(filtered_dataR4');
% Calculate Correlation
% FP1 and FP2 are channel 1 and 3;
% compute correlation between FP1 and FP2;
[RHO1,PVAL1] = corr(filtered_dataR4(:,1),icasigR');
[RHO3,PVAL3] = corr(filtered_dataR4(:,3),icasigR');
[B1,I1]=sort(abs(RHO1),'descend');[B3,I3]=sort(abs(RHO3),'descend');
ComponentsExamR=unique([I1(1) I3(1)]);
% topoplot to examine them
figure;
for i=1:length(ComponentsExamR)
    subplot(length(ComponentsExamR),1,i);
    topoplot(AR(:,ComponentsExamR(i)),chaninfo,'nosedir','+X');title(['component' num2str(ComponentsExamR(i))]);colorbar;
end

% Remove components and mix back the signal and display
% Left player
ComponentRemoveL=ComponentsExamL;
AL(:,ComponentRemoveL)=0; icasigL(ComponentRemoveL,:)=0;
mixedsigL=AL*icasigL;
% Plot before and after
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,1,1);
% before ICA
plotx(timeL,filtered_dataL4);
hold on;hold off;%ylim([-100 100]);
title('EEG Singal Before ICA');
% after ICA
subplot(2,1,2);
plotx(timeL,mixedsigL);
hold on;hold off;%ylim([-100 100]);
title('Mixed Signal with ICs removed');

% Right player
ComponentRemoveR=ComponentsExamR;
AR(:,ComponentRemoveR)=0; icasigR(ComponentRemoveR,:)=0;
mixedsigR=AR*icasigR;
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
[x, ~] = ginput(2); % read two mouse clicks on the plot % x were index, y were real values
% get the proximate index
string(x)
ind1=round(x(1))
ind2=round(x(2))
% zoom in and plot again
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,1,1);
% before ICA
plotx(timeR,filtered_dataR4);
hold on;hold off;ylim([-100 100]);xlim([ind1 ind2]);
title('EEG Singal Before ICA');
% after ICA
subplot(2,1,2);
plotx(timeR,mixedsigR);
hold on;hold off;ylim([-100 100]);xlim([ind1 ind2]);
title('Mixed Signal with ICs removed');


%% Segment EEG
mixedsigL=mixedsigL';
mixedsigR=mixedsigR';
% PacersL=PacerTimeIndL([1 2 3 4 5 34 35 36 37 66 67 68 69 98 99 100 101 130 131 132 133 162 163 164 165 194]);
% PacersR=PacerTimeIndR([1 2 3 4 5 34 35 36 37 66 67 68 69 98 99 100 101 130 131 132 133 162 163 164 165 194]);
OpenEyeRestingL=mixedsigL(PacersL(1):PacersL(2),:);
CloseEyeRestingL=mixedsigL(PacersL(3):PacersL(4),:);
Condi1L=mixedsigL(PacersL(6):PacersL(7),:);
Resting2L=mixedsigL(PacersL(7):PacersL(8),:);
Condi2L=mixedsigL(PacersL(10):PacersL(11),:);
Resting3L=mixedsigL(PacersL(11):PacersL(12),:);
Condi3L=mixedsigL(PacersL(14):PacersL(15),:);
Resting4L=mixedsigL(PacersL(15):PacersL(16),:);
Condi4L=mixedsigL(PacersL(18):PacersL(19),:);
Resting5L=mixedsigL(PacersL(19):PacersL(20),:);
Condi5L=mixedsigL(PacersL(22):PacersL(23),:);
Resting6L=mixedsigL(PacersL(23):PacersL(24),:);
Condi6L=mixedsigL(PacersL(26):end,:);



