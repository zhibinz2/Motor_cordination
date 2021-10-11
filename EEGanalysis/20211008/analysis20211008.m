
%% load behaviral data
clear;
cd('/home/zhibin/Documents/GitHub/Motor_cordination/data')
load('matlab2021100802.mat');
%% Reorgainize behavioral data into BehavData and condition
subjectnumber=behaviraldata.subjectnumber;
dataBlock=behaviraldata.dataBlock;

% Just to try out with one trial and plot
% dataTrialL=dataBlock.dataTrialL;  
% dataTrialR=dataBlock.dataTrialR;  
% dataTrialJ=dataBlock.dataTrialJ;

% xLyL=dataTrialL.xLyL; plot(xLyL(:,1),xLyL(:,2),'r'); hold on;
% xRyR=dataTrialR.xRyR; plot(xRyR(:,1),xRyR(:,2),'b');
% xJyJ=dataTrialJ.xJyJ; plot(xJyJ(:,1),xJyJ(:,2),'g'); hold off;
% set(gca, 'YDir', 'reverse');xlim([0,screenXpixels]);ylim([0,screenYpixels]);


blockNumber=dataBlock.blockNumber;
xLyL=dataTrialL.xLyL;
condition=dataTrialL.condition;

numBlock=length(dataBlock);
BehavData={};
CondiData=[]; 

for b=1:numBlock
    for t=1:length(behaviraldata.dataBlock(b).dataTrialL)
        dataxLyLtemp=behaviraldata.dataBlock(b).dataTrialL(t).xLyL;
        dataxRyRtemp=behaviraldata.dataBlock(b).dataTrialR(t).xRyR;
        dataxJyJtemp=behaviraldata.dataBlock(b).dataTrialJ(t).xJyJ;
        dataxLyLxRyRxJyJtemp=[dataxLyLtemp dataxRyRtemp dataxJyJtemp];
        dataConditiontemp=behaviraldata.dataBlock(b).dataTrialL(t).condition;
        BehavData=[BehavData; {dataxLyLxRyRxJyJtemp}];
        CondiData=[CondiData; dataConditiontemp];% same as allPerm
    end
end

% ans=BehavData{1};
% CondiData

NumTrialsRecorded=length(CondiData);

%% Check EEG on 20210916
% load EEG all as one
cd('/home/zhibin/Documents/Acquisition/zhibin20211008')
EEG=loadcurry([pwd '/Bimanual_Reach_zhibin_20211008.cdt']);
dataEEG=EEG.data;
datatimes=EEG.times;
Fs=EEG.srate;
% plot(datatimes,dataEEG);

%% photocell processing to find trial start and end index
% plot photocell analog signal
analog1=133; 
analog1data=dataEEG(analog1,:);
% plot(datatimes,analog1data,'b'); % might need to mannually cut off the start and end of block

% hold on;
% analog2=134; 
% analog2data=dataEEG(analog2,:);
% plot(datatimes,analog2data,'r'); % weird start and end of block signal
% hold off;

% close all
Halfhigh1=1/2*(max(analog1data)-min(analog1data));
% Halfhigh2=1/2*(max(analog2data)-min(analog2data));

%***************************** Use the analog1data 
% locate the trial sessions % pks=value of the peak % locs=time of the peak
[pks,locs] = findpeaks(analog1data,datatimes,'MinPeakProminence',Halfhigh1,'Annotate','extents');

locsDiff=diff(locs);% time distances between peaks
% plot(locsDiff,'ro');% look at distribution of these distances
ind_locsDiff=find([locsDiff, inf]>5*mean(locsDiff)); % pks distance indices at the end of each trial 
trialLength=diff([0 ind_locsDiff]); % number of peaks in each trial
% plot(trialLength,'ro'); % look at trial length
ind_trial_end=cumsum(trialLength); % index of the last peak in each trial (index in locs)
ind_trial_start=ind_trial_end-trialLength+1; % index of the first peak in each trial (index in locs)

% % checking starting peak of the trials (index in locs)
% plot(ind_trial_start,zeros(length(ind_trial_start)),'go'); 
% hold on;
% plot(ind_trial_end,ones(length(ind_trial_end)),'ro');
% ylim([-20 20]);

%***************************** Use the analog2data 
% % locate the trial sessions % pks=value of the peak % locs=time of the peak
% [pks,locs] = findpeaks(analog2data,datatimes,'MinPeakProminence',Halfhigh2,'Annotate','extents');
% 
% locsDiff=diff(locs);% time distances between peaks
% % plot(locsDiff,'ro'); % look at distribution of these distances
% ind_locsDiff=find([locsDiff, inf]>(mean(locsDiff)+mean(locsDiff)/5)); % indices of pks distances at the end of each trial 
% trialLength=diff([0 ind_locsDiff]); % number of peaks in each trial
% % plot(trialLength,'ro'); % look at trial length 
% ind_trial_end=cumsum(trialLength); % index of the last peak in each trial (index in locs)
% ind_trial_start=ind_trial_end-trialLength+1; % index of the first peak in each trial (index in locs)

%% detrend the data

addpath(genpath('/home/zhibin/Documents/GitHub'));
% data=dataEEG(1:128,2000000:2002000);
EEGdata=dataEEG(1:128,:);

% tranpose the data so that time is the first dimension, channel is 2nd dimension, and trial is the 3rd dimension
transpose_data=permute(EEGdata,[2,1]);
% plotx(mean(transpose_data,3));
% plotx(transpose_data);

% detrend the data
detrend_data=ndetrend(transpose_data,1); 
% plotx(mean(detrend_data,2)); xlim([5000 10000]);
% plotx(detrend_data); xlim([5000 10000]);

%% Broadband filter
% apply high pass filter to filter out frequencies lower than 1;
% load eeglab toolbox to avoid filtered_data becoming NaN
eeglab;
close all;
Hd = makefilter(Fs,1,0.5,3,20,0); % Astop 6 or 20;
filtered_data=filtfilthd(Hd,double(detrend_data));

% remove channel 91 44
% filtered_data(:,91)=[];filtered_data(:,44)=[];
%  plotx(filtered_data); xlim([5000 10000]);

% apply low pass filter to filter out frequencies higher than 50;
Hd = makefilter(Fs,50,55,3,20,0);
filtered_data=filtfilthd(Hd,filtered_data);

% plotx(mean(filtered_data,2));  
% plotx(filtered_data);xlim([5000 10000]);

%% fiter5bands 
run filter5bands.m

% filtered_data=filtered_delta;
% filtered_data=filtered_theta;
% filtered_data=filtered_alpha;
% filtered_data=filtered_beta;
% filtered_data=filtered_gamma;

% plot(filtered_beta(5000:10000,8));
% figure;
% plot(filtered_gamma(5000:10000,8));
%% organize filtered_data into data_trials=sample_timepoints x (139+4 channels) x  trial  
% Fs=2000;ifi=0.0167;
NumTrialsRecorded=length(ind_trial_start);
NumTrialtimepoints=round(max(trialLength)*2*ifi*Fs);%number of samples in a trial,ideally
NumEEGChannels=size(filtered_data,2);
NumAllChannels=NumEEGChannels+6; % +4;

% Initialize the data_trials matrix
% Zeros at the end of the trial session if trials shorter than is should be, ideally
% Reality is, sometimes this trial length do not cover the length of some trial
data_trials=zeros(NumTrialtimepoints,NumAllChannels,NumTrialsRecorded);

% Initialize the indices for trial start and end time point
trialstart=zeros(1,NumTrialsRecorded);trialend=zeros(1,NumTrialsRecorded);

% timing (trial length) might mess up
for ntr=1:NumTrialsRecorded % take too much memory, matlab might quit
    trialstart(ntr)=find(datatimes==locs(ind_trial_start(ntr)));
    trialend(ntr)  =find(datatimes==locs(ind_trial_end(ntr)));
    data_trials(1:(trialend(ntr)-trialstart(ntr)+1),1:NumEEGChannels,ntr)=filtered_data(trialstart(ntr):trialend(ntr),1:NumEEGChannels); 
end

% examine 
% plot(ind_trial_end-ind_trial_start,'ro');ylabel('trial peaks');xlabel('trial');
% locsstart=locs(ind_trial_start);
% locsend=locs(ind_trial_end);
% plot(trialend-trialstart,'ro');ylabel('trial length');xlabel('trial');
% plot(trialstart,ones(1,length(trialstart)),'go');hold on;plot(trialend,ones(1,length(trialend)),'ro'); % check 1 second inter-trial interval
% plot(data_trials(:,1:128,5));
% plot(datatimes)

%% Resample BehavData
RSBehavData={};
% block=1;
% for ntr=NumTrialsRecorded*(block-1)+1:NumTrialsRecorded*(block-1)+NumTrialsRecorded
for ntr=1:NumTrialsRecorded
    xLyLxRyRxJyJtemp=BehavData{ntr};
    shift=xLyLxRyRxJyJtemp(1,:); % shift to the origin (start from zero to minimized edge artifact)
    xLyLxRyRxJyJtempshift=xLyLxRyRxJyJtemp-shift.*ones(size(xLyLxRyRxJyJtemp,1),1);
    
    padstart=xLyLxRyRxJyJtempshift(1,:).*ones(size(xLyLxRyRxJyJtempshift,1),1); % padding before the start of trial
    padtend=xLyLxRyRxJyJtempshift(end,:).*ones(size(xLyLxRyRxJyJtempshift,1),1); % padding after
    xLyLxRyRxJyJtempshiftpad=[padstart;xLyLxRyRxJyJtempshift;padtend]; % concatenate
    RSxLyLxRyRxJyJtempshiftpad=resample(xLyLxRyRxJyJtempshiftpad,size(data_trials,1),length(xLyLxRyRxJyJtemp)); % resample
    RSxLyLxRyRxJyJtempshift=RSxLyLxRyRxJyJtempshiftpad(size(data_trials,1)+1:2*size(data_trials,1),:);% remove the padding
    RSxLyLxRyRxJyJtemp=RSxLyLxRyRxJyJtempshift+shift.*ones(size(RSxLyLxRyRxJyJtempshift,1),1); % shift back to the original position
    
    RSxLyLxRyRxJyJtemp=RSxLyLxRyRxJyJtemp'; % transpose to fit into trial EEG matrix
    RSBehavData{ntr}={RSxLyLxRyRxJyJtemp};
end

% % examine RSBehavData
% ans=RSBehavData{1};
% ans=RSBehavData{51};


%% integrate BehavData and CondiData into data_trials: trials x (139+4 channels) x sample_timepoints;

% for ntr=NumTrialsRecorded*(block-1)+1:NumTrialsRecorded*(block-1)+NumTrialsRecorded
for ntr=1:NumTrialsRecorded
    data_trials(:,NumEEGChannels+1:NumEEGChannels+6,ntr)=[cell2mat(RSBehavData{ntr})]'; % replace the last four channels
end

% % Just to check out
CondiData;
figure;plot(CondiData);
% yticks([1 2 3 4]);
% yticklabels({'A1&A2','A1&A1','A2&A1','A2&A2'});title('task sequence'); xlabel('trial');
% unique(CondiData)
% allPerm;
% conditionfuctions;
% conditionNames={'A1&A2','A1&A1','A2&A1','A2&A2'};
conditionNames=conditions;
% conditionNames{1};

%% draw condition trace

% plot the mouse trace
UniCondi=unique(CondiData);

for u=1:length(UniCondi)% u=2
    
    indtemp=find(CondiData==UniCondi(u));
    
    % pick a condition from randomized set allPerm
    conditionSelected = UniCondi(u);
    % produced the position parameters
        rad_ang=conditions(conditionSelected);
        [x,y] = drawReach(radius,rad_ang, xCenter, yCenter);
            
    
    subplot(1,length(UniCondi),u);
    
    plot(x,y,'-', 'color', grey, 'LineWidth',10);
    hold on;

    for indt=1:length(indtemp)
    TracexL=squeeze(data_trials(:,NumEEGChannels+1,indtemp(indt)));
    TraceyL=squeeze(data_trials(:,NumEEGChannels+2,indtemp(indt)));
    TracexR=squeeze(data_trials(:,NumEEGChannels+3,indtemp(indt)));
    TraceyR=squeeze(data_trials(:,NumEEGChannels+4,indtemp(indt)));
    TracexJ=squeeze(data_trials(:,NumEEGChannels+5,indtemp(indt)));
    TraceyJ=squeeze(data_trials(:,NumEEGChannels+6,indtemp(indt)));
    plot(TracexL,TraceyL,'-r',TracexR,TraceyR,'-b',TracexJ,TraceyJ,'-g');
    end
    
    set(gca, 'YDir', 'reverse');
    hold off;
    
    xlim([0 screenXpixels]);ylim([0 screenYpixels]);
    %title(['condition' conditionNames{UniCondi(u)}]);
    title(['condition: ' num2str(conditionNames(UniCondi(u))) ' (radian)']);
  
end


%% take a look at the EEG
% plotx(data_trials(:,1:128,2));
% EEG_trial1=squeeze(data_trials(:,1:128,50));
% EEG_trial1mean=mean(EEG_trial1,2).*ones(1,128);
% EEG_trial1=EEG_trial1-EEG_trial1mean;
% plotx(EEG_trial1);

%% clean artifact
filtered_data=data_trials(:,1:128,:);

threshold=60; % set a amplitude threshold as criteria to remove bad channels and epochs later

% set good standard percentages
prctepoch = 70; % percentage of epochs that should be good
prctchan = 70; % percentage of channels that should be good

% set bad standard percentage
cirtepoch =  1-prctepoch/100; % the percentage of bad epochs allowed
critchan = 1-prctchan/100; % the percentage of bad channels allowed

nSamps = size(filtered_data,1);
nChans = size(filtered_data,2);
nEpochs = size(filtered_data,3);

 % find the the empty channels
sumovertime=sum(filtered_data,1);
sumovertrials=sum(sumovertime,3);
zerochans = find(sumovertrials== 0);

% design a matrix with highest value for each channel on each trial(each single channel in each epoch). 
trialdatamax = squeeze(max(abs(filtered_data),[],1));  
% examine the maximal abs trial matrix
% imagesc(trialdatamax);colorbar;ylabel('chan');xlabel('trial');

% create a bad channel list containing the last 4 channels to be removed
badchans = [129 130 131 132 133 134]; 

% initiate a good channel list by removing the badchans
goodchans = setdiff(1:nChans,badchans); 

% create the bad trials matrix
trialpass = zeros(size(trialdatamax));
trialpass(find(trialdatamax > threshold)) = 1;   % locates the 1 (bad), i.e. the bad trials
% examine the bad trials
% imagesc(trialpass);colorbar;

% find the bad epochs/trials
trialpass = sum(trialpass(goodchans,:),1)/length(goodchans); % the percentage of bad channels on each epoch
badepochs = find(trialpass > cirtepoch);  %determines if an epoch has too many of the bad channels. 

% create a good epoch list 
goodepochs = setdiff(1:nEpochs,badepochs);

% create the bad trials matrix again
chanpass = zeros(size(trialdatamax)); % create the bad trial matrix
chanpass(find(trialdatamax >threshold)) = 1; %locates the 1 (bad), i.e.bad trials
% examine the bad channel matrix
% imagesc(chanpass);colorbar;

% find the bad channels
chanpass = sum(chanpass(:,goodepochs),2)/length(goodepochs); % the percentage of bad epochs in each channel
newbadchan = goodchans(find(chanpass(goodchans) > critchan)); %  a new list of bad channels

badchans = [badchans newbadchan];% update the list of bad channels 
goodchans = setdiff(1:nChans,badchans); % update the list of good channels

artifact.goodchans = goodchans;
artifact.goodepochs = goodepochs;
artifact.editmatrix = trialdatamax;
artifact.threshhold = threshold;
artifact.zerochans = zerochans;


%% Re reference
% This function pass filtered_data and get back average referenced data.
[reRef_data] = reRef(filtered_data,goodchans);
% erpdata = mean(reRef_data(:,:,goodepochs),3); % This includes only good epochs
% plotx(erpdata(:,goodchans));% plot only the good channels
% xlabel(['time (ms)']);ylabel(['amplitude']);title(['ERP'])

%% power spetra and coorelation
% fs=2000; maxfreq=50;
% [pow,freqs,df,eppow,corr,cprod,fcoef] = allspectra(reRef_data,fs,maxfreq,goodepochs);


%% using laplacian_perrinX
load('chanlocs.mat')
chanlocs=chanlocs(1:128); 

X = [chanlocs.X];
Y = [chanlocs.Y];
Z = [chanlocs.Z];

plot3(X,Y,Z)

point_id={chanlocs.labels};point_id=point_id(1:128);
scatter3(X,Y,Z, 'b*');
text(X,Y,Z, point_id);

% [X,Y] = meshgrid(1:0.5:10,1:20);
% Z = sin(X) + cos(Y);
% surf(X,Y,Z)
% surf(X,Y,Z)

% data to use is all channels in trial 1
data1=reRef_data(2000,1:128,1);
data1tranpose=data1';

% apply laplacian_perrinx
n=10;
[surf_lap,G,H] = laplacian_perrinX(data1tranpose,X,Y,Z,n,1e-5);

% plot the 1000th time point
subplot(1,2,1);
topoplot(data1tranpose,chanlocs);
title('Original Data');
colorbar;
% set(gca,'clim',[-150 150]);

subplot(1,2,2);
topoplot(surf_lap,chanlocs);
title('Surface Laplacian');
colorbar;
% set(gca,'clim',[-150 150]);

% https://mikexcohen.com/lecturelets/
% http://mikexcohen.com/lecturelets/laplacian/laplacian.html
%%
% plot the 1000th time point
subplot(1,2,1);
topoplot(data1tranpose,chanlocs);
title('Original Data');
colorbar;
set(gca,'clim',[-150 150]);

subplot(1,2,2);
topoplot(surf_lap,chanlocs);
title('Surface Laplacian');
colorbar;
set(gca,'clim',[-150 150]);
%% GIF of Xcorr from trial 1 to 50 

% run plotXcorr.m instead

h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
titlename='Xcoor-theta';
filename = [titlename '.gif'];

% conditionNames={'A1&A2','A1&A1','A2&A1','A2&A2'};
% allPerm
% Perm=[2 3 1 4];
% allPerm
% conditionNames{2}

% for i=1:50 % trial number
% for i=51:100 
% for i=101:150
for i=151:200
    C3EEG=reRef_data(:,28,i);
    C4EEG=reRef_data(:,32,i);

    [r,lags]=xcorr(C3EEG, C4EEG);
    plot(lags./2,r);xlabel('time [ms]');ylabel('cross correlation');xlim([-250 250]);%ylim([-2*1e8 2*1e8]);
    title(titlename);

    drawnow 
    % Capture the plot as an image 
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if i == 1 %51 %101 %151
      imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
    else 
      imwrite(imind,cm,filename,'gif','WriteMode','append'); 
    end 

end

% average 
rAll=[];

% for i=1:50 % trial number
% for i=51:100 
% for i=101:150
for i=151:200
    C3EEG=reRef_data(:,28,i);
    C4EEG=reRef_data(:,32,i);

    [r,lags]=xcorr(C3EEG, C4EEG);
    rAll=[rAll r];
 
end
rAllmean=mean(rAll,2);
plot(lags./2,rAllmean);xlabel('time [ms]');ylabel('cross correlation');xlim([-250 250]);%ylim([-2*1e8 2*1e8]);
title(['mean ' titlename]);



