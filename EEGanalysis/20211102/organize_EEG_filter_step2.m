%% Check EEG on 2021018
% load EEG all as one
% cd /home/zhibin/Documents/Acquisition/bimanual_Reach_zhibin_20211018

% filenames=dir; EEGfileName=filenames(3).name;
% [EEGfileName]=uigetfile('*.cdt');

% EEG=loadcurry([pwd '/' 'bimanual_Reach_Jack_20211118.cdt']);

%%
dataEEG=EEG.data;
datatimes=EEG.times;
Fs=EEG.srate;
% plot(datatimes,dataEEG);
% plot(datatimes,dataEEG(1:128,:));

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
Halfhigh1=3/4*(max(analog1data)-min(analog1data));
% Halfhigh2=1/2*(max(analog2data)-min(analog2data));

% Check if need to adjust the Halfhigh cutoff
findpeaks(analog1data,datatimes,'MinPeakProminence',Halfhigh1,'Annotate','extents');
yline(Halfhigh1);


%***************************** Use the analog1data 
% locate the trial sessions % pks=value of the peak % locs=time of the peak
[pks,locs] = findpeaks(analog1data,datatimes,'MinPeakProminence',Halfhigh1,'Annotate','extents');

locsDiff=diff(locs);% time distances between peaks
% plot(locsDiff,'ro');% look at distribution of these distances
cutoff=5*mean(locsDiff);% cutoff standard of the peak distance between trials
% hold on; yline(cutoff,'m--');hold off; % examine the cutoff line
ind_locsDiff=find([locsDiff, inf]>cutoff); % pks distance indices at the end of each trial 
trialLength=diff([0 ind_locsDiff]); % number of peaks in each trial
% plot(trialLength,'ro'); % look at trial length (number of peaks in each trial)
ind_trial_end=cumsum(trialLength); % index of the last peak in each trial (index in locs)
ind_trial_start=ind_trial_end-trialLength+1; % index of the first peak in each trial (index in locs)

% find end point of data
dataEnd=locs(ind_trial_end(end));

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
% % plot(locsDiff,'go'); ylabel('locsDiff'); % look at distribution of these distances
% cutoff=mean(locsDiff)+mean(locsDiff)/5; % cutoff standard of the peak distance between trials
% hold on; yline(cutoff,'m--');hold off; % examine the cutoff line
% ind_locsDiff=find([locsDiff, inf]>cutoff); % indices of pks distances at the end of each trial 
% trialLength=diff([0 ind_locsDiff]); % number of peaks in each trial
% % plot(trialLength,'ro'); % look at trial length (number of peaks in each trial)
% ind_trial_end=cumsum(trialLength); % index of the last peak in each trial (index in locs)
% ind_trial_start=ind_trial_end-trialLength+1; % index of the first peak in each trial (index in locs)


%% detrend the data

addpath(genpath('/home/zhibin/Documents/GitHub'));
% data=dataEEG(1:128,2000000:2002000);
EEGdata=dataEEG(1:128,:);

% tranpose the data so that time is the first dimension, channel is 2nd dimension, and trial is the 3rd dimension
transpose_data=permute(EEGdata,[2,1]);
% plotx(mean(transpose_data,3));
% plotx(transpose_data);hold on;xline(dataEnd,'m--');

% detrend the data
detrend_data=ndetrend(transpose_data,1); 
% plotx(mean(detrend_data,2)); xlim([5000 10000]);
% plotx(detrend_data); xlim([5000 10000]);
% plotx(detrend_data(5000:10000,:));

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
% plotx(filtered_data);hold on;xline(dataEnd,'m--');
% plotx(filtered_data(5000:10000,:));

%%
filtered_broadband=filtered_data;
% filtered_data=filtered_broadband;

%% fiter5bands 
% tic;
% run filter5bands.m
% toc;

% filtered_data=filtered_delta;
% filtered_data=filtered_theta;
% filtered_data=filtered_alpha;
% filtered_data=filtered_beta;
% filtered_data=filtered_gamma;

% plot(filtered_beta(5000:10000,8));
% figure;
% plot(filtered_gamma(5000:10000,8));

% tic;
% % cd /home/zhibin/Documents/Acquisition/Bimanual_reach_zhibin_20211102
% filename='filtered_data.mat';
% % save(filename); % save all variables in the workspace
% save(filename, 'filtered_broadband', 'filtered_delta', 'filtered_theta', 'filtered_alpha', 'filtered_beta', 'filtered_gamma', '-v7.3');
% toc;