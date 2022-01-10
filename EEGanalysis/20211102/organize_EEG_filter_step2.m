%% Check EEG on 2021018
% load EEG all as one
% cd /home/zhibin/Documents/Acquisition/bimanual_Reach_zhibin_20211018

% filenames=dir; EEGfileName=filenames(3).name;
% [EEGfileName]=uigetfile('*.cdt');

% cd /home/zhibin/Documents/Acquisition/
% EEG=loadcurry([pwd '/' 'bimanual_Reach_Jack_20211118.cdt']);

%%
dataEEG=EEG.data;
datatimes=EEG.times;
Fs=EEG.srate;
event=EEG.event;

% plot(datatimes,dataEEG);
% plot(datatimes,dataEEG(1:128,:));

%% using EEG.events to find tiral start and end index
type=[event.type];
Unitype=unique(type);
latency=[event.latency];

% examine the events
plot(latency,type,'ro');xlabel('time(ms)');ylabel('event type');

ind_trial_start=latency(find(type==Unitype(3)))-500;
ind_trial_end=latency(find(type==Unitype(3)))+1000;


%% photocell processing to find trial start and end index (too complicative and unreliable,use event instead)

    % plot photocell analog signal
    analog1=133; 
    analog1data=dataEEG(analog1,:);

    % Examine analog1 and analog2
    close;figure;
    plot(datatimes,analog1data,'b'); % might need to mannually cut off the start and end of block
    hold on;
    analog2=134; 
    analog2data=dataEEG(analog2,:);
    plot(datatimes,analog2data,'r'); % weird start and end of block signal
    hold off;

    % Set start and End of data
    % Startpoint=50000;
    % Endpoint=2340000;
    % analog1data=analog1data(Startpoint,Endpoint);
    % datatimes=datatimes(Start,Endpoint);

    %***************************** Use the analog1data 
    
    % Examine peaks detection in analog1
    Halfhigh1=3/4*(max(analog1data)-min(analog1data)); % value for 'MinPeakProminence'
    % Check if need to adjust the Halfhigh cutoff
    close;figure;
    findpeaks(analog1data,datatimes,'MinPeakProminence',Halfhigh1,'Annotate','extents');
    yline(Halfhigh1,'m');

    % locate the trial sessions % pks=value of the peak % locs=time of the peak
    [pks,locs] = findpeaks(analog1data,datatimes,'MinPeakProminence',Halfhigh1,'Annotate','extents');

    locsDurations=diff(locs);% time durations between peaks
    close;figure;
    plot(locsDurations,'ro');% look at distribution of these durations   
    xlabel('each peak');ylabel('locsDurations = Between-peaks time duration (ms)');
    Lowcutoff=5*mean(locsDurations);% cutoff standard of the between-peak duarations to separate between trials
    Highcutoff=180*mean(locsDiff);
    hold on; yline(Lowcutoff,'m--','lowcut');hold off; % examine the cutoff line
    hold on; yline(Highcutoff,'m--','highcut');hold off;
    
    % find indices of these large between-pks duartions at the end of each trial (index in locsDurations) 
    ind_locsDiff=find([locsDurations, inf]>Lowcutoff); % add a number of "inf" at the end to include the last trial
    % add "Highcutoff" at the end to include the last trial and remove artifact of the first peak
    ind_locsDiff=find([locsDurations, Highcutoff]>Lowcutoff & [locsDurations, Highcutoff]<=Highcutoff); 
    
    % number of peaks in each trial
    trialPeaks=diff([0 ind_locsDiff]); 
    
    % plot(trialPeaks,'ro'); % look at trial length (number of peaks in each trial)
    ind_trial_end=cumsum(trialPeaks); % index of the last peak in each trial (index in locs)
    ind_trial_start=ind_trial_end-trialPeaks+1; % index of the first peak in each trial (index in locs)

    % find end point of data
    dataEnd=locs(ind_trial_end(end));

    % % checking starting peak of the trials (index in locs)
    % plot(ind_trial_start,zeros(length(ind_trial_start)),'go'); 
    % hold on;
    % plot(ind_trial_end,ones(length(ind_trial_end)),'ro');
    % ylim([-20 20]);

    %***************************** Use the analog2data 
    
    % Examine peaks detection in analog2
    Halfhigh2=prctile(analog2data,95); % value for 'MinPeakProminence'
    % Check if need to adjust the Halfhigh cutoff
    close;figure;
    findpeaks(analog2data,datatimes,'MinPeakProminence',Halfhigh2,'Annotate','extents');
    yline(Halfhigh2,'m');
    
    % locate the trial sessions % pks=value of the peak % locs=time of the peak
    [pks,locs] = findpeaks(analog2data,datatimes,'MinPeakProminence',Halfhigh2,'Annotate','extents');
    
    locsDiff=diff(locs);% time distances between peaks
    % plot(locsDiff,'go'); ylabel('locsDiff'); % look at distribution of these distances
    cutoff=mean(locsDiff)+mean(locsDiff)/5; % cutoff standard of the peak distance between trials
    hold on; yline(cutoff,'m--');hold off; % examine the cutoff line
    ind_locsDiff=find([locsDiff, inf]>cutoff); % indices of pks distances at the end of each trial 
    trialLength=diff([0 ind_locsDiff]); % number of peaks in each trial
    % plot(trialLength,'ro'); % look at trial length (number of peaks in each trial)
    ind_trial_end=cumsum(trialLength); % index of the last peak in each trial (index in locs)
    ind_trial_start=ind_trial_end-trialLength+1; % index of the first peak in each trial (index in locs)


%% This could be used as one event channel to import into EEGLAB
%     EventTimes=locs(ind_trial_end);% 
%     for i=1:length(EventTimes)
%         EventInd(i)=find(datatimes==EventTimes(i));
%     end
%     EventChannel=zeros(1,length(datatimes));
%     EventChannel(EventInd-1500)=1;% trajectory appear onset

%% select the EEG channels

addpath(genpath('/home/zhibin/Documents/GitHub'));
% data=dataEEG(1:128,2000000:2002000);
EEGdata=dataEEG(1:128,:);

%% Appen EventChannel for EEGLAB
% EEGdata(129,:)=EventChannel;

%% detrend the data

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
Hd = makefilter(Fs,1.5,1,3,20,0); % Astop 6 or 20;
filtered_data=filtfilthd(Hd,double(detrend_data));

% remove channel 91 44
% filtered_data(:,91)=[];filtered_data(:,44)=[];
%  plotx(filtered_data); xlim([5000 10000]);

% apply low pass filter to filter out frequencies higher than 50;
Hd = makefilter(Fs,50,55,3,20,0);
filtered_data=filtfilthd(Hd,filtered_data);

% plotx(mean(filtered_data,2));  
% plot(filtered_data(1:dataEnd,:));
% hold on;
% for i=1:length(ind_trial_start)
%     xline(locs(ind_trial_start(i)),'k');
% end
% xline(dataEnd,'m--');xlim([0 dataEnd]);
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