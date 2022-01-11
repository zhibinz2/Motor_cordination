%% Check EEG on 2021018
% load EEG all as one
% cd /home/zhibin/Documents/Acquisition/bimanual_Reach_zhibin_20211018

% filenames=dir; EEGfileName=filenames(3).name;
% [EEGfileName]=uigetfile('*.cdt');

% cd /home/zhibin/Documents/Acquisition/
% EEG=loadcurry([pwd '/' 'bimanual_Reach_Jack_20211118.cdt']);

%% Get useful variables from EEG
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

% method 1
datatimes_ind_trial_go=latency(find(type==Unitype(3)));
datatimes_ind_trial_start=datatimes_ind_trial_go-1000;
datatimes_ind_trial_end=datatimes_ind_trial_go+1000;

%     % method 2
    datatimes_ind_trial_startend=latency(find(type==Unitype(4)));
    datatimes_ind_trial_start=datatimes_ind_trial_startend(1:2:end-1);
    datatimes_ind_trial_end=datatimes_ind_trial_startend(2:2:end);

    % examine the start and end points
    figure;
    subplot(1,2,1);
    plot(datatimes(datatimes_ind_trial_end)-datatimes(datatimes_ind_trial_start),'ro');
    xlabel('each trial');ylabel('length of time (ms)');
    subplot(1,2,2);
    hist(datatimes(datatimes_ind_trial_end)-datatimes(datatimes_ind_trial_start));
    ylabel('number of trials');xlabel('length of time (ms)');


%% photocell processing to find trial start and end index (too complicative and unreliable,use event instead)
% 
%     % plot photocell analog signal
%     analog1=133; 
%     analog1data=dataEEG(analog1,:);
% 
%     % Examine analog1 and analog2
%     close;figure;
%     plot(datatimes,analog1data,'b'); % might need to mannually cut off the start and end of block
%     hold on;
%     analog2=134; 
%     analog2data=dataEEG(analog2,:);
%     plot(datatimes,analog2data,'r'); % weird start and end of block signal
%     hold off;
% 
%     % Set start and End of data
%     % Startpoint=50000;
%     % Endpoint=2340000;
%     % analog1data=analog1data(Startpoint,Endpoint);
%     % datatimes=datatimes(Start,Endpoint);
% 
%     %***************************** Use the analog1data 
%     
%     % Examine peaks detection in analog1
%     Halfhigh1=3/4*(max(analog1data)-min(analog1data)); % value for 'MinPeakProminence'
%     % Check if need to adjust the Halfhigh cutoff
%     close;figure;
%     findpeaks(analog1data,datatimes,'MinPeakProminence',Halfhigh1,'Annotate','extents');
%     yline(Halfhigh1,'m','MinPeakProminence');
% 
%     % locate the trial sessions % pks=value of the peak % locs=time of the peak
%     [pks,locs] = findpeaks(analog1data,datatimes,'MinPeakProminence',Halfhigh1,'Annotate','extents');
%     
%     % examine pks and locs (both are values of analog1data and datatimes, not indices)
%     i=1;
%     find(analog1data==pks(i)) % return the index of peak in analog1data
%     find(datatimes==locs(i)) % return the same "index" in datatimes (ie, i=1, index=4226)
%     % so, infer from above, beacause the above "find values" for the same "index" works
%     int64(pks(i)) % return the value in analog1data, pks are values in analog1data
%     int64(locs(i)) % return the value in datatimes, locs are values in datatimes ((ie, i=1, value=4225))
% 
%     figure;plot(locs,'bo');
%     figure;bar(locs);
%     figure;plot(pks,'bo');
%     figure;bar(pks);
%     figure;plot(locs,pks,'bo');ylim([min(analog1data) max(analog1data)]);
% 
%     locsDurations=diff(locs);% time durations between peaks
%     close;figure;
%     plot(locsDurations,'ro');% look at distribution of these durations   
%     xlabel('each peak');ylabel('locsDurations = Between-peaks time duration (ms)');
%     Lowcutoff=5*mean(locsDurations);% cutoff standard of the between-peak duarations to separate between trials
%     Highcutoff=180*mean(locsDurations);
%     hold on; yline(Lowcutoff,'m--','lowcut');hold off; % examine the cutoff line
%     hold on; yline(Highcutoff,'m--','highcut');hold off;
%     
%     % find indices of these large between-pks duartions at the end of each trial (index in locsDurations) 
%     locsDurations_ind_locsDurations=find([locsDurations, inf]>Lowcutoff); % add a number of "inf" at the end to include the last trial
%     % add "Highcutoff" at the end to include the last trial and remove artifact of the first peak
%     locsDurations_ind_locsDurations=find([locsDurations, Highcutoff]>Lowcutoff & [locsDurations, Highcutoff]<=Highcutoff); 
%     
%     % number of peaks in each trial
%     trialPeaks=diff([0 locsDurations_ind_locsDurations]); 
%     trialPeaks=diff([1 locsDurations_ind_locsDurations]); % if 1 peak should be omitted before the first trial start
%     plot(trialPeaks,'ro'); % look at trial length (number of peaks in each trial)
%     
%     % index start and end of trial in locs  (index in locs)
%     locs_ind_trial_end=cumsum(trialPeaks); % index of the last peak in each trial (index in locs)
%     locs_ind_trial_end=cumsum(trialPeaks)+1; % if 1 peak should be omitted before the first trial start
%     locs_ind_trial_start=locs_ind_trial_end-trialPeaks+1; % index of the first peak in each trial (index in locs)
%     
%     % use the indices of start and end of trial in locs to get values in datatimes
%     clear datatimes_ind_trial_end datatimes_ind_trial_start
%     datatimes_ind_trial_end=locs(locs_ind_trial_end)+1; % plus one to get indices in datatimes, because datatimes starts from 0
%     datatimes_ind_trial_start=locs(locs_ind_trial_start)+1;
% 
%     % find end point of data  (index in datatimes)
%     datatimes_dataEnd=locs(locs_ind_trial_end(end));
%     
% 
%     %*****************************  checking starting peak of the trials (index in locs)
%     close;figure;
%     subplot(1,2,1);
%     plot(locs_ind_trial_start,zeros(length(ind_trial_start)),'go'); 
%     xlabel('indices in the peak locations (locs)')
%     hold on;
%     plot(locs_ind_trial_end,ones(length(ind_trial_end)),'ro');
%     ylim([-20 20]);
%     subplot(1,2,2);
%     plot(locs_ind_trial_end-locs_ind_trial_start+1,'ro');
%     ylabel('number of peaks');xlabel('each trial');
%     
%     % checking duration of each trial (index in datatimes)
%     close;figure;
%     subplot(1,2,1);
%     plot(datatimes_ind_trial_start,zeros(length(datatimes_ind_trial_start)),'go'); 
%     xlabel('datatimes (not integers?)');
%     hold on;
%     plot(datatimes_ind_trial_end,ones(length(datatimes_ind_trial_end)),'ro');
%     ylim([-20 20]);
%     subplot(1,2,2);
%     plot(datatimes_ind_trial_end-datatimes_ind_trial_start,'go'); % look at trial duration in ms
%     yline(mean(datatimes_ind_trial_end-datatimes_ind_trial_start),'m--','mean');
%     yline(median(datatimes_ind_trial_end-datatimes_ind_trial_start),'m--','median');
%     % or
%     plot(datatimes(int64(datatimes_ind_trial_end))-datatimes(int64(datatimes_ind_trial_start)),'go'); % look at trial duration in ms
%     ylabel('time(ms)');xlabel('each trial');
%     
%     % try out double and integer index ##############
%     num2str(datatimes_ind_trial_start(1:50))
%     datatimes([1 3 5])
%     
%     datatimes([1.00 3.00 5.00])
%     A=double([1.00 3.00 5.00]);
%     datatimes(A)
%     
%     datatimes(int64(datatimes_ind_trial_end))
%     
%     X = 10+randn(16,16);
%     X = uint8(X);
%     
%     x_double = [1 2 3];
%     x_int = int64(x_double);
%     
%     a = randi(100, 4);
%     a = int64(a);
% 
% 
%     %***************************** Use the analog2data
%     
%     % Examine peaks detection in analog2
%     Halfhigh2=mean(analog2data); % value for 'MinPeakProminence'
%     % Check if need to adjust the Halfhigh cutoff
%     close;figure;
%     findpeaks(analog2data,datatimes,'MinPeakProminence',Halfhigh2,'Annotate','extents');
%     yline(Halfhigh2,'m');
%     
%     % elect from the plot above
%     cut1=47492;
%     cut2=2351708;
%     hold on; xline(cut1,'m','cut1');xline(cut2,'m','cut2');hold off;
%     
%     % locate the trial sessions % pks=value of the peak % locs=time of the peak
%     [pks,locs] = findpeaks(analog2data(cut1:cut2),datatimes(cut1:cut2),'MinPeakProminence',Halfhigh2,'Annotate','extents');
%     
%     locsDurations=diff(locs);% time duarations between peaks
%     close;figure;
%     plot(locsDurations,'go'); ylabel('locsDurations'); % look at distribution of these distances
%     cutoff=mean(locsDurations)+mean(locsDurations)/5; % cutoff standard of the peak distance between trials
%     hold on; yline(cutoff,'m--');hold off; % examine the cutoff line
%     locsDurations_ind_locsDurations=find([locsDurations, inf]>cutoff); % indices of pks distances at the end of each trial 
%     trialPeaks=diff([0 locsDurations_ind_locsDurations]); % number of peaks in each trial
%     close;figure;
%     plot(trialPeaks,'ro'); % look at trial length (number of peaks in each trial)
%     
%     % index start and end of trial in locs  (index in locs)
%     locs_ind_trial_end=cumsum(trialPeaks); % index of the last peak in each trial (index in locs)
%     locs_ind_trial_start=locs_ind_trial_end-trialPeaks+1; % index of the first peak in each trial (index in locs)
%     
%     % the indices of start and end of trial in datatimes  
%     % method 1
%     datatimes_ind_trial_end=int64(locs(locs_ind_trial_end)); % (locs(locs_ind_trial_end) are values in datatimes(cut1:cut2),not indices)
%     datatimes_ind_trial_start=int64(locs(locs_ind_trial_start));
%     datatimes_ind_trial_start(end)
%     int64(datatimes_ind_trial_start(end))
%     % or method 2 
%     datatimes_ind_trial_end=int64(datatimes(int64(locs(locs_ind_trial_end))+1)); % plus one, because datatimes starts from 0
%     datatimes_ind_trial_start=datatimes(int64(locs(locs_ind_trial_start))+1); 
%     datatimes_ind_trial_start(end)
%     int64(datatimes_ind_trial_start(end))
%     
%     %*****************************  checking starting peak of the trials (index in locs)
%     close;figure;
%     subplot(1,2,1);
%     plot(locs_ind_trial_start,zeros(length(ind_trial_start)),'go'); 
%     xlabel('indices in the peak locations (locs)')
%     hold on;
%     plot(locs_ind_trial_end,ones(length(ind_trial_end)),'ro');
%     ylim([-20 20]);
%     subplot(1,2,2);
%     plot(locs_ind_trial_end-locs_ind_trial_start+1,'ro');
%     ylabel('number of peaks');xlabel('each trial');
%     
%     % checking duration of each trial (index in datatimes)
%     close;figure;
%     subplot(1,2,1);
%     plot(datatimes_ind_trial_start,zeros(length(datatimes_ind_trial_start)),'go'); 
%     xlabel('datatimes (not integers?)');
%     hold on;
%     plot(datatimes_ind_trial_end,ones(length(datatimes_ind_trial_end)),'ro');
%     ylim([-20 20]);
%     subplot(1,2,2);
%     plot(datatimes_ind_trial_end-datatimes_ind_trial_start,'go'); % look at trial duration in ms
%     yline(mean(datatimes_ind_trial_end-datatimes_ind_trial_start),'m--','mean');
%     yline(median(datatimes_ind_trial_end-datatimes_ind_trial_start),'m--','median');
%     % or
%     plot(datatimes(int64(datatimes_ind_trial_end))-datatimes(int64(datatimes_ind_trial_start)),'go'); % look at trial duration in ms
%     ylabel('time(ms)');xlabel('each trial');

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