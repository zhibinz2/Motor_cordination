%% organize filtered_data into data_trials=sample_timepoints x (139+4 channels) x  trial  
% Fs=1000;ifi=0.0167;
NumTrialsRecorded=length(datatimes_ind_trial_start);
% NumTrialtimepoints=round(max(trialPeaks)*2*ifi*Fs);% number of samples in a trial,ideally
NumTrialtimepoints=round(mean((datatimes(datatimes_ind_trial_end)-datatimes(datatimes_ind_trial_start))));
NumEEGChannels=size(detrend_data,2);
NumAllChannels=NumEEGChannels+6; % +4;

% Initialize the data_trials matrix
% Zeros at the end of the trial session if trials shorter than is should be, ideally
% Reality is, sometimes this trial length do not cover the length of some trial
% data_trials=zeros(500+NumTrialtimepoints,NumAllChannels,NumTrialsRecorded);% add the 500 ms before trial as baseline
data_trials=zeros(NumTrialtimepoints,NumAllChannels,NumTrialsRecorded);
% (You need to go 500 ms before the joint dot turned green up to 1.5 secs after the Go signal.
% I like to have a little padding because of the way the spectrogram works.
% IF you have this, we can compute a spectrogram.)
% I also add the 1 second of bonus display after movement
% So the 4000 timepoints are structured as 
% "500ms padding + 500ms green dot baseline + 500ms plan phase + 1000ms movement + 1000ms show bonus + 500ms padding"


% Initialize the indices for trial start and end time point
trialstart=zeros(1,NumTrialsRecorded);
trialend=zeros(1,NumTrialsRecorded);

%     % Method 1 (use locs_ind_trial_end)
%     for ntr=1:NumTrialsRecorded % take too much memory, matlab might quit
%         trialend(ntr)  =find(datatimes==locs(locs_ind_trial_end(ntr)));
%         data_trials((1:(500+NumTrialtimepoints)),1:NumEEGChannels,ntr)=filtered_data((trialend(ntr)-(500+NumTrialtimepoints)+1):trialend(ntr),1:NumEEGChannels); 
%     end
%     % Keep in mind that after the filter, the last 10 seconds or so of the data is not usable)
% 
%     % Method 2 -- timing(trial length) might mess up
%     for ntr=1:NumTrialsRecorded % take too much memory, matlab might quit
%         trialstart(ntr)=find(datatimes==locs(locs_ind_trial_start(ntr)));
%         trialend(ntr)  =find(datatimes==locs(locs_ind_trial_end(ntr)));
%         data_trials(1:(trialend(ntr)-trialstart(ntr)+1),1:NumEEGChannels,ntr)=filtered_data(trialstart(ntr):trialend(ntr),1:NumEEGChannels); 
%     end
%     % Keep in mind that after the filter, the last 10 seconds or so of the data is not usable)

%     % Method 3 -- use datatimes_ind_trial_start and datatimes_ind_trial_end
%     for ntr=1:NumTrialsRecorded 
%         data_trials(1:(datatimes_ind_trial_end(ntr)-datatimes_ind_trial_start(ntr)+1),1:NumEEGChannels,ntr)=...
%             filtered_data(datatimes_ind_trial_start(ntr):datatimes_ind_trial_end(ntr),1:NumEEGChannels); 
%     end
%     % examine
%     data_trials;

% Method 4 -- use datatimes_ind_trial_go
for ntr=1:NumTrialsRecorded 
    data_trials(1:NumTrialtimepoints,1:NumEEGChannels,ntr)=...
        filtered_data(datatimes_ind_trial_go(ntr)-1000:datatimes_ind_trial_go(ntr)+1000-1,1:NumEEGChannels); 
end
% examine
data_trials;


% examine 
% plot(ind_trial_end-ind_trial_start,'ro');ylabel('trial peaks');xlabel('trial');
% 
% plot(datatimes_ind_trial_end-datatimes_ind_trial_start,'ro');ylabel('trial peaks');xlabel('trial');

% locsstart=locs(ind_trial_start);
% locsend=locs(ind_trial_end);
% plot(trialend-trialstart,'ro');ylabel('trial length');xlabel('trial');
% plot(trialstart,ones(1,length(trialstart)),'go');hold on;plot(trialend,ones(1,length(trialend)),'ro'); % check 1 second inter-trial interval
% plot(data_trials(:,1:128,5));
% plot(datatimes);
