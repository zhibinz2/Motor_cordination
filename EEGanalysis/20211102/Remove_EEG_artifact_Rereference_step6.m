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
