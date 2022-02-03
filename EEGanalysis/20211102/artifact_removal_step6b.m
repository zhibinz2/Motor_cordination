%% take a look at the EEG
% figure; plotx(data_trials(:,1:128,2));
% EEG_trial1=squeeze(data_trials(:,1:128,50));
% EEG_trial1mean=mean(EEG_trial1,2).*ones(1,128);
% EEG_trial1=EEG_trial1-EEG_trial1mean;
% figure; plotx(EEG_trial1);

%% clean artifact
data_trials_EEG=data_trials(:,1:128,:);

nSamps = size(data_trials_EEG,1);
nChans = size(data_trials_EEG,2);
nEpochs = size(data_trials_EEG,3);

% find the the empty channels
sumovertime=sum(data_trials_EEG,1);
sumovertrials=sum(sumovertime,3);
zerochans = find(sumovertrials== 0);

% design a matrix with highest value for each channel on each trial(each single channel in each epoch). 
trialdatamax = squeeze(max(abs(data_trials_EEG),[],1));  
% examine the maximal abs trial matrix
% figure;imagesc(trialdatamax);colorbar;ylabel('chan');xlabel('trial');
% figure;imagesc(log10(trialdatamax));colorbar;

% first round of artifact labelling **************#########################
threshold=500; % set a amplitude threshold as criteria to remove bad channels and epochs later

% Create a label matrix that labels bad chan and epoch as 1 
binarythreshhold = zeros(size(trialdatamax));
binarythreshhold(trialdatamax > threshold) = 1;
% figure; imagesc(binarythreshhold); colorbar;title(['sum of Ones: ' num2str(sum(binarythreshhold,'all'))]);

chansum = sum(binarythreshhold,2); % how many bad epochs in each chan
epochsum = sum(binarythreshhold,1); % how many bad chans in each epoch

% bar(chansum); % examine the chans 
badchancrit = prctile(chansum,0.001); % criteria to determine bad chans
badchanedit = find(chansum > badchancrit); % get the initial bad chans
goodchanedit = setdiff([1:nChans],badchanedit); % get the initial good chans

% bar(epochsum); % examine the epochs
badepochcrit = prctile(epochsum,0.001); % criteria to determine bad epochs
badepochedit = find(epochsum > badepochcrit); % get the initial bad epochs
goodepochedit = setdiff([1:nEpochs],badepochedit); % get the initial good epochs

% second round of labelling **************#################################
threshold2=100;
% update the label matrix that labels bad chan and epoch as 1 
binarythreshhold = zeros(size(trialdatamax));
binarythreshhold(trialdatamax > threshold2) = 1;
% figure; imagesc(binarythreshhold); colorbar;
% figure; imagesc(binarythreshhold(goodchanedit,goodepochedit));colorbar;

epochsum2 = sum(binarythreshhold(goodchanedit,goodepochedit),1); % how many bad chans in each good epoch of the upated label matrix % whos % to look at the workspace
% figure;bar(goodepochedit,epochsum2);
newbadepochcrit = prctile(epochsum2,90); % criteria to determine bad epochs in 2nd round
newbadepochedit = goodepochedit(find(epochsum2 > newbadepochcrit)); % get the bad epochs in 2nd round
badepochedit = [badepochedit newbadepochedit]; % add the bad epochs in 1s round
goodepochedit_final = setdiff([1:nEpochs],badepochedit);  % get the final good epochs

chansum2 = sum(binarythreshhold(goodchanedit,goodepochedit),2); % how many bad epochs in each good chan of the upated label matrix
% figure;bar(goodchanedit,chansum2);
newbadchancrit = prctile(chansum2, 90); % criteria to determine bad chans in 2nd round
newbadchanedit = goodchanedit(find(chansum2 > newbadchancrit)); % get the bad chans in 2nd round
badchanedit = [badchanedit; newbadchanedit']; % add the bad chans in 1s round
goodchanedit_final = setdiff([1:nChans],badchanedit);  % get the final good chans

% examine the updated label matrix after the 2nd round
% figure;imagesc(trialdatamax(goodchanedit_final,goodepochedit_final));colorbar;
% figure;imagesc(log10(trialdatamax(goodchanedit_final,goodepochedit_final)));colorbar;
% compare with previous one round artifact labelling method
% figure;imagesc(log10(trialdatamax(goodchans,goodepochs)));colorbar;
% figure;imagesc(log10(trialdatamax));colorbar;

% create a bad channel list containing the last 4 channels to be removed (skip if not including behaviral data)
% badchans = [129 130 131 132 133 134]; 
% add to the final list of bad chans
% badchans = [badchans badchanedit'];% update the list of bad channels 

badchans = badchanedit;
badepochs = badepochedit;
goodchans = goodchanedit_final;
goodepochs = goodepochedit_final;

artifact.badchans = badchans;
artifact.badepochs = badepochs;
artifact.goodchans = goodchans;
artifact.goodepochs = goodepochs;
artifact.editmatrix = trialdatamax;
artifact.threshhold = threshold2;
artifact.zerochans = zerochans;

%% goodepochs in boolean
Bgoodepochs=zeros(1,nEpochs);
Bgoodepochs(goodepochs)=1;

%% Re reference
% This function pass filtered_data and get back average referenced data.
[reRef_data] = reRef(data_trials_EEG,goodchans);

% figure; % to examine reRef_data
% erpdata = mean(reRef_data(:,:,goodepochs),3); % This includes only good epochs
% plotx(erpdata(:,goodchans));% plot only the good channels
% xlabel(['time (ms)']);ylabel(['amplitude']);title(['ERP'])


% figure; % to compare with data before re reference
% erpdata = mean(data_trials_EEG(:,:,goodepochs),3); % This includes only good epochs
% plotx(erpdata(:,goodchans));% plot only the good channels

%% Remove bad epochs and connect them back to one time serie and plot only the good channels
% reRef_data_goodepochs=[];
% for i=1:length(goodepochs)
%     reRef_data_goodepochs=[reRef_data_goodepochs;reRef_data(:,:,goodepochs(i))];
% end

% Fast Method
split_reRef_data = num2cell(reRef_data(:,:,goodepochs), [1 2]); %split A keeping dimension 1 and 2 intact
reRef_data_goodepochs=vertcat(split_reRef_data{:});

% figure;plot(filtered_data(1:dataEnd,:));
figure;plot(reRef_data_goodepochs(:,goodchans));
title(['Number of good trials: ' num2str(length(goodepochs))]);
