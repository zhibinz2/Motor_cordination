function [artifact] = artifact_extraction(data_trials_EEG)

% data_trials_EEG=stimTrials_synchronization_b;
% 
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
% figure;imagesc(log10(trialdatamax));colorbar;ylabel('chan(log)');xlabel('trial(log)');

% first round of artifact labelling **************#########################
threshold=0.9*max(trialdatamax,[],'all'); % set a amplitude threshold as criteria to remove bad channels and epochs later

% Create a label matrix that labels bad chan and epoch as 1 
binarythreshhold = zeros(size(trialdatamax));
binarythreshhold(trialdatamax > threshold) = 1;
% figure; imagesc(binarythreshhold); colorbar; ylabel('chan');xlabel('trial');
% title(['sum of Ones: ' num2str(sum(binarythreshhold,'all'))]);

chansum = sum(binarythreshhold,2); % how many bad epochs in each chan
epochsum = sum(binarythreshhold,1); % how many bad chans in each epoch

% bar(chansum); % examine the chans 
badchancrit = prctile(chansum,90); % criteria to determine bad chans
badchanedit = find(chansum > badchancrit); % get the initial bad chans
goodchanedit = setdiff([1:nChans],badchanedit); % get the initial good chans

% bar(epochsum); % examine the epochs
badepochcrit = prctile(epochsum,85); % criteria to determine bad epochs
badepochedit = find(epochsum > badepochcrit); % get the initial bad epochs
goodepochedit = setdiff([1:nEpochs],badepochedit); % get the initial good epochs

% second round of labelling **************#################################
threshold2=100;
% update the label matrix that labels bad chan and epoch as 1 
binarythreshhold = zeros(size(trialdatamax));
binarythreshhold(trialdatamax > threshold2) = 1;
% figure; imagesc(binarythreshhold); colorbar; ylabel('chan');xlabel('trial');
% figure; imagesc(binarythreshhold(goodchanedit,goodepochedit));colorbar; ylabel('chan');xlabel('trial');

epochsum2 = sum(binarythreshhold(goodchanedit,goodepochedit),1); % how many bad chans in each good epoch of the upated label matrix % whos % to look at the workspace
% figure;bar(goodepochedit,epochsum2); xlabel('good trials'); ylabel('number of bad channels');
newbadepochcrit = prctile(epochsum2,80); % criteria to determine bad epochs in 2nd round
newbadepochedit = goodepochedit(find(epochsum2 > newbadepochcrit)); % get the bad epochs in 2nd round
badepochedit = [badepochedit newbadepochedit]; % add the bad epochs in 1s round
goodepochedit_final = setdiff([1:nEpochs],badepochedit);  % get the final good epochs

chansum2 = sum(binarythreshhold(goodchanedit,goodepochedit),2); % how many bad epochs in each good chan of the upated label matrix
% figure;bar(goodchanedit,chansum2);xlabel('good channels');ylabel('number of bad trials');
newbadchancrit = prctile(chansum2, 90); % criteria to determine bad chans in 2nd round
newbadchanedit = goodchanedit(find(chansum2 > newbadchancrit)); % get the bad chans in 2nd round
badchanedit = [badchanedit; newbadchanedit']; % add the bad chans in 1s round
goodchanedit_final = setdiff([1:nChans],badchanedit);  % get the final good chans

% examine the updated label matrix after the 2nd round
% figure;imagesc(trialdatamax(goodchanedit_final,goodepochedit_final));colorbar;
% figure;imagesc(log10(trialdatamax(goodchanedit_final,goodepochedit_final)));colorbar;
% compare with previous one round artifact labelling method
% figure;imagesc(trialdatamax(goodchanedit,goodepochedit));colorbar;
% figure;imagesc(log10(trialdatamax(goodchanedit,goodepochedit)));colorbar;
% figure;imagesc(log10(trialdatamax));colorbar;

badchans = badchanedit;
badepochs = badepochedit;
goodchans = goodchanedit_final;
goodepochs = goodepochedit_final;

artifact.badchans = badchans';
artifact.badepochs = badepochs;
artifact.goodchans = goodchans;
artifact.goodepochs = goodepochs;
artifact.editmatrix = trialdatamax;
artifact.threshhold = threshold;
artifact.threshhold2 = threshold2;
artifact.zerochans = zerochans;
end