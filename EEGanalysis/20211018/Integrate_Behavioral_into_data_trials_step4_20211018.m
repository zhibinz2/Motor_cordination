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
% CondiData;
% figure;plot(CondiData);
% yticks([1 2 3 4]);
% yticklabels({'A1&A2','A1&A1','A2&A1','A2&A2'});title('task sequence'); xlabel('trial');
% unique(CondiData)
% allPerm;
% conditionfuctions;
% conditionNames={'A1&A2','A1&A1','A2&A1','A2&A2'};
% conditionNames=conditions;
% conditionNames{1};
