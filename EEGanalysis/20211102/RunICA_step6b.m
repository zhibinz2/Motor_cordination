% Run ICA on reRef_data
% this should be step 6b

%% Recombine all trials into one big long time series

% % Slow Method
% EEGdataShort=reRef_data(:,:,1);
% for i=2:size(reRef_data,3)
%     %EEGdataShortTemp=reRef_data(:,:,i);
%     EEGdataShort=cat(1,EEGdataShort,reRef_data(:,:,i));
% end

% Fast Method
split_reRef_data = num2cell(reRef_data, [1 2]); %split A keeping dimension 1 and 2 intact
EEGdataShort=vertcat(split_reRef_data{:});

%% Recalculate the event array (for display and for EEGLAB)
% So the 4000 timepoints are structured as 
% "500ms padding + 500ms green dot baseline + 500ms plan phase + 1000ms movement + 1000ms show bonus + 500ms padding"

NumTrials=size(reRef_data,3); % NumTrials=350;

% Index for start of trials
% IndEvents=round(linspace(1000,4000*NumTrials,NumTrials));
IndEvents=[1000:4000:4000*NumTrials]
IndEnds=IndEvents+1500;

% Index for EEGLAB
IndEventsEEGLAB=zeros(1,size(reRef_data,1));
IndEventsEEGLAB(IndEvents)=1;
IndEventsEEGLAB(IndEnds)=2;

