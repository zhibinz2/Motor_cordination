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

%% Recalculate the events
