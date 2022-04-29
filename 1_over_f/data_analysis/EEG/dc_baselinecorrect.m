function [trial_data_b] = dc_baselinecorrect(trial_data,sr)
% baselinesamps = 1:2*sr; % correct DC offset using the mean of the whole trials
% stimTrials_synchronization_b=trials_baselinecorrect(stimTrials_synchronization,baselinesamps);
% baselinesamps = 1:sr; % baseline correct
% stimTrials_synchronization_b=trials_baselinecorrect(stimTrials_synchronization_b,baselinesamps);
baselinesamps = 1:2*sr; % correct DC offset using the mean of the whole trials
trial_data_dc=trials_baselinecorrect(trial_data,baselinesamps);
baselinesamps = 1:sr; % baseline correct
trial_data_b=trials_baselinecorrect(trial_data_dc,baselinesamps);
end