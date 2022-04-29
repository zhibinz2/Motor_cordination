function [baselinecorrected_trials] = trials_baselinecorrect(trials_data,baselinesamples)
% baseline=1:sr;
% stimTrials_synchronization-mean(stimTrials_synchronization(baseline,:,:),1)*ones(size(stimTrials_synchronization,1),1);
baselinecorrected_trials=trials_data-mean(trials_data(baselinesamples,:,:),1).*ones(size(trials_data,1),1);
end