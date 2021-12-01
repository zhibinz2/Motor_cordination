%% Apply Laplacian
cd /home/zhibin/Documents/GitHub/Motor_cordination/EEGanalysis
load Neuroscan_spherical_chans.mat

Neuroscan_spherical_goodchans=Neuroscan_spherical_chans(goodchans,:);
mat=matnlap(10, Neuroscan_spherical_goodchans, Neuroscan_spherical_chans);

% laplacian (after the filter, the last 10 seconds of data is not usable)
filtered_data=filtered_data(:, goodchans)*mat;


%% pick the behaviorally good trials
% TrialScores
