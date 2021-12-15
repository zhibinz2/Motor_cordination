%% Apply Laplacian
cd /home/zhibin/Documents/GitHub/Motor_cordination/EEGanalysis
load Neuroscan_spherical_chans.mat

Neuroscan_spherical_goodchans=Neuroscan_spherical_chans(goodchans,:);
mat=matnlap(10, Neuroscan_spherical_goodchans, Neuroscan_spherical_chans);

% laplacian (after the filter, the last 10 seconds of data is not usable)
% Interplate the bad channels as well
mixedsig=mixedsig';
laplacian_data=mixedsig(:, goodchans)*mat;

%% Plot laplacian_data
plotx(TimesdataShort(Duration),laplacian_data(Duration,:));
hold on;
for i=PlotStart:PlotEnd
    xline(IndEvents(i),'r',{'Trial Start'});
    xline(IndEnds(i),'g',{'Trial End'});
end
hold off;
title('Laplacian with ICs removed');


%% pick the behaviorally good trials
% TrialScores
