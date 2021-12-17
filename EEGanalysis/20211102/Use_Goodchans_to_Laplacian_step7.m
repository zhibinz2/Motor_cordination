%% Apply Laplacian
cd /home/zhibin/Documents/GitHub/Motor_cordination/EEGanalysis
load Neuroscan_spherical_chans.mat

Neuroscan_spherical_goodchans=Neuroscan_spherical_chans(goodchans,:);
mat=matnlap(10, Neuroscan_spherical_goodchans, Neuroscan_spherical_chans);

% laplacian (after the filter, the last 10 seconds of data is not usable)
% Interplate the bad channels as well
mixedsig=mixedsig';
laplacian_data=mixedsig(:, goodchans)*mat;

%% just to examine whether edge artifact exists. Answer: No
plot(TimesdataShort,laplacian_data);

%% Plot laplacian_data
plotx(TimesdataShort(Duration),laplacian_data(Duration,:));
hold on;
for i=PlotStart:PlotEnd
    xline(IndEvents(i)+500,'r',{'Target'});
    xline(IndEnds(i)-1000,'g',{'Go'});
end
hold off;
title('Laplacian with ICs removed');

%% Plot before and after laplacian
figure;
subplot(2,1,1);
% PlotStart=95;PlotEnd=100;
% Duration=(IndEvents(PlotStart)-1000):(IndEnds(PlotEnd)+1000);
% TimesdataShort=1:size(EEGdataShort,1);
plotx(TimesdataShort(1,Duration),mixedsig(Duration,:));
hold on;
for i=PlotStart:PlotEnd
    xline(IndStart(i)+500,'r',{'Target'});
    xline(IndEnds(i)-1000,'g',{'Go'});
end
hold off;
title('Mixed Signal with ICs removed');

subplot(2,1,2);
plotx(TimesdataShort(Duration),laplacian_data(Duration,:));
hold on;
for i=PlotStart:PlotEnd
    xline(IndEvents(i)+500,'r',{'Target'});
    xline(IndEnds(i)-1000,'g',{'Go'});
end
hold off;
title('Laplacian with ICs removed and bad channels interpolated');

%% pick the behaviorally good trials
% TrialScores
