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

NumTrials=size(reRef_data,3); % NumTrials=350; or numtotal or NumTrialsRecorded

% Index for start of trials
% IndEvents=round(linspace(1000,4000*NumTrials,NumTrials));
IndEvents=[1000:4000:4000*NumTrials];
IndEnds=IndEvents+1500;

% Index for EEGLAB
IndEventsEEGLAB=zeros(1,size(reRef_data,1));
IndEventsEEGLAB(IndEvents)=1;
IndEventsEEGLAB(IndEnds)=2;
plot(IndEventsEEGLAB,'ro');

%% Run ICA
plot(1:size(EEGdataShort,1), EEGdataShort);hold on;
for i=1:length(IndEvents)
    xline(IndEvents(i),'g');
    xline(IndEnds(i),'r');
end
hold off;

tic
[icasig, A, W] = fastica(EEGdataShort');
toc

%% Plot ICA component
for i=1:size(A,2)
    SqA(i)=sumsqr(A(:,i));
end
figure;
plot(1:size(A,2),SqA,'ro');ylabel('sum of square of column in A');xlabel('ICs');
[B,I]=sort(SqA,'descend');


% Component 56 68 69 have the highest weights in the mixing matrics A
% topoplot to examine them
cd /home/zhibin/Documents/GitHub/Motor_cordination/EEGanalysis
load('Neuroscan_spherical__topoplot_chans.mat')
figure;
subplot(3,1,1);topoplot(A(:,69),test,'nosedir','+Y');title('component 69');colorbar;
subplot(3,1,2);topoplot(A(:,56),test,'nosedir','+Y');title('component 56');colorbar;
subplot(3,1,3);topoplot(A(:,68),test,'nosedir','+Y');title('component 68'); colorbar;

figure;
subplot(3,1,1);topoplot(A(:,91),test,'nosedir','+Y');title('component 91');colorbar;
subplot(3,1,2);topoplot(A(:,76),test,'nosedir','+Y');title('component 76');colorbar;
subplot(3,1,3);topoplot(A(:,19),test,'nosedir','+Y');title('component 19'); colorbar;

figure;
subplot(3,1,1);topoplot(A(:,11),test,'nosedir','+Y');title('component 11');colorbar;
subplot(3,1,2);topoplot(A(:,61),test,'nosedir','+Y');title('component 61');colorbar;
subplot(3,1,3);topoplot(A(:,111),test,'nosedir','+Y');title('component 111'); colorbar;

figure;
subplot(3,1,1);topoplot(A(:,3),test,'nosedir','+Y');title('component 3');colorbar;
subplot(3,1,2);topoplot(A(:,53),test,'nosedir','+Y');title('component 53');colorbar;
subplot(3,1,3);topoplot(A(:,101),test,'nosedir','+Y');title('component 101'); colorbar;

%% Display signal before ICA
% addpath /home/zhibin/Documents/GitHub/
addpath(genpath('/home/zhibin/Documents/GitHub/'));

figure;
PlotStart=95;PlotEnd=100;Duration=(IndEvents(PlotStart)-1000):(IndEnds(PlotEnd)+1000);
TimesdataShort=1:size(EEGdataShort,1);
plotx(TimesdataShort(Duration),EEGdataShort(Duration,:));
hold on;
for i=PlotStart:PlotEnd
    xline(IndEvents(i),'r',{'Trial Start'});
    xline(IndEnds(i),'g',{'Trial End'});
end
hold off;
title('EEG Singal Before ICA');

%%

