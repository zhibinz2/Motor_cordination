%% select the EEG channels

addpath(genpath('/home/zhibin/Documents/GitHub'));
% data=dataEEG(1:128,2000000:2002000);
EEGdata=dataEEG(1:128,:);

%% Appen EventChannel for EEGLAB
% EEGdata(129,:)=EventChannel;

%% detrend the data

% tranpose the data so that time is the first dimension, channel is 2nd dimension, and trial is the 3rd dimension
transpose_data=permute(EEGdata,[2,1]);
% plotx(mean(transpose_data,3));
% plotx(transpose_data);hold on;xline(dataEnd,'m--');

% detrend the data
detrend_data=ndetrend(transpose_data,1); 
% plotx(mean(detrend_data,2)); xlim([5000 10000]);
% plotx(detrend_data); xlim([5000 10000]);
% plotx(detrend_data(5000:10000,:));

%% Broadband filter
% apply high pass filter to filter out frequencies lower than 1;
% load eeglab toolbox to avoid filtered_data becoming NaN
tic
eeglab;
close all;
Hd = makefilter(Fs,1.5,1,3,20,0); % Astop 6 or 20;
filtered_data=filtfilthd(Hd,double(detrend_data));

% remove channel 91 44
% filtered_data(:,91)=[];filtered_data(:,44)=[];
%  plotx(filtered_data); xlim([5000 10000]);

% apply low pass filter to filter out frequencies higher than 50;
Hd = makefilter(Fs,50,55,3,20,0);
filtered_data=filtfilthd(Hd,filtered_data);
toc

% plotx(mean(filtered_data,2));  
% plot(filtered_data(1:dataEnd,:));
% hold on;
% for i=1:length(ind_trial_start)
%     xline(locs(ind_trial_start(i)),'k');
% end
% xline(dataEnd,'m--');xlim([0 dataEnd]);
% plotx(filtered_data(5000:10000,:));


%%
% filtered_broadband=filtered_data;
% filtered_data=filtered_broadband;

%% fiter5bands 
% tic;
% run filter5bands.m
% toc;

% filtered_data=filtered_delta;
% filtered_data=filtered_theta;
% filtered_data=filtered_alpha;
% filtered_data=filtered_beta;
% filtered_data=filtered_gamma;

% plot(filtered_beta(5000:10000,8));
% figure;
% plot(filtered_gamma(5000:10000,8));

% tic;
% % cd /home/zhibin/Documents/Acquisition/Bimanual_reach_zhibin_20211102
% filename='filtered_data.mat';
% % save(filename); % save all variables in the workspace
% save(filename, 'filtered_broadband', 'filtered_delta', 'filtered_theta', 'filtered_alpha', 'filtered_beta', 'filtered_gamma', '-v7.3');
% toc;