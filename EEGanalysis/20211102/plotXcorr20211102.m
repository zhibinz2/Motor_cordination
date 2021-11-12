cd /home/zhibin/Documents/GitHub/Motor_cordination/EEGanalysis/20211102

run organize_Behavioral_step1.m

cd /home/zhibin/Documents/GitHub/Motor_cordination/EEGanalysis/20211102
%  get detrend_data and filtered_broadband
run organize_EEG_filter_step2.m

% run the above script five times to get different filtered band data and
% save as filtered_data.mat in the the Acquisition folder

% load all the filtered data
cd /home/zhibin/Documents/Acquisition/bimanual_Reach_zhibin_20211018
load filtered_data.mat

%% select a filtered data to process
filtered_data=filtered_broadband;
filtered_data=filtered_delta;
filtered_data=filtered_theta;
filtered_data=filtered_alpha;
filtered_data=filtered_beta;
filtered_data=filtered_gamma;

%% get the good channels
cd /home/zhibin/Documents/GitHub/Motor_cordination/EEGanalysis/20211102

run integrate_EEG_into_data_trials_step3.m
run Remove_EEG_artifact_Rereference_step6.m

%% apply matnlap spline laplacian to the good channels
addpath /home/zhibin/Documents/GitHub


cd /home/zhibin/Documents/GitHub/Motor_cordination/EEGanalysis
load Neuroscan_spherical_chans.mat

Neuroscan_spherical_goodchans=Neuroscan_spherical_chans(goodchans,:);
mat=matnlap(10, Neuroscan_spherical_goodchans, Neuroscan_spherical_chans);
% mat=matnlap(10, Neuroscan_spherical_chans, Neuroscan_spherical_chans);

% laplacian
filtered_data=filtered_data(:, goodchans)*mat;

% reorganzied into data_trials again
cd /home/zhibin/Documents/GitHub/Motor_cordination/EEGanalysis/20211102
run integrate_EEG_into_data_trials_step3.m

% re-reference with the good chans
[reRef_data] = reRef(data_trials(:,1:NumEEGChannels,:),goodchans);


%% plot xcorr for all conditions
UniCondi=unique(CondiData);

figure;tic;
for u=1:length(UniCondi)% u=2
    
    indtemp=find(CondiData==UniCondi(u));
    
    % pick a condition from randomized set allPerm
    conditionSelected = UniCondi(u);
    
%     subplot(1,length(UniCondi),u);
    subplot(1,length(UniCondi),length(UniCondi)+1-u);
    % subplot(7,13,abs(7-conditionSelected)*13+(14-conditionSelected));
    
    rAll=[];rAllMaxLag=[];
    
    yline(0);hold on;

    for indt=1:length(indtemp)


%         FP1EEG = reRef_data(500:1500,1,indtemp(indt)); % Left prefrontal
%         FP2EEG = reRef_data(500:1500,3,indtemp(indt)); % Right prefrontal
%         LeftCortex=FP1EEG;
%         RightCortex=FP2EEG;
        
%         AF7EEG = reRef_data(500:1500,67,indtemp(indt)); % Left prefrontal
%         AF8EEG = reRef_data(500:1500,71,indtemp(indt)); % Right prefrontal
%         LeftCortex=AF7EEG;
%         RightCortex=AF8EEG;
        
%         C3EEG=reRef_data(500:1500,28,indtemp(indt)); % Left motor cortex
%         C4EEG=reRef_data(500:1500,32,indtemp(indt)); % Right motor
%         LeftCortex=C3EEG;
%         RightCortex=C4EEG;

        P5EEG=reRef_data(500:1500,47,indtemp(indt)); % Left motor cortex
        P6EEG=reRef_data(500:1500,53,indtemp(indt)); % Right motor
        LeftCortex=P5EEG;
        RightCortex=P6EEG;
        
        [r,lags]=xcorr(LeftCortex, RightCortex); % if r peak is positive, then the right hemisphere is leading
        rMaxLag=lags(find(r==max(r)));
        xline(rMaxLag);
        rAll=[rAll r];
        rAllMaxLag=[rAllMaxLag rMaxLag];

        
    end
    
    rAllmean=mean(rAll,2);
    %rPositiveMax=max(rAll,[],1);rNegativeMin=min(rAll,[],1);
    
    plot(lags,rAllmean);xlabel('time [ms]');ylabel('xcorr'); xlim([-100 100]); %ylim([-12000 12000]);% ylim([-1*1e4 1*1e4]);
    xline(0,'r');
    % title(['condition: ' conditionNames(u)]);
    title(['condition: ' num2str(conditions(u))]);
    hold off;
    
end

% suptitle(['mean xcorr']);
toc;



%% xcorr for one condition

u=6; % select one condition
    
indtemp=find(CondiData==UniCondi(u));

% find the high peformance trial index within one condition
HighInd=indtemp(find(TrialScores(indtemp)>0.6));
LowInd=indtemp(find(TrialScores(indtemp)<0.6));

% pick a condition from randomized set allPerm
conditionSelected = UniCondi(u)

rAll=[];rAllMaxLag=[];

yline(0);hold on;

for indt=1:length(HighInd)

%     C3EEG=reRef_data(500:1500,28,HighInd(indt)); % Left motor cortex
%     C4EEG=reRef_data(500:1500,32,HighInd(indt)); % Right motor
%     LeftCortex=C3EEG;
%     RightCortex=C4EEG;

%         FP1EEG = reRef_data(500:1500,1,indtemp(indt)); % Left prefrontal
%         FP2EEG = reRef_data(500:1500,3,indtemp(indt)); % Right prefrontal
%         LeftCortex=FP1EEG;
%         RightCortex=FP2EEG;

        AF7EEG = reRef_data(500:1500,67,indtemp(indt)); % Left prefrontal
        AF8EEG = reRef_data(500:1500,71,indtemp(indt)); % Right prefrontal
        LeftCortex=AF7EEG;
        RightCortex=AF8EEG;

    [r,lags]=xcorr(LeftCortex, RightCortex); % if r peak is positive, then the right hemisphere is leading
    rMaxLag=lags(find(r==max(r)));
    xline(rMaxLag);
    rAll=[rAll r];
    rAllMaxLag=[rAllMaxLag rMaxLag];


end

rAllmean=mean(rAll,2);
%rPositiveMax=max(rAll,[],1);rNegativeMin=min(rAll,[],1);

plot(lags,rAllmean);xlabel('time [ms]');ylabel('xcorr'); xlim([-100 100]); %ylim([-12000 12000]);% ylim([-1*1e4 1*1e4]);
xline(0,'r');
% title(['condition: ' conditionNames(u)]);
title(['condition: ' num2str(conditions(u))]);
hold off;
    
