clear
cd /home/zhibin/Documents/Artinis_NIRS/zhibin/2022032801
cfg             = [];
cfg.dataset = '2022032801.oxy3';
data_raw        = ft_preprocessing(cfg);

% (demean and view the data, which is noisy) (work!)
cfg                = [];
cfg.preproc.demean = 'yes'; % (subtract the mean value only for plotting, data remains the same)
cfg.viewmode       = 'vertical';
cfg.continuous     = 'no';
cfg.ylim           = [ -0.01   0.01 ];
cfg.channel        = 'Rx*'; % only show channels starting with Rx
ft_databrowser(cfg, data_raw);

%% paddings
trial=cell2mat(data_raw.trial);
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(trial(1:44,:)');
[x, ~] = ginput(2); % read two mouse clicks on the plot % x were index, y were real values
% get the proximate index
string(x)
ind1=round(x(1))
ind2=round(x(2))
% extend
trial_paddings=[repmat(trial(:,ind1),1,ind1-1+size(trial,2)) trial(:,ind1:ind2) repmat(trial(:,ind2),1,size(trial,2)-ind2+size(trial,2))];
plot(trial_paddings(1:44,:)')
% recplace variables in data_raw
fsample=data_raw.fsample;
time=[0:1:size(trial_paddings,2)-1]*1/fsample;
data_raw.time={time};
data_raw.trial={trial_paddings};
data_raw.sampleinfo=[1 size(trial_paddings,2)];
data_raw.hdr.nSamples=size(trial_paddings,2);
data_raw.cfg.trl=[1 size(trial_paddings,2)];


% (demean and view the data, which is noisy) (work!)
cfg                = [];
cfg.preproc.demean = 'yes'; % (subtract the mean value only for plotting, data remains the same)
cfg.viewmode       = 'vertical';
cfg.continuous     = 'no';
cfg.ylim           = [ -0.01   0.01 ];
cfg.channel        = 'Rx*'; % only show channels starting with Rx
ft_databrowser(cfg, data_raw);

%% (transform optical densities to oxy- and deoxy-hemoglobin concentration changes) (work!)
cfg                 = [];
cfg.dpf             = 5.8;
cfg.target          = {'O2Hb', 'HHb'};
cfg.channel         = 'nirs'; % e.g., one channel incl. wildcards, you can also use ?all? to select all NIRS channels
data_conc           = ft_nirs_transform_ODs(cfg, data_raw);
% view data
cfg = [];
cfg.channel = {'Rx*'}; 
cfg.preproc.demean = 'yes'; % (subtract the mean value only for plotting, data remains the same)
cfg.viewmode       = 'vertical';
cfg.continuous     = 'no';
cfg.ylim           = [ -0.3   0.3 ];
ft_databrowser(cfg, data_conc);

%% (band-pass filter: [0.05 0.06] Hz)
% (low-pass filtering) (work!)
cfg                   = [];
cfg.lpfilter          = 'yes';
cfg.lpfreq            = 0.2;
data_lpf              = ft_preprocessing(cfg, data_conc);
% view data
cfg = [];
cfg.channel = {'Rx*'}; 
cfg.preproc.demean = 'yes'; % (subtract the mean value only for plotting, data remains the same)
cfg.viewmode       = 'vertical';
cfg.continuous     = 'no';
cfg.ylim           = [ -0.3   0.3 ];
ft_databrowser(cfg, data_lpf);

% % (high-pass filtering)(large edge artifact at the start)
% cfg                 = [];
% cfg.hpfilter        = 'yes';
% cfg.hpfreq          = 0.05;
% data_flt            = ft_preprocessing(cfg,data_lpf); 
% % view data
% cfg = [];
% cfg.channel = {'Rx*'}; 
% cfg.preproc.demean = 'yes'; % (subtract the mean value only for plotting, data remains the same)
% cfg.viewmode       = 'vertical';
% cfg.continuous     = 'no';
% cfg.ylim           = [ -0.3   0.3 ];
% ft_databrowser(cfg, data_flt);

%% Pick the signal from fieldtrip result
filtered_data=cell2mat(data_lpf.trial);
filtered_data=filtered_data(1:44,:)';
Fs=data_flt.fsample;

%% load and plot singal from oxy3.mat to examine ylim and events
% oxysoft2matlab % select oxy3 and oxyproj files
% load('2022042104oxy3.mat')

%% detect events from photocells
% % Use photocell in oxy3
% time=nirs_data.time;
% ADlabel=nirs_data.ADlabel;
% ADvalues=nirs_data.ADvalues;
% figure;
% for i=1:16
%     subplot(4,4,i)
%     plot(ADvalues(:,i));
%     title(num2str(i));
% end
% Photocell=ADvalues(:,15);plot(Photocell);
% Photocell=ADvalues(:,16);plot(Photocell);% this channel is better

% Use photocell from fieldtrip
trial=cell2mat(data_raw.trial);
figure;
for i=1:16
    subplot(4,4,i)
    plot(trial(44+i,:)');
    title(num2str(i));
end
Photocell=trial(44+15,:)';plot(Photocell);

% find peaks
% view the time course of photocell signals
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(Photocell);xlabel('time');ylabel('photocell signal');
% click and select the start and end point for peak extraction
[x, y] = ginput(2); % read two mouse clicks on the plot % x were index, y were real values
Startpoint=round(x(1));Endpoint=round(x(2)); % Startpoint and Endpoint are sample number or index in time
hold on;xline(x(1),'r');xline(x(2),'r');hold off;

% replace the beginning and end with baseline value
Photocell(1:Startpoint)=mean(y);Photocell(Endpoint:end)=mean(y); % plot(Photocell');
plot(Photocell,'b'); 

% Examine peaks detection in analog1
Halfhigh1=1/4*(max(Photocell)-min(Photocell));% +min(Photocell); % value for 'MinPeakProminence'
% Check if need to adjust the Halfhigh cutoff
close;figure;
findpeaks(Photocell,time,'MinPeakProminence',Halfhigh1,'Annotate','extents');
yline(Halfhigh1,'m','MinPeakProminence');
ylim([0 max(Photocell)+1000]);

% locate the trial sessions % pks=value of the peak % locs=time of the peak
[pks,locs] = findpeaks(Photocell,time,'MinPeakProminence',Halfhigh1,'Annotate','extents');

% find the event indices in time
EventIndices=zeros(length(locs),1);
for i=1:length(locs)
    EventIndices(i)=find(time==locs(i));
end
%% examine EventIndices on top of signal
figure('units','normalized','outerposition',[0 0 1 0.3]);
plotx(time,filtered_data(:,1:2:43),'r');
hold on;
plotx(time,filtered_data(:,2:2:44),'b');
% title('oxyval & dxyval 0.2~0.8 Hz'); 
hold on;
for i=1:length(locs)/2
    xline(time(EventIndices(2*i-1)),'m',{allPerm(i)});
end
hold off;
ylim([-1 1])
title('2022032801 only low pass filt')

%% organized into trials based on photocells
data_trials=zeros(round(fsample*40),44,length(locs)/2);
for i=1:length(locs)/2
    ind1=EventIndices(2*i-1);
    data_trials(:,:,i)=filtered_data((ind1-round(fsample*10)):ind1+(round(fsample*30)-1),:);
end
% deduct the avagerage from trial
data_trials=data_trials-mean(data_trials,1).*ones(size(data_trials,1),1);
% deduct the baseline average 
data_trials=data_trials-mean(data_trials(1:round(fsample*10),:,:),1).*ones(size(data_trials,1),1);
% divided by the baseline average 
data_trials=data_trials./(mean(data_trials(1:round(fsample*10),:,:),1).*ones(size(data_trials,1),1));

% channels
oxychanL=1:2:21;
dxychanL=2:2:22;
oxychanR=23:2:43;
dxychanR=24:2:44;

% Only the occipital lobe
oxychanL=[1 5 7 9 11 17];
dxychanL=[2 6 8 10 12 18];
oxychanR=[23 27 29 31 33 35];
dxychanR=[24 28 30 32 34 36];

% Only the most likely visual cortex
oxychanL=[7 9 11];
dxychanL=[8 10 12];
oxychanR=[29 31 33];
dxychanR=[30 32 34];

% plot all trials as subplots
timeTrials=[0:1:round(fsample*40)-1]*1/fsample-10;
figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:72
    subplot(8,9,i)
    plot(timeTrials,data_trials(:,1:2:43,i),'r');hold on;
    plot(timeTrials,data_trials(:,2:2:44,i),'b');hold off;
    xlabel('time (s)');
%     ylim([-2.5 2.5]);
    title('i');
end
suptitle('recording 2022032801')
% plot the average
data_trials_m=mean(data_trials,3);
figure;
subplot
plotx(timeTrials,data_trials_m(:,oxychanL),'r');hold on;
plotx(timeTrials,data_trials_m(:,dxychanL),'b');hold off;
xlabel('time (s)');
% ylim([-2.5 2.5]);
title('average across 72 trials');
% plot 2 conditions
trialsL=find(allPerm==1);
trialsR=find(allPerm==2);
data_trials_mL=mean(data_trials(:,:,trialsL),3);
data_trials_mR=mean(data_trials(:,:,trialsR),3);
figure;
subplot(2,2,1);
plotx(timeTrials,data_trials_mL(:,oxychanL),'r');hold on;
plotx(timeTrials,data_trials_mL(:,dxychanL),'b');hold off;
title('Left stim, Left Brain');
hold on;xline(0,'m',{'stim onset'});xline(10,'m',{'stim end'});hold off;
ylabel('concentration of oxygenated (red) and deoxygenated (blue) Hb');
xlabel('time (s)');ylim([-0.3 0.5]);
subplot(2,2,2);
plotx(timeTrials,data_trials_mL(:,oxychanR),'r');hold on;
plotx(timeTrials,data_trials_mL(:,dxychanR),'b');hold off;
title('Left stim, Right Brain');
hold on;xline(0,'m',{'stim onset'});xline(10,'m',{'stim end'});hold off;
ylabel('concentration of oxygenated (red) and deoxygenated (blue) Hb');
xlabel('time (s)');ylim([-0.3 0.5]);
subplot(2,2,3);
plotx(timeTrials,data_trials_mR(:,oxychanL),'r');hold on;
plotx(timeTrials,data_trials_mR(:,dxychanL),'b');hold off;
title('Right stim, Left Brain');
hold on;xline(0,'m',{'stim onset'});xline(10,'m',{'stim end'});hold off;
ylabel('concentration of oxygenated (red) and deoxygenated (blue) Hb');
xlabel('time (s)');ylim([-0.3 0.5]);
subplot(2,2,4);
plotx(timeTrials,data_trials_mR(:,oxychanR),'r');hold on;
plotx(timeTrials,data_trials_mR(:,dxychanR),'b');hold off;
title('Right stim, Right Brain');
hold on;xline(0,'m',{'stim onset'});xline(10,'m',{'stim end'});hold off;
ylabel('concentration of oxygenated (red) and deoxygenated (blue) Hb');
xlabel('time (s)');ylim([-0.3 0.5]);
suptitle('Recording2022032801')


%% 2 subplots of all trials
figure;
Ylim=3;
subplot(2,1,1); 
plot(time,filtered_data(:,oxychanL),'r',time,filtered_data(:,dxychanL),'b');
title('Left Hemisphere condition ');ylim([-1*Ylim Ylim]);
subplot(2,1,2); 
plot(time,filtered_data(:,oxychanR),'r',time,filtered_data(:,dxychanR),'b');
title('Right Hemisphere condition ');ylim([-1*Ylim Ylim]);
for i=1:numTrials
    subplot(2,1,1); 
    hold on;xline(time(EventIndices(2*i-1)),'m',{allPerm(i)});xline(time(EventIndices(2*i)),'k');hold off;
    subplot(2,1,2); 
    hold on;xline(time(EventIndices(2*i-1)),'m',{allPerm(i)});xline(time(EventIndices(2*i)),'k');hold off;
end