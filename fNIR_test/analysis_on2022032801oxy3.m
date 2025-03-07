clear
% cd C:\Users\NIRS\Documents\zhibin\2022032801
cd /home/zhibin/Documents/Artinis_NIRS/zhibin/20220410
% addpath /home/zhibin/Documents/GitHub/fieldtrip/fieldtrip-20220304
%% single channel demo: load data with fieldtrip
cfg = [];
% cfg.dataset = '2022032801.oxy3'; % only oxy3 save the events 
cfg.dataset = '2022041001.oxy3'; % only oxy3 save the events 
[data]=ft_preprocessing(cfg);

cfg = [];
cfg.ylim = 'maxabs';
% cfg.ylim = 'maxmin';
ft_databrowser(cfg, data);

cfg = [];
cfg.ylim = 'maxmin';
cfg.channel = {'Rx*'}; 
ft_databrowser(cfg, data);

cfg                = [];
cfg.preproc.demean = 'yes'; % (subtract the mean value only for plotting, data remains the same)
% cfg.viewmode       = 'vertical';
cfg.continuous     = 'no';
cfg.ylim           = [ -0.003   0.003 ];
cfg.channel        = {'Rx*'}; % only show channels starting with Rx
ft_databrowser(cfg, data);
%% % extra all variables from the data (too many!)
hdr=data.hdr;
    Fs=hdr.Fs;
    nChans=hdr.nChans;
    label=hdr.label;
    nSamples=hdr.nSamples;
    nTrials=hdr.nTrials;
    chantype=hdr.chantype;
    chanunit=hdr.chanunit;
    opto=hdr.opto;
        DPF=opto.DPF;
        chanpos=opto.chanpos;
        chanlabel=opto.label;
        optolabel=opto.optolabel;
        optopos=opto.optopos;
        tra=opto.tra;
        type=opto.type;
        unit=opto.unit;
        wavelength=opto.wavelength;
label=data.label;
time=cell2mat(data.time);
trial=cell2mat(data.trial);
fsample=data.fsample;
sampleinfo=data.sampleinfo;
opto=data.opto;
    DPF=opto.DPF;
    chanpos=opto.chanpos;
    chanlabel=opto.label;
    optolabel=opto.optolabel;
    optopos=opto.optopos;
    tra=opto.tra;
    type=opto.type;
    unit=opto.unit;
    wavelength=opto.wavelength;
cfg=data.cfg;
    dataset=cfg.dataset;
    checkpath=cfg.checkpath;
    showlogo=cfg.showlogo;
    outputfilepresent=cfg.outputfilepresent;
    tracktimeinfo=cfg.tracktimeinfo;
    trackmeminfo=cfg.trackmeminfo;
    toolbox=cfg.toolbox;
        images=toolbox.images;
        stats=toolbox.stats;
        signal=toolbox.signal;
        cleanup=toolbox.cleanup;
    callinfo=cfg.callinfo;
        usercfg=callinfo.usercfg;
        fieldtrip=callinfo.fieldtrip;
        matlab=callinfo.matlab;
        computer=callinfo.computer;
        hostname=callinfo.hostname;
        user=callinfo.user;
        pwd=callinfo.pwd;
        calltime=callinfo.calltime;
        proctime=callinfo.proctime;
    version=cfg.version;
        name=version.name;
        id=version.id;
    method=cfg.method;
    channel=cfg.channel;
    removemcg=cfg.removemcg;
    removeeog=cfg.removeeog;
    precision=cfg.precision;
    padding=cfg.padding;
    paddir=cfg.paddir;
    montage=cfg.montage;
    updatesens=cfg.updatesens;
    dataformat=cfg.dataformat;
    dftfilter=cfg.dftfilter;
    lpfilter=cfg.lpfilter;
    hpfilter=cfg.hpfilter;
    bpfilter=cfg.bpfilter;
    bsfilter=cfg.bsfilter;
    medianfilter=cfg.medianfilter;
    padtype=cfg.padtype;
    reref=cfg.reref;
    refchannel=cfg.refchannel;
    refmethod=cfg.refmethod;
    groupchans=cfg.groupchans;
    implicitref=cfg.implicitref;
    feedback=cfg.feedback;
    datafile=cfg.datafile;
    headerfile=cfg.headerfile;
    headerformat=cfg.headerformat;
    continuous=cfg.continuous;
    trl=cfg.trl;
    leadfield=cfg.leadfield;
    polyremoval=cfg.polyremoval;
    polyorder=cfg.polyorder;
    detrend=cfg.detrend;
    demean=cfg.demean;
    baselinewindow=cfg.baselinewindow;
    dftfreq=cfg.dftfreq;
    lpfiltord=cfg.lpfiltord;
    hpfiltord=cfg.hpfiltord;
    bpfiltord=cfg.bpfiltord;
    bsfiltord=cfg.bsfiltord;
    lpfilttype=cfg.lpfilttype;
    hpfilttype=cfg.hpfilttype;
    bpfilttype=cfg.bpfilttype;
    bsfilttype=cfg.bsfilttype;
    lpfiltdir=cfg.lpfiltdir;
    hpfiltdir=cfg.hpfiltdir;
    bpfiltdir=cfg.bpfiltdir;
    bsfiltdir=cfg.bsfiltdir;
    lpinstabilityfix=cfg.lpinstabilityfix;
    hpinstabilityfix=cfg.hpinstabilityfix;
    bpinstabilityfix=cfg.bpinstabilityfix;
    bsinstabilityfix=cfg.bsinstabilityfix;
    lpfiltdf=cfg.lpfiltdf;
    hpfiltdf=cfg.hpfiltdf;
    bpfiltdf=cfg.bpfiltdf;
    bsfiltdf=cfg.bsfiltdf;
    lpfiltwintype=cfg.lpfiltwintype;
    hpfiltwintype=cfg.hpfiltwintype;
    bpfiltwintype=cfg.bpfiltwintype;
    bsfiltwintype=cfg.bsfiltwintype;
    lpfiltdev=cfg.lpfiltdev;
    hpfiltdev=cfg.hpfiltdev;
    bpfiltdev=cfg.bpfiltdev;
    bsfiltdev=cfg.bsfiltdev;
    plotfiltresp=cfg.plotfiltresp;
    usefftfilt=cfg.usefftfilt;
    medianfiltord=cfg.medianfiltord;
    hilbert=cfg.hilbert;
    derivative=cfg.derivative;
    rectify=cfg.rectify;
    boxcar=cfg.boxcar;
    absdiff=cfg.absdiff;
    conv=cfg.conv;
    dftinvert=cfg.dftinvert;
    standardize=cfg.standardize;
    denoise=cfg.denoise;
    subspace=cfg.subspace;
    custom=cfg.custom;
    resample=cfg.resample;
    previous=cfg.previous;
%% % Exercise 1 (View Artifact detection) (not good)
cfg = [];
cfg.artfctdef.zvalue.channel = {'Rx1-Tx1 [760nm]', 'Rx1-Tx1 [842nm]'};
% cfg.artfctdef.zvalue.channel = {'Rx2*'};
cfg.artfctdef.zvalue.cutoff = 5;
cfg.artfctdef.zvalue.hpfilter = 'yes';
cfg.artfctdef.zvalue.hpfreq = 0.1;
cfg.artfctdef.zvalue.rectify = 'yes';
cfg.artfctdef.zvalue.artpadding = 2;
cfg.artfctdef.zvalue.interactive = 'yes'; % the interactive display makes more sense after segmentating data in trials
[cfg, artifact] = ft_artifact_zvalue(cfg, data);
% detected 8 artifacts, call ft_refectartifact to remove after filtering
% and segmenting the data into epochs

%% % Exercise 2 (Covert OD to concentration) (work!)
cfg = [];
cfg.dpf = 6.2;
cfg.channel = {'Rx*'};
data_conc = ft_nirs_transform_ODs(cfg, data);
    % examine data_conc
    conc_trial=data_conc.trial;

%% % Exercise 3 (bandpass filter) % not working
% cfg = [];
% cfg.ylim = 'maxmin';
% cfg.channel = {'Rx*'}; 
% ft_databrowser(cfg, data_conc);

cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = [0.01 0.8]; % (bandpass filter)
data_filtered = ft_preprocessing(cfg, data_conc); 

%% % (Define epochs of interest) not work
help ft_definetrial

cfg = [];
cfg.dataset = '2022032801.oxy3';
cfg.trialdef.eventtype = '?';

ft_definetrial(cfg); % no event detected in cfg

% [trl, event] = ft_trialfun_general(cfg) 

cfg.trialdef.eventtype  = 'event';
cfg.trialdef.eventvalue = 'R';
cfg.trialdef.prestim    = 10;
cfg.trialdef.poststim   = 35;
cfg = ft_definetrial(cfg);

cfg.channel = {'Rx*'};
data_epoch = ft_redefinetrial(cfg, data_filtered); % We have now selected one pair of channels and cut the data in 12 trials. Check them out using the databrowser, but let us use some settings to make the plots look neater and also visualize the artifacts that we identified earlier:
%% multi channels demo: load data with fieldtrip
clear;
addpath  /home/zhibin/Documents/GitHub/fieldtrip/fieldtrip-20220304
cfg             = [];
cd /home/zhibin/Documents/Artinis_NIRS/zhibin/
% cfg.dataset = '2022032801.oxy3';
% cfg.dataset = '2022041001.oxy3';
% cfg.dataset = '2022041702.oxy3';
cfg.dataset = '2022042104.oxy3';
data_raw        = ft_preprocessing(cfg);

% % (To retrieve the layout from the data file)
% cfg           = [];
% cfg.opto      = '2022032801.oxy3';
% ft_layoutplot(cfg);

% % all 14 channels
% label={'ADC001','ADC002','ADC003','ADC004','ADC005','ADC006','ADC007','ADC008','ADC009','ADC010','ADC011','ADC012','ADC013','ADC014','ADC015'};
% % (Detecting triggers)
% find(strcmp(data_raw.label,'ADC014'))
% find(strcmp(data_raw.label,'ADC015')) % good for photocell;
% % plot the voltage of ADC001 to ADC015
% for i=1:15 % 45:59=1+44:15+44;
%     subplot(4,4,i);
%     plot(data_raw.time{1}, data_raw.trial{1}(i+44,:)*1, 'r');
%     title(['ADC ' num2str(i)])
% end

%% % Exercise 1 (Detect onset in ADC channel by finding upward going flank as event)
% event = ft_read_event('2022032801.oxy3')
% ADC014 and ADC015 are the channels for events
% adc015 = find(strcmp({event.type}, 'ADC015')); % number not right

% (check sampling rate)
data_raw.fsample

% (demean and view the data, which is noisy) (work!)
cfg                = [];
cfg.preproc.demean = 'yes'; % (subtract the mean value only for plotting, data remains the same)
cfg.viewmode       = 'vertical';
cfg.continuous     = 'no';
cfg.ylim           = [ -0.01   0.01 ];
cfg.channel        = 'Rx*'; % only show channels starting with Rx
ft_databrowser(cfg, data_raw);

% % (high pass filtering, throw away very low-frequency) (work! only when > 0.05 Hz)
% cfg                 = [];
% cfg.hpfilter        = 'yes';
% cfg.hpfreq          = 0.05;
% data_hpf            = ft_preprocessing(cfg,data_raw);

% % (Epoch) (not working)
% event = ft_read_event('2022032801.oxy3');
% adc014 = find(strcmp({event.type}, 'ADC014'));
% adc015 = find(strcmp({event.type}, 'ADC015'));

% % get the sample number in the original data
% % note that we transpose them to get columns
% smp014 = [event(adc014).sample]';
% smp015 = [event(adc015).sample]';

% % (remove bad channels:poor contact with the skin of the scalp) not working
% cfg      = [];
% data_sci = ft_nirs_scalpcouplingindex(cfg, data_raw);
% % you can see that we throw away some channels in data_sci.label, it now
% % has 86 instead of 104 channels

% % (View Artifact detection) (not good)
% cfg = [];
% cfg.artfctdef.zvalue.channel = {'Rx*'};
% cfg.artfctdef.zvalue.cutoff = 5;
% cfg.artfctdef.zvalue.hpfilter = 'yes';
% cfg.artfctdef.zvalue.hpfreq = 0.1;
% cfg.artfctdef.zvalue.rectify = 'yes';
% cfg.artfctdef.zvalue.artpadding = 2;
% cfg.artfctdef.zvalue.interactive = 'yes'; % the interactive display makes more sense after segmentating data in trials
% [cfg, artifact] = ft_artifact_zvalue(cfg, data_raw);
% % detected 8? artifacts, call ft_refectartifact to

% (transform optical densities to oxy- and deoxy-hemoglobin concentration changes) (work!)
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

%% (0.05~0.8 Hz band-pass filter: >0.05 high pass filt on <0.8 data_lpf)
% (low-pass filtering) (work!)
cfg                   = [];
cfg.lpfilter          = 'yes';
cfg.lpfreq            = 0.8;
data_lpf              = ft_preprocessing(cfg, data_conc);
% view data
cfg = [];
cfg.channel = {'Rx*'}; 
cfg.preproc.demean = 'yes'; % (subtract the mean value only for plotting, data remains the same)
cfg.viewmode       = 'vertical';
cfg.continuous     = 'no';
cfg.ylim           = [ -0.3   0.3 ];
ft_databrowser(cfg, data_lpf);

% (high-pass filtering)(large edge artifact at the start)
cfg                 = [];
cfg.hpfilter        = 'yes';
cfg.hpfreq          = 0.05;
data_flt            = ft_preprocessing(cfg,data_lpf); 
% view data
cfg = [];
cfg.channel = {'Rx*'}; 
cfg.preproc.demean = 'yes'; % (subtract the mean value only for plotting, data remains the same)
cfg.viewmode       = 'vertical';
cfg.continuous     = 'no';
cfg.ylim           = [ -0.3   0.3 ];
ft_databrowser(cfg, data_flt);

%% load and plot singal from oxy3.mat to examine ylim and events
oxysoft2matlab % select oxy3 and oxyproj files
load('2022041001oxy3.mat')

% signal 
oxyvals=nirs_data.oxyvals;
dxyvals=nirs_data.dxyvals;
time=nirs_data.time;
Fs=nirs_data.Fs;
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(time,oxyvals,'r',time,dxyvals,'b');title('raw oxyval & dxyval'); 
% detrend
detrend_oxyvals=detrend(oxyvals,1);detrend_dxyvals=detrend(dxyvals,1);
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(time,detrend_oxyvals,'r',time,detrend_dxyvals,'b');title('detrend oxyval & dxyval'); 
ylim([-10 10]);
% high pass (large edge artifact at the start)
Hd = makefilter(Fs,0.05,0.04,6,20,0); % xlim([0 1]);
which filtfilt.m -all
% cd 'C:\Program Files\MATLAB\R2021b\toolbox\signal\signal\'
cd /usr/local/MATLAB/R2022a/toolbox/signal/signal/
which filtfilt.m
filtered_oxyvals1=filtfilthd(Hd,detrend_oxyvals);
filtered_dxyvals1=filtfilthd(Hd,detrend_dxyvals);
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(time,filtered_oxyvals1,'r',time,filtered_dxyvals1,'b');title('oxyval & dxyval > 0.05Hz'); 
ylim([-20 20]);ylim([-3 3]);
% low pass (this will create short edge artifact, no padding for now)
Hd = makefilter(Fs,0.8,0.9,6,20,0);
which filtfilt.m -all
% cd 'C:\Program Files\MATLAB\R2021b\toolbox\signal\signal\'
which filtfilt.m
filtered_oxyvals2=filtfilthd(Hd,filtered_oxyvals1);
filtered_dxyvals2=filtfilthd(Hd,filtered_dxyvals1);
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(time,filtered_oxyvals2,'r',time,filtered_dxyvals2,'b');title('oxyval & dxyval 0.05~0.8 Hz'); 
ylim([-20 20]);ylim([-1.5 1.5]);

%% plot channel position based on oxy3.mat
chanlabel=nirs_data.label;
chanPos=nirs_data.chanPos;
figure;
plot(chanPos(1,:),chanPos(2,:),'r.');
text(chanPos(1,:),chanPos(2,:),chanlabel,'VerticalAlignment','bottom','HorizontalAlignment','right');
hold on;
text(chanPos(1,:),chanPos(2,:),num2cell([1:22]),'VerticalAlignment','top','HorizontalAlignment','left','Color',[1 0 0]);
ylim([-7 4]);xlim([-2 13]);

chanSequence=[];
for i=1:22
    Appendchan=[1 2]+2*(i-1);
    chanSequence=[chanSequence; Appendchan];
end

chanSequenceLabel={};
for i=1:22
    AppendchanLabel=num2str(chanSequence(i,:));
    chanSequenceLabel{i}=AppendchanLabel;
end
    
text(chanPos(1,:),chanPos(2,:),chanSequenceLabel,'VerticalAlignment','top','HorizontalAlignment','right','Color',[0 1 0.5]);
hold off;

% AD channel
ADvalues=nirs_data.ADvalues;
ADlabel=nirs_data.ADlabel;
for i=1:size(ADvalues,2)
    subplot(4,4,i);
    plot(time, ADvalues(:,i));
    title(ADlabel{i});
end 
%% use event onset in nirs_data from oxy3.mat
events=nirs_data.events;
onsets=events.onsets{1,1};
ADlabel=nirs_data.ADlabel;
ADvalues=nirs_data.ADvalues;

figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(time, ADvalues(:,16),'b');
hold on;
for i=1:length(onsets)
    plot(time(onsets(i)),10000*ones(length(onsets),1),'r.');
end
hold off;
legend('photocell signal','event onsets generated by Oxysoft');

% manually picked from the auto-generated events
onsets([1:3 end-2:end])=NaN; % remove the first 3 and last 3, judged from the above plot
onsets([1:3 end-2:end])=[]; % 72 = only 9 blocks x 4 trials x 2 photocells (start and end of stimulus)

%% detect events from photocells in oxy3
ADlabel=nirs_data.ADlabel;
ADvalues=nirs_data.ADvalues;
figure;
for i=1:16
    subplot(4,4,i)
    plot(ADvalues(:,i));
    title(num2str(i));
end
Photocell=ADvalues(:,15);plot(Photocell);
Photocell=ADvalues(:,16);plot(Photocell);% this channel is better

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
plot(time,Photocell,'b'); 

% Examine peaks detection in analog1
Halfhigh1=1/4*(max(Photocell)-min(Photocell));% +min(Photocell); % value for 'MinPeakProminence'
% Check if need to adjust the Halfhigh cutoff
close;figure;
findpeaks(Photocell,time,'MinPeakProminence',Halfhigh1,'Annotate','extents');
yline(Halfhigh1,'m','MinPeakProminence');
ylim([0 max(Photocell)+1000]);

% locate the trial sessions % pks=value of the peak % locs=time of the peak
[pks,locs] = findpeaks(Photocell,time,'MinPeakProminence',Halfhigh1,'Annotate','extents');

% examine locs on top of signal plot
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(time,filtered_oxyvals2,'r',time,filtered_dxyvals2,'b');title('oxyval & dxyval 0.05~0.8 Hz'); 
% ylim([-20 20]);ylim([-1.5 1.5]);
hold on;
for i=1:length(locs)
    xline(locs(i),'m');
end

find(time==locs(2)) % has a return value
% find the event indices in time
EventIndices=zeros(length(locs),1);
for i=1:length(locs)
    EventIndices(i)=find(time==locs(i));
end

% examine EventIndices on top of signal
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(filtered_oxyvals2,'r');
hold on;
plot(filtered_dxyvals2,'b');
title('oxyval & dxyval 0.05~0.8 Hz'); 
hold on;
for i=1:length(EventIndices)
    xline(EventIndices(i),'m');
end

%% Pick the signal from fieldtrip result
filtered_data=cell2mat(data_flt.trial);
filtered_data=filtered_data(1:44,:)';
Fs=data_flt.fsample;
time=cell2mat(data_flt.time);

%% Use the event onset from oxy3.mat to organize data into time x chan x trials
numTrials=length(onsets)/2;
% numTrials=12;
numChans=44;
trial_length=round(Fs*40); % 10s baseline + 10s stimulus + 20 s rest
% trial_length=round(Fs*30); % 10s stimulus + 20 s rest
time(onsets(i))
TrialTime=1/Fs*[1:trial_length]-10;
% TrialTime=1/Fs*[1:trial_length];
data_trials=zeros(trial_length,numChans,numTrials);
for i=1:numTrials % i=numTrials
    data_trials(:,:,i)=filtered_data((onsets(2*i-1)-round(Fs*10)):(onsets(2*i-1)+round(Fs*30)-1),:);
end
% examine one trial
figure;plotx(data_trials(:,:,20));
figure;
for i=1:numTrials
    plot(TrialTime,data_trials(:,1:2:43,i));hold on;
end
hold off;

%% Use photocell detection generated EventIndices to organize data into time x chan x trials
numTrials=12;
numChans=44;
% trial_length=round(Fs*40); % 10s baseline + 10s stimulus + 20 s rest
trial_length=round(Fs*30); % 10s stimulus + 20 s rest
time(EventIndices(1))
% TrialTime=1/Fs*[1:trial_length]-10;
data_trials=zeros(trial_length,numChans,numTrials);
% for i=1:numTrials % i=numTrials
%     data_trials(:,:,i)=filtered_data((EventIndices(2+2*i-1)-round(Fs*10)):(EventIndices(2+2*i-1)+round(Fs*30)-1),:);
% end
for i=1:numTrials % i=numTrials
    data_trials(:,:,i)=filtered_data(EventIndices(2*i-1):(EventIndices(2*i-1)+round(Fs*30)-1),:);
end


%% baseline normalization (skip)
normalized_trials=zeros(trial_length,numChans,numTrials);
for i=1:numTrials
    baselineMean=mean(abs(data_trials(1:round(Fs*10),:,i)),1);
    normalized_trials(:,:,i)=data_trials(:,:,i)./baselineMean;
end
% examine every trial
figure;plot(normalized_trials(:,:,10));
figure;
for i=1:numTrials
    plot(TrialTime,normalized_trials(:,1:2:43,i));hold on;
end
hold off;

% Recombine all trials into one big long time series
% cd /home/zhibin/Documents/GitHub/Motor_cordination/EEGanalysis/20211102
% open RunICA_step6b.m
% Fast Method
split_data = num2cell(normalized_trials, [1 2]); %split A keeping dimension 1 and 2 intact
Combined_norm_data=vertcat(split_data{:});

%% extract the conditions and channels
% load stimulus data
cd C:\Users\NIRS\Documents\zhibin\2022032801
cd /home/zhibin/Documents/Artinis_NIRS/zhibin/20220415
load('2022032801.mat'); allPerm;
numTrials=12;
allPerm=allPerm(1:numTrials);
trialsL=find(allPerm==1);
trialsR=find(allPerm==2);

% separate left and right stimulus condition
data_trialsL=data_trials(:,:,trialsL);
data_trialsR=data_trials(:,:,trialsR);

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

%% 2 subplots trial by trial
% examine trial to trial
figure;
Ylim=3;
for i=1:numTrials
    subplot(2,1,1); 
    plot(TrialTime,data_trials(:,oxychanL,i),'r',TrialTime,data_trials(:,dxychanL,i),'b');
    title(['Left Hemisphere condition ' num2str(allPerm(i))]);xline(0,'m');ylim([-1*Ylim Ylim]);
    subplot(2,1,2); 
    plot(TrialTime,data_trials(:,oxychanR,i),'r',TrialTime,data_trials(:,dxychanR,i),'b');
    title(['Right Hemisphere condition ' num2str(allPerm(i))]);xline(0,'m');ylim([-1*Ylim Ylim]);
    pause(2);
end

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

%% average across the middle 10 trials
data_trials10=data_trials(:,:,2:11);
data_trials10m=mean(data_trials10,3);
figure;
plot(TrialTime,data_trials10m(:,1:2:43),'r',TrialTime,data_trials10m(:,2:2:44),'b');
xline(10,'k','stimulus end');
xlabel('time (s)')

allPerm10=allPerm(2:11);


%% 2 subplots based on a combined time series
% filtered_data;filtered_oxyvals2;filtered_dxyvals2;
Fs=data_flt.fsample;
time=cell2mat(data_flt.time);
time(onsets(i))

% load stimulus data
% cd C:\Users\NIRS\Documents\zhibin\2022032801
% load('2022032801.mat'); allPerm;num2str(seed)

Ylim=1.5;
Xlim=[40000 60000];

figure('units','normalized','outerposition',[0 0 1 0.3]);
subplot(2,1,1);% left channels
% plot(time,filtered_data(:,1:2:21),'r');% oxy values on the left hemisphere
% hold on;
% plot(time,filtered_data(:,2:2:22),'b');% dxy values on the left hemisphere
% plot(time,filtered_oxyvals2(:,1:11),'r');hold on;
% plot(time,filtered_dxyvals2(:,1:11),'b');
plot(Combined_norm_data(:,1:2:21),'r');% oxy values on the left hemisphere
hold on;
plot(Combined_norm_data(:,2:2:22),'b'); % dxy values on the left hemisphere
% for i=1:length(allPerm)
%     if allPerm(i)==1
%         xline(time(onsets(2*i-1)),'r',{'L'})
%     else
%         xline(time(onsets(2*i-1)),'b',{'R'})
%     end
% end
for i=1:length(allPerm)
    if allPerm(i)==1
        xline(trial_length*(i-1)+1,'r',{'L'})
    else
        xline(trial_length*(i-1)+1,'b',{'R'})
    end
end
title('Left cortex');hold off;
% ylim([-Ylim Ylim]);
 xlim(Xlim);

subplot(2,1,2);% right channels
% plot(time,filtered_data(:,23:2:43),'r'); % oxy values on the left hemisphere
% hold on;
% plot(time,filtered_data(:,24:2:44),'b');  % dxy values on the left hemisphere
% plot(time,filtered_oxyvals2(:,12:22),'r');hold on;
% plot(time,filtered_dxyvals2(:,12:22),'b');
plot(Combined_norm_data(:,23:2:43),'r');% oxy values on the left hemisphere
hold on;
plot(Combined_norm_data(:,24:2:44),'b'); % dxy values on the left hemisphere
% for i=1:length(allPerm)
%     if allPerm(i)==1
%         xline(time(onsets(2*i-1)),'r',{'L'})
%     else
%         xline(time(onsets(2*i-1)),'b',{'R'})
%     end
% end
for i=1:length(allPerm)
    if allPerm(i)==1
        xline(trial_length*(i-1)+1,'r',{'L'})
    else
        xline(trial_length*(i-1)+1,'b',{'R'})
    end
end
title('Right cortex');hold off;
% ylim([-Ylim Ylim]);
 xlim(Xlim);


%% Four subplots
% ******************************************************************
cd C:\Users\NIRS\Documents\GitHub\Motor_cordination\fNIR_test
% plot 
figure;% keep all channels, but average across all trials
Chanselected=[1:6];
Trialselected=[30];
TrialTime=1/Fs*[1:trial_length]-10;
Ylim=3;
subplot(2,2,1);
TrialsLO2ChansL=data_trialsL(:,oxychanL(Chanselected),Trialselected);  
TrialsLHHChansL=data_trialsL(:,dxychanL(Chanselected),Trialselected);
plot(TrialTime,TrialsLO2ChansL,'r');hold on;
plot(TrialTime,TrialsLHHChansL,'b');hold off;grid;ylim([-Ylim Ylim]);
title('Left cortex - Left stim');

subplot(2,2,2);
TrialsLO2ChansR=data_trialsL(:,oxychanR(Chanselected),Trialselected); 
TrialsLHHChansR=data_trialsL(:,dxychanR(Chanselected),Trialselected);
plot(TrialTime,TrialsLO2ChansR,'r');hold on;
plot(TrialTime,TrialsLHHChansR,'b');hold off;grid;ylim([-Ylim Ylim]);
title('Right cortex - Left stim');

subplot(2,2,3);
TrialsRO2ChansL=data_trialsR(:,oxychanL(Chanselected),Trialselected);  
TrialsRHHChansL=data_trialsR(:,dxychanL(Chanselected),Trialselected); 
plotx(TrialTime,TrialsRO2ChansL,'r');hold on;
plotx(TrialTime,TrialsRHHChansL,'b');hold off;grid;ylim([-Ylim Ylim]);
title('Left cortex - Right stim');

subplot(2,2,4);
TrialsRO2ChansR=data_trialsR(:,oxychanR(Chanselected),Trialselected);  
TrialsRHHChansR=data_trialsR(:,dxychanR(Chanselected),Trialselected); 
plot(TrialTime,TrialsRO2ChansR,'r');hold on;
plot(TrialTime,TrialsRHHChansR,'b');hold off;grid;ylim([-Ylim Ylim]);
title('Right cortex - Right stim');

% suptitle(['channels ' num2str(dxychanR(Chanselected)) ' trials ' num2str(Trialselected)]);

%******************************************************************
figure;
Chanselected=[1:6];
subplot(2,2,1);
TrialsLO2ChansL=squeeze(mean(data_trialsL(:,oxychanL(Chanselected),:),2))'; % keep all trials, but average across all right oxy channels 
TrialsLHHChansL=squeeze(mean(data_trialsL(:,dxychanL(Chanselected),:),2))';
plotSEM(TrialTime,TrialsLO2ChansL,'r');hold on;
plotSEM(TrialTime,TrialsLHHChansL,'b');hold off;grid;ylim([-Ylim Ylim]);
title('Left cortex - Left stim');

subplot(2,2,2);
TrialsLO2ChansR=squeeze(mean(data_trialsL(:,oxychanR(Chanselected),:),2))'; % keep all trials, but average across all right oxy channels 
TrialsLHHChansR=squeeze(mean(data_trialsL(:,dxychanR(Chanselected),:),2))';
plotSEM(TrialTime,TrialsLO2ChansR,'r');hold on;
plotSEM(TrialTime,TrialsLHHChansR,'b');hold off;grid;ylim([-Ylim Ylim]);
title('Right cortex - Left stim');

subplot(2,2,3);
TrialsRO2ChansL=squeeze(mean(data_trialsR(:,oxychanL(Chanselected),:),2))'; % keep all trials, but average across all right oxy channels 
TrialsRHHChansL=squeeze(mean(data_trialsR(:,dxychanL(Chanselected),:),2))';
plotSEM(TrialTime,TrialsRO2ChansL,'r');hold on;
plotSEM(TrialTime,TrialsRHHChansL,'b');hold off;grid;ylim([-Ylim Ylim]);
title('Left cortex - Right stim');

subplot(2,2,4);
TrialsRO2ChansR=squeeze(mean(data_trialsR(:,oxychanR(Chanselected),:),2))'; % keep all trials, but average across all right oxy channels 
TrialsRHHChansR=squeeze(mean(data_trialsR(:,dxychanR(Chanselected),:),2))';
plotSEM(TrialTime,TrialsRO2ChansR,'r');hold on;
plotSEM(TrialTime,TrialsRHHChansR,'b');hold off;grid;ylim([-Ylim Ylim]);
title('Right cortex - Right stim');












%% load xdf from LSL 
% check out the script analysis_on20220415xdf.m 
cd /home/zhibin/Documents/GitHub/Motor_cordination/fNIR_test
open analysis_on20220415xdf.m

% cd C:\Users\NIRS\Documents\zhibin\2022032801
cd /home/zhibin/Documents/Artinis_NIRS/zhibin/20220410
fileList = dir('*.xdf');
fileList.name
[streams,fileheader] = load_xdf(fileList.name);
[streams,fileheader] = load_xdf('20220415.xdf');

info=streams{1,2}.info;
info.name % to comfirm it is the stream of OxySoft
time_series=streams{1,2}.time_series;
time_stamps=streams{1,2}.time_stamps;
segments=streams{1,2}.segments;
srate=info.effective_srate;
clock_offsets=info.clock_offsets;

% open LSL mapping and find the oxy doxy channels
open 20220415.xlsx

figure('units','normalized','outerposition',[0 0 1 0.6]);
plotx(time_stamps',time_series(17:2:59,:)','r');
hold on;
plotx(time_stamps',time_series(18:2:60,:)','b');
xlabel('time-stamps from LSL')
xlabel('time')

% examine the time of ambient light noise
signal3=sort(time_series(3,:));
find(time_series(3,:)==signal3(end))
% examine werid timestamps (don't use LSL time_stamps)
plot(time_stamps','r.') % time stamps values is high for a short time first then dropped

% event marker stream from LSL
info=streams{1,1}.info;
info.name % confirm that it is OxySoft Event Marker
time_series2=streams{1,1}.time_series; % noting here
time_stamps2=streams{1,1}.time_stamps; % noting here

% plot Markers on top the signal (not work on event markers from LSL in OxySoft Event Markers)
hold on;
for i=1:length(time_stamps2)
    xline(time_stamps2(i),'m');
end
hold off;

% event marker stream from LSL
info=streams{1,3}.info;
info.name % confirm that it is MyMarkerStream
time_series3=streams{1,3}.time_series; % C C A B A B A B ...
time_stamps3=streams{1,3}.time_stamps; % 

% plot Markers on top the signal (work on event markers from LSL in MyMarkerStream)
hold on;
for i=1:length(time_stamps3)
    xline(time_stamps3(i),'m');
end
hold off;

