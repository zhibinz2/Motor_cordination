clear
cd C:\Users\NIRS\Documents\zhibin\2022032801

%% single channel demo: load data with fieldtrip
cfg = [];
cfg.dataset = '2022032801.oxy3'; % only oxy3 save the events 
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
cfg.dpf = 5.9;
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
cfg             = [];
cfg.dataset = '2022032801.oxy3';
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

% (high pass filtering, throw away very low-frequency) (work! only when > 0.05 Hz)
cfg                 = [];
cfg.hpfilter        = 'yes';
cfg.hpfreq          = 0.05;
data_hpf            = ft_preprocessing(cfg,data_raw);

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

%% (0.05~0.8 Hz band-pass filter: >0.05 high pass filt on <0.8 data_lpf)
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
%% plot singal from oxy3.mat to examine ylim
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
% high pass (no paddings needed? for now)
Hd = makefilter(Fs,0.05,0.04,6,20,1); xlim([0 1]);
which filtfilt.m -all
cd 'C:\Program Files\MATLAB\R2021b\toolbox\signal\signal\'
which filtfilt.m
filtered_oxyvals1=filtfilthd(Hd,detrend_oxyvals);
filtered_dxyvals1=filtfilthd(Hd,detrend_dxyvals);
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(time,filtered_oxyvals1,'r',time,filtered_dxyvals1,'b');title('oxyval & dxyval > 0.05Hz'); 
ylim([-20 20]);ylim([-3 3]);
% low pass (this will create short edge artifact, no padding for now)
Hd = makefilter(Fs,0.8,0.9,6,20,1);
which filtfilt.m -all
cd 'C:\Program Files\MATLAB\R2021b\toolbox\signal\signal\'
which filtfilt.m
filtered_oxyvals2=filtfilthd(Hd,filtered_oxyvals1);
filtered_dxyvals2=filtfilthd(Hd,filtered_dxyvals1);
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(time,filtered_oxyvals2,'r',time,filtered_dxyvals2,'b');title('oxyval & dxyval 0.05~0.8 Hz'); 
ylim([-20 20]);ylim([-1.5 1.5]);

% channel position
chanlabel=nirs_data.label;
chanPos=nirs_data.chanPos;
plot(chanPos(1,:),chanPos(2,:),'r.');
text(chanPos(1,:),chanPos(2,:),chanlabel,'VerticalAlignment','bottom','HorizontalAlignment','right');
ylim([-7 4]);xlim([-2 13]);

% AD channel
ADvalues=nirs_data.ADvalues;
ADlabel=nirs_data.ADlabel;
for i=1:15
    subplot(4,4,i);
    plot(time, ADvalues(:,i));
    title(ADlabel{i});
end
%% use event onset in nirs_data from oxy3.mat
events=nirs_data.events;
onsets=events.onsets{1,1};

figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(time, ADvalues(:,15),'b');
hold on;
for i=1:length(onsets)
    plot(time(onsets(i)),10000*ones(length(onsets),1),'r.');
end
hold off;
legend('photocell signal','event onsets generated by Oxysoft');

%% load xdf from LSL
cd C:\Users\NIRS\Documents\zhibin\2022032801
fileList = dir('*.xdf');
fileList.name
[streams,fileheader] = load_xdf(fileList.name);

info=streams{1,2}.info;
info.name % to comfirm it is the stream of OxySoft
time_series=streams{1,2}.time_series;
time_stamps=streams{1,2}.time_stamps;
segments=streams{1,2}.segments;
srate=info.effective_srate;
clock_offsets=info.clock_offsets;

figure('units','normalized','outerposition',[0 0 1 0.6]);
plotx(time_stamps',time_series(3:2:45,:)','r');
hold on;
plotx(time_stamps',time_series(4:2:46,:)','b');
xlabel('time-stamps from LSL')
xlabel('time')

% examine the time of ambient light noise
signal3=sort(time_series(3,:));
find(time_series(3,:)==signal3(end))
% examine werid timestamps (don't use LSL time_stamps)
plot(time_stamps','r.') % time stamps values is high for a short time first then dropped

% event marker stream from LSL
info=streams{1,1}.info;
info.name % confirm that it is OxySoft Event Marker / MyMarkerStream
time_series2=streams{1,1}.time_series; % noting here
time_stamps2=streams{1,1}.time_stamps; % noting here


% event marker stream from LSL
info=streams{1,3}.info;
info.name % confirm that it is MyMarkerStream
time_series3=streams{1,3}.time_series; % C C A B A B A B ...
time_stamps3=streams{1,3}.time_stamps; % 