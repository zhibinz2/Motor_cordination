% load EEG all as one
EEG=loadcurry([pwd '/test_new_photocell_positions_sr1k_20210916.cdt']);
dataEEG=EEG.data;
datatimes=EEG.times;
Fs=EEG.srate;
% chanlocs=EEG.chanlocs;
% chanlocs128=chanlocs(1:128);
load('chanlocs.mat');
plot(datatimes,dataEEG);

%%
event=EEG.event;
% urevent=EEG.urevent;
event=struct2cell(event);
event=squeeze(event);
event=cell2mat(event);

eventtype=event(1,:);
eventtime=event(2,:);
eventnum=event(3,:);

% Upper left photocell event (start and end of trial)
ULevent=1200002;
% Upper right photocell event (end of planning)
URevent=1200001;

% find the index 
IndexTrialStartEnd=find(eventtype==ULevent);
IndexTrialStartEndLatency=eventtime(IndexTrialStartEnd);
diff(IndexTrialStartEndLatency) % trial length -- inter-trial interval -- trial length --
IndexTrialStartLatency=IndexTrialStartEndLatency(1:2:end);
diff(IndexTrialStartLatency) % trial length between start
ind_trial_start % index in locs of findpeaks
trialstart % index in time or latency

IndexTrialEndLatency=IndexTrialStartEndLatency(2:2:end);
diff(IndexTrialEndLatency) % trial length between end
ind_trial_end % index in locs of findpeaks
trialend % index in time or latency


IndexBetweenPlanEnd=find(eventtype==URevent);
IndexBetweenPlanEndLatency=eventtime(IndexBetweenPlanEnd);
diff(IndexBetweenPlanEndLatency) % trial length + inter-trial interval



%%
% plot photocell analog signal
analog1=133; %analog2=134;trigger=139;
analog1data=dataEEG(133,:);
plot(datatimes,analog1data);

% plot fft series
F133=fft(analog1data);
A133=abs(F133);
N=length(analog1data);
fV=linspace(0,Fs/2,N/2+1);

plot(fV,A133(1:length(fV)));
xlim([0 70]);
xlabel('frequency');
ylabel('amplitude (uV)');
title('Spetra');
