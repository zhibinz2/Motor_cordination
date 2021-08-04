%% This code plot the scalp and source Coherence for all subjects in one 200 ms window
%% In 50 ms steps with 150ms overlay
%%
cd 'G:\concurrent_TMS_EEG'
clear; clc; 
eeglab; %load eeglab files, so that frequency filter function will work properly later
close all;

fs=2048;
fileName=[17 26 32 34 35 36 41 71 73 74 75 79 80 81 87 88]; % index for each subject's data file
NumSub=length(fileName); % number of subjects
conditions=[3 6]; % conditions to analyze

% select 200 ms interval before TMS as the baseline
stim=1045;% index at the time point of TMS
sample_start_prestim=stim-round(0.4*fs);
interval_50ms=round(50/1000*fs);
interval_200ms=round(200/1000*fs);
prestim=[sample_start_prestim sample_start_prestim+interval_200ms];

% select 200 ms windows post TMS with 150 ms overlap
t_start=0.013; % Set 13 ms after TMS as the starting point
sample_start_postim=stim+round(t_start*fs);
nWin=25; % number of the 200 ms windows to investigate
poststim=[sample_start_postim+[1:interval_50ms:(interval_50ms*(nWin-1)+1)]' ...
    sample_start_postim+interval_200ms+[1:interval_50ms:(interval_50ms*(nWin-1)+1)]']; % create an array of the windows

dataEnd=sample_start_postim+interval_200ms+(interval_50ms*(nWin-1)+1)+200; % add 200 at the end for the edge artifact to be cut out later after filter
win =1:dataEnd;
maxfreq = 50; 
load('TMSEEG_headmodel.mat');

% return

%% plot all subjects in for loop

for s=1:NumSub %s=14; % select a subject
    subject=fileName(s);
    [cond,time_stamp,channel_info,labels] = loadsubject(subject);
    window_stamp=time_stamp(win);
        for c=1:2 % c=2; % select a conditon
            condition=conditions(c); 
            condition_data=cell2mat(cond(condition));
            condition_data=condition_data(:,win,:);
 
%			This function pass condition_data and get back artifact and goodepochs and goodchanels.
            [filtered_data,artifact,goodchans,goodepochs] = Clean_artifact(condition_data,window_stamp);
%			This function pass filtered_data and get back average referenced data.
            [reRef_data] = reRef(filtered_data,goodchans);
            % erpdata = mean(reRef_data(:,:,goodepochs),3); % This includes only good epochs
            % plotx(window_stamp,erpdata(:,goodchans));% plot only the good channels
            % timeEnd=time_stamp(dataEnd)
            % sampleEnd=time_stamp(sample_start_postim+interval_200ms+(interval_50ms*(nWin-1)+1))

            
%           to do source localization
%             gainmat=Gain(goodchans,:); %you need to identify the good channels for each run.  
%           WE NEED TO DO WEIGHTED L2 NORM INVERSE LATER
			%gainmat = Gain(goodchan,:)./(ones(62,1)*sqrt(sum(abs(Gain(goodchan,:)).^2,1))
            gainmat = Gain(goodchans,:)./((ones(length(goodchans),1)*sqrt(sum(abs(Gain(goodchans,:)).^2,1))));
            regmethod='prctile';
            regparam=15;
            [inversemat, stat, reconstructed] = inversemodel(gainmat,regmethod,regparam);
            
            % project to the source
            for m=1:size(reRef_data,3)
                source_temp_trial=reRef_data(:,goodchans,m)*inversemat';
                source_data(:,:,m)=source_temp_trial;
            end
            
            % allspectra the scalpe electrode before stim
            win_prestim=prestim(1):prestim(2);
            [pow,freqs,df,eppow,corr,cprod,fcoef] = allspectra(reRef_data,fs,maxfreq,goodepochs,win_prestim);
            coh_before_electrode{s,c} = abs(corr).^2;
            pow_before_electrode{s,c} = pow;
            % allspectra the source coherence before stim
            [pow,freqs,df,eppow, corr,cprod,sourceorient,fcoefsource] = allspectrasource(source_data,fs,maxfreq,goodepochs,win_prestim);  
            coh_before_source{s,c} =abs(corr).^2;
            pow_before_source{s,c} = pow; 	
            
            % allspectra post stim
            for w = 1:size(poststim,1) 
                win_poststim=poststim(w,1):poststim(w,2);
                % allspectra the electrode before stim
            	[pow,freqs,df,eppow,corr,cprod,fcoef] = allspectra(reRef_data,fs,maxfreq,goodepochs,win_poststim);
                coh_poststim_electrode{s,c,w} = abs(corr).^2;
                pow_poststim_electrode{s,c,w} = pow;
                % allspectra the source post stim
                [pow,freqs,df,eppow, corr,cprod,sourceorient,fcoefsource] = allspectrasource(source_data,fs,maxfreq,goodepochs,win_poststim);  
                coh_poststim_source{s,c,w} = abs(corr).^2;
                pow_poststim_source{s,c,w} = pow;
        	end;
            
        end
end

%%
% run sum_average_change
% run plot_coh_average
% run plot_coh_change
% run plot_coh_over_time

