%% analyse TMSi recorded EEG on the checkerboard experiment
% to find lateralization in TMSi EEG
cd /home/zhibin/Documents/Artinis_NIRS/zhibin/20220417
cd /home/zhibin/Documents/Artinis_NIRS/zhibin/2022042105


%% load TMSi data
clear;close all;
% Step 1: Open Poly5 file.
addpath /home/zhibin/Documents/GitHub/TMSi_SAGA_matlab % hnlb
d = TMSiSAGA.Poly5.read('/home/zhibin/Documents/Artinis_NIRS/zhibin/20220417/2022041704-20220417T111613.DATA.Poly5');
d = TMSiSAGA.Poly5.read('/home/zhibin/Documents/Artinis_NIRS/zhibin/2022042105/2022042105.Poly5');


samples=d.samples;
sr=d.sample_rate;
channels=d.channels;

% Organize the channels
numbers=num2str([1:length(channels)]');
labels=strings(length(channels),1);
units=strings(length(channels),1);
for i=1:length(channels)
    labels(i)=channels{i}.alternative_name;
    units(i)=channels{i}.unit_name;
end
channels_info=table(numbers,labels,units)

% Create time stamps
num2str(d.time)
time=[1/sr:1/sr:d.time]';

% Plot channels of Key presses, photocells, EMG
% look for the second ISO aux channel for the photocell 
ISOauxind=find(labels=='ISO aux');
plot(samples(ISOauxind(1),:),'b'); %  ISO aux = analog
plot(samples(ISOauxind(2),:),'b'); %  ISO aux = analog
plot(-1*samples(ISOauxind(2),:),'b'); %  ISO aux = analog
ISOauxind=ISOauxind(1);
ISOauxind=ISOauxind(2);
% Save Channel of photocells
Photocell=samples(ISOauxind,:)';
Photocell=-1*samples(ISOauxind,:)';
% Save Channel of EEG
EEG=samples(1:32,:)';
numChan=32;
% View EEG
plot(time,EEG);

%% load TMSi for EEGLAB
[file, pathname] = uigetfile({'*.Poly5';'*.S00';'*.TMS32'},'Pick your file');
[path,filename,extension] = fileparts(file);
% open Poly5toEEGlab

% Step 2: Plot a single channel.
plot((0:(d.num_samples - 1)) / d.sample_rate, d.samples(2, :));

% Step 3: Save dataset in the same directory as the *.Poly5 file.
% Open EEGlab
eeglab

% Load TMSi's channel location file
load('EEGChannels64TMSi.mat', 'ChanLocs');

% Transform data to eeglab. 
eegdataset = toEEGLab(d, ChanLocs);

% Save dataset in the same directory as the *.Poly5 file.
pop_saveset(eegdataset,'filename',filename,'filepath',pathname)
disp(['Data saved as EEGlab dataset (.set) in this folder: ',pathname])

%% load stimulus data
cd /home/zhibin/Documents/Artinis_NIRS/zhibin/20220417
load('2022041703.mat')

%% Covert Photocell signals to time points
time;

% for photocell
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
Halfhigh1=1/4*(max(Photocell)-min(Photocell)); % value for 'MinPeakProminence'
% Check if need to adjust the Halfhigh cutoff
close;figure;
findpeaks(Photocell,time,'MinPeakProminence',Halfhigh1,'Annotate','extents');
yline(Halfhigh1,'m','MinPeakProminence');
% ylim([0 5]);
ylim([-5 5]);

% locate the trial sessions % pks=value of the peak % locs=time of the peak
[pks,locs] = findpeaks(Photocell,time,'MinPeakProminence',Halfhigh1,'Annotate','extents');

% examine locs on top of signal plot
EEG=samples(1:32,:)';
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(time,EEG); 
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
plot(EEG); 
hold on;
for i=1:length(EventIndices)
    xline(EventIndices(i),'m');
end

%% extract the conditions 
load('2022041703.mat')
allPerm;
numTrials=12;
allPerm=allPerm(1:numTrials);

%% EEG process
% cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis
% open EEG_preprocessing.m
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(EEG);
[x, y] = ginput(2); % read two mouse clicks on the plot % x were index, y were real values
% get the proximate index
string(x)
ind1=round(x(1))
ind2=round(x(2))
xline(ind1,'b');xline(ind2,'b');
hold off;
% cut off artifact at the beginning and end, then extended the two sides to the original length
EEG2=[repmat(EEG(ind1,:),ind1-1,1); EEG(ind1:ind2,:); repmat(EEG(ind2,:),length(time)-ind2,1)];
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(EEG2);
%% detrend the EEG data (no padding needed)
detrend_data=detrend(EEG2,2); % even better (When n = 2, detrend removes the quadratic trend.)
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(detrend_data);
%% filtfilthd method (hnl) high pass first
which filtfilt -all;
cd /usr/local/MATLAB/R2022a/toolbox/signal/signal/

% high pass (no paddings needed)
Hd = makefilter(sr,0.2,0.15,6,20,0); 
filtered_data1=filtfilthd(Hd,detrend_data);
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(filtered_data1);

% add padding
padding=zeros(round(size(filtered_data1,1)/10), size(filtered_data1,2));
filtered_data2=cat(1,padding,filtered_data1,padding);
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(filtered_data2);

% low pass 
% (this will create short edge artifact)
% (if added zero paddings, edge artifact disappear)
% (remove existing filtered_data variable from workspace might fasten)
tic
filtered_data3=[];
Hd = makefilter(sr,50,51,6,20,0);  
filtered_data3=filtfilthd(Hd,filtered_data2);
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(filtered_data3);
% ylim([-100 100]);
toc

% remove padding
filtered_data4=filtered_data3((size(padding,1)+1):(size(padding,1)+size(detrend_data,1)),:);
figure('units','normalized','outerposition',[0 0 1 0.3]);
plotx(filtered_data4);
% ylim([-100 100]);

clearvars samples filtered_data1 filtered_data2 filtered_data3 d EEG2 EEG padding

%% examine EventIndices on top of signal
figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(filtered_data4); 
hold on;
for i=1:length(EventIndices)
    xline(EventIndices(i),'m');
end

%% ICA
% run ICA
[icasig, A, W] = fastica(filtered_data4');
% Plot ICA component
for i=1:size(A,2)
    SqA(i)=sumsqr(A(:,i));
end
figure;
plot(1:size(A,2),SqA,'ro');ylabel('sum of square of column in A');xlabel('ICs');
[B,I]=sort(SqA,'descend');

% topoplot to examine them
cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/channels_info
load('chaninfo.mat')
ComponentsExam=I(1:10);
figure;
for i=1:length(ComponentsExam)
    subplot(5,2,i);
    topoplot(A(:,ComponentsExam(i)),chaninfo,'nosedir','+X');title(['component' num2str(ComponentsExam(i))]);colorbar;
end

% check amplitude for compoent 10 and 15
figure;
plot(icasig(10,:)','r');hold on;plot(icasig(15,:)','b');legend('component 5','component 15');
hold on;
indtemp=find(allPerm==UniCondi(1));
xline(EventIndices(2+2*indtemp-1),'m',{'L'});
indtemp=find(allPerm==UniCondi(2));
xline(EventIndices(2+2*indtemp-1),'m',{'R'});
yline(0,'m')
hold off;

% Display signal before ICA 
figure('units','normalized','outerposition',[0 0 1 0.5]);
plot(time,filtered_data4);
hold on; hold off;
title('EEG Singal Before ICA');

% Display ICA signals
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,1,1);
plot(time,icasig');
hold on;hold off;
title('all IC time series');ylabel('uV');
subplot(2,1,2);
plot(time,icasig(ComponentsExam,:)');
hold on;hold off;
legend(strsplit(num2str(ComponentsExam)));
title('first 10 IC time series');xlabel('time samples');ylabel('uV');

% Calculate Correlation
% FP1 and FP2 are channel 1 and 3;
% compute correlation between FP1 and FP2;
[RHO1,PVAL1] = corr(filtered_data4(:,1),icasig');
figure
subplot(1,2,1);
plot(1:length(RHO1),RHO1,'bo');hold on;
yyaxis left;ylabel('correlation coefficient','color','b');set(gca,'ycolor','b');
plot(1:length(RHO1),PVAL1,'ro');
yyaxis right; ylabel('p value','color','r');set(gca,'ycolor','r'); xlabel('ICs');
legend({'correlation coefficient','p-values '});title('correlation with FP1');hold off;
subplot(1,2,2);
[RHO3,PVAL3] = corr(filtered_data4(:,3),icasig');
plot(1:length(RHO1),RHO3,'bo');hold on;
yyaxis left;ylabel('correlation coefficient','color','b');set(gca,'ycolor','b');
plot(1:length(RHO1),PVAL3,'ro');
yyaxis right; ylabel('p value','color','r');set(gca,'ycolor','r'); xlabel('ICs');
legend({'correlation coefficient','p-values '});title('correlation with FP2');hold off;

[B1,I1]=sort(abs(RHO1),'descend');[B3,I3]=sort(abs(RHO3),'descend');
ComponentsExam=unique([I1(1) I3(1)]);
figure;
for i=1:length(ComponentsExam)
    subplot(length(ComponentsExam),1,i);
    topoplot(A(:,ComponentsExam(i)),chaninfo,'nosedir','+X');title(['component' num2str(ComponentsExam(i))]);colorbar;
end

%% run ICA separately on the two conditions
% cut data into trials and reconnect them back according to condition
%  organize EEG into trials using EventIndices
numTrials=12;
numChans=32;
trial_length=round(sr*3); % 2 rest + 2s stimulus + 2 rest 

TrialTime=1/sr*[1:trial_length];
data_trials=zeros(trial_length,numChans,numTrials);
for i=1:numTrials % i=numTrials
    data_trials(:,:,i)=filtered_data4(EventIndices(2+2*i-1):(EventIndices(2+2*i-1)+round(sr*3)-1),:);
end
% separate the two condition
% conditions
allPerm;
EventIndices;
UniCondi=unique(allPerm);
indtemp=find(allPerm==UniCondi(1));
datatrials1=data_trials(:,:,indtemp);
indtemp=find(allPerm==UniCondi(2));
datatrials2=data_trials(:,:,indtemp);
% connect two set of data separately
% Fast Method
split_data = num2cell(datatrials1, [1 2]); %split A keeping dimension 1 and 2 intact
datatrials1connected=vertcat(split_data{:});
split_data = num2cell(datatrials2, [1 2]); %split A keeping dimension 1 and 2 intact
datatrials2connected=vertcat(split_data{:});

% run ICA separtely
[icasig, A, W] = fastica(datatrials1connected');
[icasig, A, W] = fastica(datatrials2connected');
% Plot ICA component
for i=1:size(A,2)
    SqA(i)=sumsqr(A(:,i));
end
figure;
plot(1:size(A,2),SqA,'ro');ylabel('sum of square of column in A');xlabel('ICs');
[B,I]=sort(SqA,'descend');

% topoplot to examine them
cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/channels_info
load('chaninfo.mat')
ComponentsExam=I(1:10);
figure;
for i=1:length(ComponentsExam)
    subplot(5,2,i);
    topoplot(A(:,ComponentsExam(i)),chaninfo,'nosedir','+X');title(['component' num2str(ComponentsExam(i))]);colorbar;
end


%% Deside which components to remove and mix back the signal and display
% ComponentRemove=I(1);
ComponentRemove=I(2);

Abackup=A;icasigbackup=icasig;
% A=Abackup;icasig=icasigbackup;
A(:,ComponentRemove)=0; icasig(ComponentRemove,:)=0;

mixedsig=A*icasig;

% Plot before and after
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,1,1);
% before ICA
plotx(time,filtered_data4);
hold on;hold off;
title('EEG Singal Before ICA');
% after ICA
subplot(2,1,2);
plotx(time,mixedsig);
hold on;hold off;
title('Mixed Signal with ICs removed');

%% Only Keep the two visual component 
ComponentRemove=I([2 4:32]);

Abackup=A;icasigbackup=icasig;
% A=Abackup;icasig=icasigbackup;
A(:,ComponentRemove)=0; icasig(ComponentRemove,:)=0;

mixedsig=A*icasig;

% Plot before and after
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,1,1);
% before ICA
plotx(time,filtered_data4);
hold on;hold off;
title('EEG Singal Before ICA');
% after ICA
subplot(2,1,2);
plotx(time,mixedsig);
hold on;hold off;
title('Mixed Signal with ICs removed');
%% examine
close;figure;
plot(mixedsig');
title('Mixed Signal with ICs removed');
hold on;
for i=1:length(EventIndices)
    xline(EventIndices(i),'m');
end

%% organize EEG into trials using EventIndices
mixedsig=mixedsig';
numTrials=12;
numChans=32;
% trial_length=round(sr*13); % 3s baseline + 10s stimulus 
% trial_length=round(sr*30); % 10s stimulus + 20 rest 
trial_length=round(sr*6); % 2 rest + 2s stimulus + 2 rest 

time(EventIndices(1))
% TrialTime=1/sr*[1:trial_length];
TrialTime=1/sr*[1:trial_length]-2;
data_trials=zeros(trial_length,numChans,numTrials);
for i=1:numTrials % i=numTrials
%     data_trials(:,:,i)=mixedsig((EventIndices(2+2*i-1):(EventIndices(2+2*i-1)+round(sr*30)-1),:);
%     data_trials(:,:,i)=mixedsig(EventIndices(2+2*i-1):(EventIndices(2+2*i-1)+round(sr*30)-1),:);
    data_trials(:,:,i)=mixedsig((EventIndices(2+2*i-1)-round(sr*2)):(EventIndices(2+2*i-1)+round(sr*4)-1),:);
end
% examine
for i=1:numTrials % i=numTrials
    plot(TrialTime,data_trials(:,:,i));pause(1);
end 
%% baseline correction
% baselinesamps = 1:round(3*sr); % use the first 3s as baseline
% baselinesamps = 40000:60000;
baselinesamps = 1:round(2*sr); % use the first 5s as baseline

baselinecorrected_trial=zeros(size(data_trials)); % initialize this matric
% Then loop through each trial for baseline correction
for i=1:size(baselinecorrected_trial,3) % loop through trials
    trialdata=data_trials(:,:,i);
    newtrialdata = baselinecorrect(trialdata,baselinesamps);
    baselinecorrected_trial(:,:,i)=newtrialdata;
end

% examine
for i=1:numTrials % i=numTrials
    plot(TrialTime,baselinecorrected_trial(:,:,i));pause(1);
end 
%% artifact removal
data_trials_EEG=baselinecorrected_trial;clear baselinecorrected_trial;

nSamps = size(data_trials_EEG,1);
nChans = size(data_trials_EEG,2);
nEpochs = size(data_trials_EEG,3);

% find the the empty channels
sumovertime=sum(data_trials_EEG,1);
sumovertrials=sum(sumovertime,3);
zerochans = find(sumovertrials== 0);

% design a matrix with highest value for each channel on each trial(each single channel in each epoch). 
trialdatamax = squeeze(max(abs(data_trials_EEG),[],1));  
% examine the maximal abs trial matrix
% figure;imagesc(trialdatamax);colorbar;ylabel('chan');xlabel('trial');
% figure;imagesc(log10(trialdatamax));colorbar;ylabel('chan(log)');xlabel('trial(log)');


% first round of artifact labelling **************#########################
threshold=100; % set a amplitude threshold as criteria to remove bad channels and epochs later

% Create a label matrix that labels bad chan and epoch as 1 
binarythreshhold = zeros(size(trialdatamax));
binarythreshhold(trialdatamax > threshold) = 1;
% figure; imagesc(binarythreshhold); colorbar; ylabel('chan');xlabel('trial');
% title(['sum of Ones: ' num2str(sum(binarythreshhold,'all'))]);

chansum = sum(binarythreshhold,2); % how many bad epochs in each chan
epochsum = sum(binarythreshhold,1); % how many bad chans in each epoch

% bar(chansum); % examine the chans 
badchancrit = prctile(chansum,90); % criteria to determine bad chans
badchanedit = find(chansum > badchancrit); % get the initial bad chans
goodchanedit = setdiff([1:nChans],badchanedit); % get the initial good chans

% bar(epochsum); % examine the epochs
badepochcrit = prctile(epochsum,85); % criteria to determine bad epochs
badepochedit = find(epochsum > badepochcrit); % get the initial bad epochs
goodepochedit = setdiff([1:nEpochs],badepochedit); % get the initial good epochs

% second round of labelling **************#################################
threshold2=50;
% update the label matrix that labels bad chan and epoch as 1 
binarythreshhold = zeros(size(trialdatamax));
binarythreshhold(trialdatamax > threshold2) = 1;
% figure; imagesc(binarythreshhold); colorbar; ylabel('chan');xlabel('trial');
% figure; imagesc(binarythreshhold(goodchanedit,goodepochedit));colorbar; ylabel('chan');xlabel('trial');

epochsum2 = sum(binarythreshhold(goodchanedit,goodepochedit),1); % how many bad chans in each good epoch of the upated label matrix % whos % to look at the workspace
% figure;bar(goodepochedit,epochsum2); xlabel('good trials'); ylabel('number of bad channels');
newbadepochcrit = prctile(epochsum2,80); % criteria to determine bad epochs in 2nd round
newbadepochedit = goodepochedit(find(epochsum2 > newbadepochcrit)); % get the bad epochs in 2nd round
badepochedit = [badepochedit newbadepochedit]; % add the bad epochs in 1s round
goodepochedit_final = setdiff([1:nEpochs],badepochedit);  % get the final good epochs

chansum2 = sum(binarythreshhold(goodchanedit,goodepochedit),2); % how many bad epochs in each good chan of the upated label matrix
% figure;bar(goodchanedit,chansum2);xlabel('good channels');ylabel('number of bad trials');
newbadchancrit = prctile(chansum2, 90); % criteria to determine bad chans in 2nd round
newbadchanedit = goodchanedit(find(chansum2 > newbadchancrit)); % get the bad chans in 2nd round
badchanedit = [badchanedit; newbadchanedit']; % add the bad chans in 1s round
goodchanedit_final = setdiff([1:nChans],badchanedit);  % get the final good chans

% examine the updated label matrix after the 2nd round
% figure;imagesc(trialdatamax(goodchanedit_final,goodepochedit_final));colorbar;
% figure;imagesc(log10(trialdatamax(goodchanedit_final,goodepochedit_final)));colorbar;
% compare with previous one round artifact labelling method
% figure;imagesc(log10(trialdatamax(goodchans,goodepochs)));colorbar;
% figure;imagesc(log10(trialdatamax));colorbar;

badchans = badchanedit;
badepochs = badepochedit;
goodchans = goodchanedit_final;
goodepochs = goodepochedit_final;

artifact.badchans = badchans;
artifact.badepochs = badepochs;
artifact.goodchans = goodchans;
artifact.goodepochs = goodepochs;
artifact.editmatrix = trialdatamax;
artifact.threshhold = threshold2;
artifact.zerochans = zerochans;

%% Re reference
% This function pass filtered_data and get back average referenced data.
cd /home/zhibin/Documents/GitHub/Motor_cordination/EEGanalysis/20211102
[reRef_data] = reRef(data_trials_EEG,goodchans);
% examine
for i=1:numTrials % i=numTrials
    plot(TrialTime,reRef_data(:,goodchans,i)); 
    ylim([-60 60]);
    pause(1);
end 

%% Remove bad epochs and connect them back to one time serie and plot only the good channels
% Fast Method
split_reRef_data = num2cell(reRef_data(:,:,goodepochs), [1 2]); %split A keeping dimension 1 and 2 intact
reRef_data_goodepochs=vertcat(split_reRef_data{:});
% figure;plot(filtered_data(1:dataEnd,:));
figure;plot(reRef_data_goodepochs(:,goodchans));
title(['Number of good trials: ' num2str(length(goodepochs))]);

%% Baseline normalization
% baselinesamps = 1:round(2*sr); % use the first 2 s as baseline
basenormalized_trial=zeros(size(reRef_data)); % initialize this matric
% Then loop through each trial for baseline correction
for i=1:size(reRef_data,3) % loop through trials
    trialdata=reRef_data(:,:,i); % plot(trialdata)
    baselineMeanAbs=mean(abs(trialdata(baselinesamps,:)),1);
    newtrialdata = trialdata./(baselineMeanAbs.*ones(size(trialdata,1),1));
    basenormalized_trial(:,:,i)=newtrialdata;
end
% Plot one trial just to examine
subplot(2,1,1);plotx(TrialTime,reRef_data(:,:,1));title('before normalization');
subplot(2,1,2);plotx(TrialTime,basenormalized_trial(:,:,1));title('baseline normalized');

%% view the channel locations
cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis
open channels_info.m

%% conditions
allPerm;
EventIndices;
UniCondi=unique(allPerm);

%% ERP
maxfreq=30;
win=sr*2:sr*2.5;
Ymax=15;

allPermGood=allPerm(goodepochs);

for u=[1 2]; % sp=1:2
    indtemp=find(allPermGood==UniCondi(u));
    basenormalized_ERP=mean(data_trials(:,:,indtemp),3);
    figure('units','normalized','outerposition',[0 0 0.8 0.4]);
    for chan=1:32
        subplot('Position',[XXPLOT(chan) YYPLOT(chan) 0.05 0.05]); % not showing, why
        plot(TrialTime(win),basenormalized_ERP(win,chan));
        hold on
        if ~isempty(find([1:30 32]==chan))
        set(gca,'XTick',[]); set(gca,'YTick',[]); 
        end
        if chan==31
            xlabel('time (s)');ylabel('uV');
        end
        if any(goodchans(:) == chan) % Display good chan in green
            title([labels{chan}],'Color','green');
        else
            title([labels{chan}],'Color','red'); % Display bad chan in red
        end
% 
        ylim([-1*Ymax Ymax]);
%         xlim([1 25]);
        xline(0,'m');yline(0,'m');
        
    end
    %  sgtitle(['condition ' num2str(u)],'Color','red');
    suptitle(['condition ' num2str(u)]);
end



%% Get average spectra and Plot on scalp map 
maxfreq=30;
win=sr*2:sr*2.3;
Ymax=8;

allPermGood=allPerm(goodepochs);

for u=[1 2]; % sp=1:2
    indtemp=find(allPermGood==UniCondi(u));
    [pow,freqs,df] = allspectra(data_trials(win,:,indtemp),sr,maxfreq);
    figure('units','normalized','outerposition',[0 0 0.8 0.4]);
    for chan=1:32
        subplot('Position',[XXPLOT(chan) YYPLOT(chan) 0.05 0.05]); % not showing, why
        plot(freqs,pow(:,chan));
        hold on
        if ~isempty(find([1:30 32]==chan))
        set(gca,'XTick',[]); set(gca,'YTick',[]); 
        end
        if chan==31
            xlabel('freq');ylabel('pow');
        end
        if any(goodchans(:) == chan) % Display good chan in green
            title([labels{chan}],'Color','green');
        else
            title([labels{chan}],'Color','red'); % Display bad chan in red
        end

        ylim([0 Ymax]);
        xlim([1 25]);
        xline(8,'m')
        
    end
    %  sgtitle(['condition ' num2str(u)],'Color','red');
    suptitle(['condition ' num2str(u)]);
end


%% Plot on scalp map for spectrogram
% wavelet settings
% frequency(s) of interest
wfreq = [2 4 6 8 10 14 18 24 30 40]; % wfreq = list of frequencies to obtain wavelet coefficients. 
wfreq = [6 7 8 9 10]; % wfreq = list of frequencies to obtain wavelet coefficients. 
wfc = 1.5;  % wfc = wavelet parameter - # of cycles to use, typically 1.5

% Padding 500 samples of zeros before and after baselinecorrected_trial
basenormalized_trial;
padding=zeros(500,size(basenormalized_trial,2),size(basenormalized_trial,3));
baselinecorrected_trial_paddings=cat(1,padding,basenormalized_trial,padding);

for u=[1 2]; % sp=1:2
    indtemp=find(allPermGood==UniCondi(u));
    % Compute ERP in all channels
    basedlinecorrected_ERP=mean(baselinecorrected_trial_paddings(:,:,indtemp),3);
    % plot(basedlinecorrected_ERP);

    % wavelet
    cnorm = wavelet(basedlinecorrected_ERP,sr,wfc,wfreq);
    Powcnorm = abs(cnorm).^2;

    % Power normalization
    logPowcorm=log10(Powcnorm);
    baselineMean=log10(mean(Powcnorm(:,baselinesamps,:),2));
    % one way
    % normPowcnorm = 10* (logPowcorm./baselineMean);
    % second way
    normPowcnorm = 10*(logPowcorm-(ones(1,size(logPowcorm,2),1).*baselineMean));
    
    % plot settings
    ColorLim=10;

    figure('units','normalized','outerposition',[0 0 0.8 0.4]);

    for chan=1:32
        subplot('Position',[XXPLOT(chan) YYPLOT(chan) 0.05 0.05]);
        % open sgolay; open sgolayfilt
        imagesc((sgolayfilt(squeeze(normPowcnorm(:,:,chan))',1,61))');
        colormap jet; 
        caxis([-1*ColorLim ColorLim]);
        set(gca,'ydir','normal');
        hold on;xline(2500,'k','linewidth',2);hold off;
        xlim([2500 3000]);
        if any(goodchans(:) == chan) % Display good chan in green
            title([labels{chan}],'Color','green');
        else
            title([labels{chan}],'Color','red'); % Display bad chan in red
        end
        if ~isempty(find([1:30 32]==chan))
            set(gca,'XTick',[]); set(gca,'YTick',[]); 
        end
        if chan==31
            xlabel('time(s)'); ylabel('frequencies(Hz)'); 
            yticks([1:length(wfreq)]);yticklabels({'2','4','6','8','10','14','18','24','30','40'});
            colorbar;
        end
    end
    suptitle(['condition ' num2str(u)]);
end




