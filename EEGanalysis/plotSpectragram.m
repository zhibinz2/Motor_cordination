%% BEFORE Laplacian data check
data_trials; 
% plotx(data_trials(:,1:128,1));
data_trials_EEG;
% plotx(data_trials_EEG(:,:,1));

%% After Laplacian data check
filtered_data=filtered_data(1:dataEnd+1000,:); % pad one second

% plotx(filtered_data);hold on;xline(dataEnd,'m--');
% plotx(filtered_data(1000:2000,1));

%% 
% So the laplacian is in microvolts/mm^2 right now.    
% more typical in published literature would be microvolts/cm^2.   
% in which case you should multiply by 100 the raw Laplacian values before FFT
% ok multiply by 100 and you should be fine.
filtered_data=filtered_data*100;
% figure; plotx(filtered_data(1000:2000,1));

cd /home/zhibin/Documents/GitHub/Motor_cordination/EEGanalysis/20211102
clear data_trials
run integrate_EEG_into_data_trials_step3.m
data_trials % check whether included the baseline



%% open baselinecorrect

%% open wavelet

wfreq = 1:50; % frequency(s) of interest
wfc = 1.5;
sr=1000;
cnorm = wavelet(filtered_data(1001:5000),sr,wfc,wfreq);

imagesc(cnorm);

figure;imagesc(abs(cnorm))
colormap jet
colorbar 
caxis 
xlabel('time(ms)');
ylabel('frequencies(Hz)');
set(gca,'ydir','normal')


%% Jenny's code:  https://github.com/rameshsrinivasanuci/matlab/blob/master/jenny/WaveletTransform.m
% now let's recover the time course
wfreq = 1:50; % frequency(s) of interest
wfc = 3;
cnorm = wavelet(erp,sr,wfc,wfreq);

timecourse = complex(zeros(length(wfreq),4000));
for i = 1:4000
    timepoint = squeeze(cnorm(:,i,:));
    timesignal = timepoint * weights;
    timecourse(:,i) = timesignal;
end

% smooth it, fit a line to smooth the spectrum 
sgftimecourse = sgolayfilt(abs(timecourse.'),1,31);


% now let's loop through individual trial
trialtimecourse = complex(zeros(4000,length(wfreq),360));
for i = 1:length(finalgoodtrials)
        trialdata = data(:,:,finalgoodtrials(i));
        trialwavelet = wavelet(trialdata,sr,wfc,wfreq);
    for j = 1:length(wfreq)
        trialfreq = squeeze(trialwavelet(j,:,:));
        trialsignal = trialfreq * weights;
        trialtimecourse(:,j,i) = trialsignal;
    end
end

% extract the power at the target frequency


% making plots
figure;imagesc(abs(sgftimecourse)')
colormap jet
colorbar 
caxis ([0 2])
xlabel('time(s)');
ylabel('frequencies(Hz)');

figure;imagesc(abs(timecourse))
colormap jet
colorbar 
caxis 
xlabel('time(s)');
ylabel('frequencies(Hz)');


% smooth it, fit a line to smooth the spectrum 
sgftimecourse = sgolayfilt(abs(timecourse.'),1,31);
