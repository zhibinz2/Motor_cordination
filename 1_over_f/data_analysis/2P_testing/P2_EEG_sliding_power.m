cd /ssd/zhibin/1overf/20220610_2P/Segmented_data/1_50Hz_ICAautomized
clear

samplingrate = 2000;
winsize = 10*samplingrate; % 10 Second with df of 0.1Hz

tempEEG = load('/ssd/zhibin/1overf/20220610_2P/Segmented_data/1_50Hz_ICAautomized/EEG20220610.mat'); 
tempBP = load('/ssd/zhibin/1overf/20220610_2P/Segmented_data/1_50Hz_ICAautomized/BP20220610.mat'); 

% syncopation expt
conditionNames={'uncoupled' 'L-lead' 'R-lead' 'mutual-1.3Hz'};


data(1,1).EEG = tempEEG.EEGCondi1L;
data(1,2).EEG = tempEEG.EEGCondi1R;
data(2,1).EEG = tempEEG.EEGCondi2L;
data(2,2).EEG = tempEEG.EEGCondi2R;
data(3,1).EEG = tempEEG.EEGCondi3L;
data(3,2).EEG = tempEEG.EEGCondi3R;
data(4,1).EEG = tempEEG.EEGCondi4L;
data(4,2).EEG = tempEEG.EEGCondi4R;
data(1,1).BP = tempBP.BPCondi1L;
data(1,2).BP = tempBP.BPCondi1R;
data(2,1).BP = tempBP.BPCondi2L;
data(2,2).BP = tempBP.BPCondi2R;
data(3,1).BP = tempBP.BPCondi3L;
data(3,2).BP = tempBP.BPCondi3R;
data(4,1).BP = tempBP.BPCondi4L;
data(4,2).BP = tempBP.BPCondi4R;


for j = 1:4 % four conditions
    BP_L=data(j,1).BP; BP_R=data(j,2).BP;
    IndL=find(BP_L);IndR=find(BP_R);%extract BP indices
    IndL(end);IndR(end);% last BP index
    nwin(j) = min([floor(length(BP_L(1:IndL(end)))/winsize) floor(length(BP_R(1:IndR(end)))/winsize)]);
    eeg1 = zscore(data(j,1).EEG);
    eeg2 = zscore(data(j,2).EEG);    
    
    for k = 1:nwin(j)  % something is wrong with the segmentation of some of your files. 
        samples = (k-1)*winsize+1:k*winsize;
        % interval1= find(data(j,1).BP(samples)); 
        % interval2= find(data(j,2).BP(samples));
        % nx = min([length(interval1) length(interval2)]);
        % cor(j).BP(k) = corr(diff(interval1(1:nx)),diff(interval2(1:nx))); 
        eegdata_L = eeg1(samples,:);
        eegdata_R = eeg2(samples,:);
        tempf_L = fft(eegdata_L,[],1);
        tempf_R = fft(eegdata_R,[],1);
        
        for m = 1:50 % 50 frequencies
            freqs = (m-1)*10+1:m*10; % downsample freq
            power_L = sum(abs(tempf_L(freqs,:)),1);
            power_R = sum(abs(tempf_R(freqs,:)),1);
            pow_L(j).EEG(m,:,k) = power_L; % freq x chan x time
            pow_R(j).EEG(m,:,k) = power_R;
        end
        
    end
    
end

cd /ssd/zhibin/1overf/20220610_2P/Segmented_data/Plots/EEG_power
condi=4;
for chan=1:32
    figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(2,1,1);
    imagesc(squeeze(pow_L(condi).EEG(:,chan,:)));colorbar;caxis([0 2.4e4]);
    xlabel('windows');ylabel('freq');
    title(['L-condi-' num2str(condi) '-chan-' labels{chan}]);
    subplot(2,1,2);
    imagesc(squeeze(pow_R(condi).EEG(:,chan,:)));colorbar;caxis([0 2.4e4]);
    title(['R-condi-' num2str(condi) '-chan-' labels{chan}]);
    xlabel('windows');ylabel('freq');

    figureName=['condi-' num2str(condi) '-chan-' labels{chan}];
    saveas(gcf,figureName,'jpg');
    close;
end
    
