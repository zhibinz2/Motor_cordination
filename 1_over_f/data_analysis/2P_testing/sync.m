% clear, close all 
% open 8 processors 

samplingrate = 2000;
winsize = 10*samplingrate;

% load data

% tempEEG = load('pilot/20220609_2P/Segmented_data/1_50Hz_ICAautomized/EEG20220609.mat'); 
% tempBP = load('pilot/20220609_2P/Segmented_data/1_50Hz_ICAautomized/BP20220609.mat'); 

% tempEEG = load('/ssd/zhibin/1overf/20220609_2P/Segmented_data/1_50Hz_ICAautomized/EEG20220609.mat'); 
% tempBP = load('/ssd/zhibin/1overf/20220609_2P/Segmented_data/1_50Hz_ICAautomized/BP20220609.mat'); 

% tempEEG = load('/ssd/zhibin/1overf/20220610_2P/Segmented_data/1_50Hz_ICAautomized/EEG20220610.mat'); 
% tempBP = load('/ssd/zhibin/1overf/20220610_2P/Segmented_data/1_50Hz_ICAautomized/BP20220610.mat'); 

%im only doing 4 conditions here, and I'm gonna put them in an array structure to make me be able to loop over conditions.   

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
interval1= find(data(j,1).BP(samples)); 
interval2= find(data(j,2).BP(samples));
nx = min([length(interval1) length(interval2)]);
cor(j).BP(k) = corr(diff(interval1(1:nx)),diff(interval2(1:nx))); 
eegdata = [eeg1(samples,:) eeg2(samples,:)];
tempf = fft(eegdata,[],1);
for m = 1:50 % 50 frequencies
freqs = (m-1)*10+1:m*10; % downsample freq
power = sum(abs(tempf(freqs,:)),1);
pow(j).EEG(m,:,k) = power;
end
end 
end

for j = 1:4 % four conditions
powforpls = reshape(pow(j).EEG,50*64,nwin(j))';   
%I am doing a 1 factor prediction model.  could test more factors. 
[R2,reg,ypred] = npls_pred(powforpls,cor(j).BP',1);
plsmodel(j).weights = reshape(reg{1},50,64);
plsmodel(j).R2 = R2;
plsmodel(j).ypred = ypred;
end;

return 

