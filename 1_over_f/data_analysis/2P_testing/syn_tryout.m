%%
clear, close all 

% load data
tempEEG = load('/ssd/zhibin/1overf/20220609_2P/Segmented_data/1_50Hz_ICAautomized/EEG20220609.mat'); 
tempBP = load('/ssd/zhibin/1overf/20220609_2P/Segmented_data/1_50Hz_ICAautomized/BP20220609.mat'); 

tempEEG = load('/ssd/zhibin/1overf/20220610_2P/Segmented_data/1_50Hz_ICAautomized/EEG20220610.mat'); 
tempBP = load('/ssd/zhibin/1overf/20220610_2P/Segmented_data/1_50Hz_ICAautomized/BP20220610.mat'); 

run sync.m

%% test
% first condition - 50freqs x 64chans x epoch 7
plot(pow(1).EEG(:,:,7)) 

plsmodel % all 4 conditions
plsmodel(1) % first condition

% look at the components same as the weights (1st component)
figure; imagesc(reshape(reg{1},50,64)); colorbar; title('reg -1');
figure; imagesc(reshape(reg{2},50,64)); colorbar; title('reg -2');
figure; imagesc(reshape(reg{3},50,64)); colorbar; title('reg -2');

%% plot
% look at the 1st component - the weight  in each condition
figure;imagesc(plsmodel(1).weights);colorbar;caxis([-2 2]*10E-7);
figure;imagesc(plsmodel(2).weights);colorbar;caxis([-2 2]*10E-7);
figure;imagesc(plsmodel(3).weights);colorbar;caxis([-2 2]*10E-7);
figure;imagesc(plsmodel(4).weights);colorbar;caxis([-2 2]*10E-7);

% finer labels
yticks([1:50]);

labels = {' 1-Fp1',' 2-Fpz',' 3-Fp2',' 4-F7',' 5-F3',' 6-Fz',' 7-F4',' 8-F8','9-FC5','10-FC1','11-FC2','12-FC6','13-M1','14-T7','15-C3','16-Cz','17-C4','18-T8','19-M2','20-CP5','21-CP1','22-CP2','23-CP6','24-P7','25-P3','26-Pz','27-P4','28-P8','29-POz','30-O1','31-Oz','32-O2'};

xticks([1:64]);
xticklabels([labels labels]);xtickangle(90);

title(['syncopation condition 4']);

% pick one point to plot that line
[x, y] = ginput(1); % read two mouse clicks on the plot % x were index, y were real values
% get the proximate index
string(y)
string(x) 
freq_inst=round(y)
chan_instR=round(x)-32
chan_instL=round(x)

chan_inst=25; % P3
chan_instL=27; % P4
chan_instR=13; % M1
%% load channel info for topoplots
cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/channels_info
load('chaninfo.mat')

MatPlot=plsmodel(4).weights;

figure;
subplot(1,2,1);
topoplot(MatPlot(freq_inst,1:32),chaninfo,'nosedir','+X');colorbar;
title(['subject L']);caxis([-2 2]*10E-7);
subplot(1,2,2);
topoplot(MatPlot(freq_inst,33:64),chaninfo,'nosedir','+X');colorbar;
title(['subject R']);caxis([-2 2]*10E-7);
suptitle(['condition: ' num2str(4) '  frequency: ' num2str(freq_inst) ' Hz']);