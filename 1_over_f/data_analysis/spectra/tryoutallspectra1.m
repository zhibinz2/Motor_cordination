clear
cd /ssd/zhibin/1overf/20220721_2P/Cleaned_data
load('clean_20220721.mat')
cd /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/spectra

figure;
for s=1:2
for t=1:12
    epoches=[];
    pow=[]; corr=[];ampcorr=[];fcoef=[];
    % cut data into epoches of 1 sec
    [epoches] = intoepoches(dataL{t}(:,1:32),sr,1);
    [pow,~,~,corr,ampcorr,fcoef] = allspectra1(epoches,sr,50);
    subplot(2,12,t);
    imagesc(pow)

