%% data localtion
cd /ssd/zhibin/1overf/20220515_2P/Segmented_data/Plots
cd /ssd/zhibin/1overf/20220517_2P/Segmented_data/Plots
cd /ssd/zhibin/1overf/20220518_2P/Segmented_data/Plots
cd /ssd/zhibin/1overf/20220609_2P/Segmented_data/Plots
cd /ssd/zhibin/1overf/20220610_2P/Segmented_data/Plots

%%
EEG_L(:,15);EEG_R(:,15);
minLength=min([length(EEG_L(:,15)) length(EEG_R(:,15))]);
minLength=2000;
eegdata=[EEG_L(1:minLength,15) EEG_R(1:minLength,15)]; 
figure;plot(eegdata);
samplingrate=2000;
lagstep=1;
maxlag=512;
badchans=[];

[covavg,goodchansk,lags] = getcovlag(eegdata(:,1:2),samplingrate,lagstep,maxlag,badchans);


covavgmean = zeros(2,2,length(lags));
covavgmean(goodchansk,goodchansk,:) = covavgmean(goodchansk,goodchansk,:)+covavg;

covavgmean(:,:,1) = (covavgmean(:,:,1))./minLength;

ncomponents=2;
[sobitest] = getsobi(covavgmean,ncomponents); 

eig=sobitest.eig;
corr=sobitest.corr;
pca=sobitest.pca;
pcaeig=pca.eig;
pcacorr=pca.corr;
pcaval=pca.val;

sobieeg = squeeze(eegdata)*eig;

figure;
plot(sobieeg);
