%% try out spectra
clear
cd /ssd/zhibin/1overf/20220721_2P/Cleaned_data
load('clean_20220721.mat')
addpath /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/spectra


figure('units','normalized','outerposition',[0 0 1 1]);
tic
cmin=0;cmax=1e-4;
for s=1:2
    data=[];
    if s==1; data=dataL;else data=dataR;end
for t=1:12
    epoches=[];
    pow=[]; corr=[];ampcorr=[];fcoef=[];
    % cut data into epoches of 1 sec
    [epoches] = intoepoches(data{t}(:,1:32),sr,1);
    [pow,freqs,~,corr,ampcorr,fcoef] = allspectra1(epoches,sr,50);
    subplot(2,12,(s-1)*12+t);
    imagesc(pow);colorbar;colormap('jet');
    clim([cmin cmax]);
end
end
toc

%% compute pow in all sessions 
% original order of the experiment sessions
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1); 
% get the sequence of conditions in all sessions
run get_condition_all.m
% get the pow matrix in all essions
run get_pow_all.m 

%% compute average pow for each condition
% name of the 4 states
states4names={'Uncoupled','Leading','Following','Mutual'};
run avg_state_pow.m

%% plotting
% set plot colors
run plot_settings.m
% draw the plots
run draw_plots.m

%% compute coh and ampcorr
clear
cd /ssd/zhibin/1overf/20220721_2P/Cleaned_data
load('clean_20220721.mat','dataL')
addpath /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/spectra/pow_4states

dataL;
sr=2000;
for t=1
    epoches=[];
    pow=[]; coh=[];ampcorr=[];fcoef=[];
    % cut data into epoches of 1 sec
    [epoches] = intoepoches(dataL{t}(:,1:32),sr,1);
    [pow,~,~,coh,~,ampcorr,fcoef] = allspectra1(epoches,sr,50);
    % imagesc(pow);colorbar;colormap('jet');
    imagesc(squeeze(coh(6,:,:)));colorbar;colormap('jet');
    imagesc(squeeze(ampcorr(6,:,:)));colorbar;colormap('jet');
    % clim([cmin cmax]);
end

%% compute coh between L & R
clear
cd /ssd/zhibin/1overf/20220721_2P/Cleaned_data
load('clean_20220721.mat')
cd /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/spectra
% combine dataR and dataL 
for t=1:12
    data=
% cut data into epoches of 1 sec
[epoches] = intoepoches(dataL{t}(:,1:32),sr,1);

            




