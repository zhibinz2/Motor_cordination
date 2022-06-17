%%
clear, close all 

% load data
tempEEG = load('/ssd/zhibin/1overf/20220609_2P/Segmented_data/1_50Hz_ICAautomized/EEG20220609.mat'); 
tempBP = load('/ssd/zhibin/1overf/20220609_2P/Segmented_data/1_50Hz_ICAautomized/BP20220609.mat'); 

% tempEEG = load('/ssd/zhibin/1overf/20220610_2P/Segmented_data/1_50Hz_ICAautomized/EEG20220610.mat'); 
% tempBP = load('/ssd/zhibin/1overf/20220610_2P/Segmented_data/1_50Hz_ICAautomized/BP20220610.mat'); 

tempEEG = load('/ssd/zhibin/1overf/20220515_2P/Segmented_data/1_50Hz_ICAautomized/EEG20220515.mat'); 
tempBP = load('/ssd/zhibin/1overf/20220515_2P/Segmented_data/1_50Hz_ICAautomized/BP20220515.mat'); 

tempEEG = load('/ssd/zhibin/1overf/20220517_2P/Segmented_data/1_50Hz_ICAautomized/EEG20220517.mat'); 
tempBP = load('/ssd/zhibin/1overf/20220517_2P/Segmented_data/1_50Hz_ICAautomized/BP20220517.mat'); 

tempEEG = load('/ssd/zhibin/1overf/20220518_2P/Segmented_data/1_50Hz_ICAautomized/EEG20220518.mat'); 
tempBP = load('/ssd/zhibin/1overf/20220518_2P/Segmented_data/1_50Hz_ICAautomized/BP20220518.mat'); 

%% run PLS
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
% figure;imagesc(plsmodel(1).weights);colorbar;caxis([-2 2]*10E-7);
% figure;imagesc(plsmodel(2).weights);colorbar;caxis([-2 2]*10E-7);
% figure;imagesc(plsmodel(3).weights);colorbar;caxis([-2 2]*10E-7);
figure;imagesc(plsmodel(4).weights);colorbar;caxis([-2 2]*10E-7);

% finer labels
yticks([1:50]);
labels = {' 1-Fp1',' 2-Fpz',' 3-Fp2',' 4-F7',' 5-F3',' 6-Fz',' 7-F4',' 8-F8','9-FC5','10-FC1','11-FC2','12-FC6','13-M1','14-T7','15-C3','16-Cz','17-C4','18-T8','19-M2','20-CP5','21-CP1','22-CP2','23-CP6','24-P7','25-P3','26-Pz','27-P4','28-P8','29-POz','30-O1','31-Oz','32-O2'};
xticks([1:64]);
xticklabels([labels labels]);xtickangle(90);

title(['PLS 2']);

% figureName=['plsmodel'];
% saveas(gcf,figureName,'fig');
    
    
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

figureName=['plsmodel1'];
saveas(gcf,figureName,'fig');
%% all four conditions
figure('units','normalized','outerposition',[0 0 1 0.5]);
for i=1:4
    subplot(1,4,i);
        imagesc(plsmodel(i).weights);colorbar;caxis([-2 2]*10E-7);
        yticks([1:50]);
        xticks([1:64]);
        xticklabels([labels labels]);xtickangle(90);
        title({['Condi ' conditionNames{i}],'PLS model'},'Color',condicolors(i,:)); 
end
figureName=['4plsmodels'];
saveas(gcf,figureName,'fig');
    
%% load channel info for topoplots
cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/channels_info
load('chaninfo.mat')
cd /ssd/zhibin/1overf/20220609_2P/Segmented_data/Plots/SFN_plots

MatPlot=plsmodel(4).weights;

figure;
subplot(1,2,1);
topoplot(MatPlot(freq_inst,1:32),chaninfo,'nosedir','+X');colorbar;
title(['subject L']);caxis([-2 2]*10E-7);
subplot(1,2,2);
topoplot(MatPlot(freq_inst,33:64),chaninfo,'nosedir','+X');colorbar;
title(['subject R']);caxis([-2 2]*10E-7);
suptitle(['condition: ' num2str(4) '  frequency: ' num2str(freq_inst) ' Hz']);

figureName=['LR_topoplots'];
saveas(gcf,figureName,'fig');
    
%% all four conditions
figure('units','normalized','outerposition',[0 0 1 0.5]);
for i=1:4
    MatPlot=plsmodel(i).weights;
    subplot(1,8,2*i-1);
    topoplot(MatPlot(freq_inst,1:32),chaninfo,'nosedir','+X');
    colorbar;caxis([-2 2]*10E-7);
    title({['Condi ' conditionNames{i}],['subject L   Freq  ' num2str(freq_inst) 'Hz']},...
        'Color',condicolors(i,:));
    
    
    subplot(1,8,2*i);
    topoplot(MatPlot(freq_inst,33:64),chaninfo,'nosedir','+X');
    colorbar;caxis([-2 2]*10E-7);
    title({['Condi ' conditionNames{i}],['subject R   Freq  ' num2str(freq_inst) 'Hz']},...
        'Color',condicolors(i,:));
end
figureName=['4topoplots'];
saveas(gcf,figureName,'fig');


%% PLOT 14-4:  EEG Spectra - single chan - log log smooth at 10s  plot each condition with 3 levels
% [ Pspectra(EEG_L) + Xspectra(EEG_L vs EEG_R) + Pspectra(EEG_R) ]
% each of the 6 conditions (subplots(3,6,i)) for each of the 32 channels;
% auto save the 32 figures in a folder 
cd /ssd/zhibin/1overf/20220518_2P/Segmented_data/EEGspectra_chans_coh_PLOT_14_4
% select number of conditions in the experiment
Totalcondi=4;
% Totalcondi=5;
% Totalcondi=6;

maxfreq=50;sr=2000;win=10; % win in seconds (for cutting into epochs)
HML_labels={'High corr','Medium corr','Low corr'};

for chan = 25
    figure('units','normalized','outerposition',[0 0 1 1]);
    for condi=1:Totalcondi
        if condi==1; EEG_L=EEGCondi1L; EEG_R=EEGCondi1R; BP_L=BPCondi1L; BP_R=BPCondi1R; end 
        if condi==2; EEG_L=EEGCondi2L; EEG_R=EEGCondi2R; BP_L=BPCondi2L; BP_R=BPCondi2R; end 
        if condi==3; EEG_L=EEGCondi3L; EEG_R=EEGCondi3R; BP_L=BPCondi3L; BP_R=BPCondi3R; end
        if condi==4; EEG_L=EEGCondi4L; EEG_R=EEGCondi4R; BP_L=BPCondi4L; BP_R=BPCondi4R; end 
        if condi==5; EEG_L=EEGCondi5L; EEG_R=EEGCondi5R; BP_L=BPCondi5L; BP_R=BPCondi5R; end % R's Error base on L (3Hz bidirection)
        if condi==6; EEG_L=EEGCondi6L; EEG_R=EEGCondi6R; BP_L=BPCondi6L; BP_R=BPCondi6R; end 
        
        % calculate for each condition
        [~,~,~,~,~,~,xx1,xx2]  = spectra(EEG_L(:,chan),EEG_R(:,chan),maxfreq,sr,win);
        
        % get indices of the 3 levels
        [IndHigh,IndMedium,IndLow,corrBP] = TriSeg(BP_L,BP_R,sr,win);
        % loop through 3 levels
        for i=1:3 % 3 levels of correlation in BP
            % select the level of EEG
            if i==1; epochs=IndHigh;
            elseif i==2;epochs=IndMedium;
            else i==3; epochs=IndLow; end
            
            % calculate for each level
            [P1,P2,Px,coh,phase,freqs]  = spectra2(xx1(:,epochs),xx2(:,epochs),maxfreq,sr,win);
            
            subplot(4,4,condi); % P1
                hold on;
                loglog(smoothing(freqs,win),smoothing(P1,win),'Color',HNLcolors(i,:));
                if i==3;
                    legend(HML_labels(1:3));
                    xlabel('log-freqs');ylabel('log-Power');
                    ylim([0 0.025]);xlim([1.5 6]);
                    title({['Condi ' conditionNames{condi}],['Spectra EEG-L']},'Color',condicolors(condi,:));
                    hold off;
                end

            subplot(4,4,4+condi); % P2
                hold on;
                loglog(smoothing(freqs,win),smoothing(P2,win),'Color',HNLcolors(i,:));
                if i==3;
                    legend(HML_labels(1:3));
                    xlabel('log-freqs');ylabel('log-Power');
                    ylim([0 0.025]);xlim([1.5 6]);
                    title({['Spectra EEG-R']},'Color',condicolors(condi,:));
                    hold off;
                end

            subplot(4,4,4*2+condi); % Px
                hold on;
                loglog(smoothing(freqs,win),smoothing(abs(Px),win),'Color',HNLcolors(i,:));
                if i==3;
                    legend(HML_labels(1:3));
                    xlabel('log-freqs');ylabel('log-Power');
                    ylim([0 0.0045]);xlim([1.5 6]);
                    title({['XSpectra']},'Color',condicolors(condi,:));
                    hold off;
                end

            subplot(4,4,4*3+condi); % Coh
                hold on;
                % plot(freqs,coh,'Color',condicolors(condi,:));imagesc(coh);colorbar;
                plot(freqs,coh,'Color',HNLcolors(i,:));
                if i==3;
                    legend(HML_labels(1:3),'location','northeast');
                    xlabel('freqs');ylabel('coh');
                    ylim([0 0.8]);xlim([1.5 6]);
                    title({['Coh']},'Color',condicolors(condi,:));
                    hold off;
                end
        end

    end
    suptitle(['EEG-loglogspectra-coh-10s' labels{chan} ' subject ' num2str(seed) '-PLOT-14-4']);
    figureName=['EEG-loglogspectra-coh-10s-chan' labels{chan} '-PLOT-14-4'];
    % save the figure
    saveas(gcf,figureName,'fig');
end