%% Plots - syncopation - color scheme
red   = [1 0 0];
pink  = [1 0.65 0.75];
blue  = [0 0 1];
green = [0 1 0];
darkgreen = [0 0.5 0];
grey  = [0.5 0.5 0.5];
yellow  = [1 1 0];
deepyellow  = [1 0.8 0.2];
megenta = [1 0 1];% fill([0 1 1 0],[0 0 1 1],megenta)
cyan = [0 1 1]; % fill([0 1 1 0],[0 0 1 1],cc)
purple = [0.6 0.1 0.9];
condicolors=[darkgreen;red;blue;megenta;purple;purple];
HNLcolors=[darkgreen; deepyellow; pink]
% % test color
% showcolor=pink;
% imagesc(cat(3,showcolor(1),showcolor(2),showcolor(3)));

%% conditionNames
conditionNames; % check if correct

% syncopation expt
conditionNames={'uncoupled' 'L-lead' 'R-lead' 'mutual-1.3Hz'};


%% EEG chan labels
labels = {' 1-Fp1',' 2-Fpz',' 3-Fp2',' 4-F7',' 5-F3',' 6-Fz',' 7-F4',' 8-F8','9-FC5','10-FC1','11-FC2','12-FC6','13-M1','14-T7','15-C3','16-Cz','17-C4','18-T8','19-M2','20-CP5','21-CP1','22-CP2','23-CP6','24-P7','25-P3','26-Pz','27-P4','28-P8','29-POz','30-O1','31-Oz','32-O2'};

%% PLOT 1:  hist BP intervals  
% [ hist(BP_L) + hist(FB_L) + hist(FP_R) + hist(BP_R)  ]
% 6 conditions -subplots(4,4,condi) 
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing

figure('units','normalized','outerposition',[0 0 1 1]);
for condi=1:4 % for condi=1:5
    if condi==1; BP_L=BPCondi1L; FB_L=FBCondi1L; BP_R=BPCondi1R;FB_R=FBCondi1R;end
    if condi==2; BP_L=BPCondi2L; FB_L=FBCondi2L; BP_R=BPCondi2R;FB_R=FBCondi2R;end
    if condi==3; BP_L=BPCondi3L; FB_L=FBCondi3L; BP_R=BPCondi3R;FB_R=FBCondi3R;end
    if condi==4; BP_L=BPCondi4L; FB_L=FBCondi4L; BP_R=BPCondi4R;FB_R=FBCondi4R;end
    % compute
    subplot(4,4,condi); %  hist(BP_L)
        histogram(Calinterval(BP_L')./2,600,'FaceColor',condicolors(condi,:),'EdgeColor',condicolors(condi,:));
        xlabel('BP-intervals (ms)');ylabel('number of BP intervals');xlim([0 1200]);
        legend(['total ' num2str(length(Calinterval(BP_L')))]);
        title({['Condi ' conditionNames{condi}],' BP-R interval'},'Color',condicolors(condi,:));
    subplot(4,4,4+condi); %  hist(FB_L)
        histogram(Calinterval(FB_L')./2,600,'FaceColor',condicolors(condi,:),'EdgeColor',condicolors(condi,:));
        xlabel('FB-intervals (ms)');ylabel('number of FB intervals');xlim([0 1200]);
        legend(['total ' num2str(length(Calinterval(FB_L')))]);
        title('FB-L','Color',condicolors(condi,:));
    subplot(4,4,4*2+condi); %   hist(FP_R)
        histogram(Calinterval(FB_R')./2,600,'FaceColor',condicolors(condi,:),'EdgeColor',condicolors(condi,:));
        xlabel('FB-intervals (ms)');ylabel('number of FB intervals');xlim([0 1200]);
        legend(['total ' num2str(length(Calinterval(FB_R')))]);
        title('FB-R','Color',condicolors(condi,:));
    subplot(4,4,4*3+condi); %  hist(BP_R)
        histogram(Calinterval(BP_R')./2,600,'FaceColor',condicolors(condi,:),'EdgeColor',condicolors(condi,:));
        xlabel('BP-intervals (ms)');ylabel('number of BP intervals');xlim([0 1200]);
        legend(['total ' num2str(length(Calinterval(BP_R')))]);
        title('BP-R','Color',condicolors(condi,:));
end
suptitle(['Button press intervals ' ' -- subject ' num2str(seed) '-PLOT-1']);

% cd /ssd/zhibin/1overf/20220518_2P/Segmented_data/Plots
figureName=['BP-hist -- ' num2str(seed) '-PLOT-1'];
% save the figure
saveas(gcf,figureName,'fig');

%% PLOT 1-1:  BP-FB intervals trend
% [(BP_L + FB_L) + (FP_R + BP_R)]+smooth[(BP_L + FB_L) + (FP_R + BP_R)]
% 6 conditions -subplots(2,6,condi) 
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing

figure('units','normalized','outerposition',[0 0 1 1]);
win=20;
for condi=1:4
    if condi==1; BP_L=BPCondi1L; FB_L=FBCondi1L; BP_R=BPCondi1R; FB_R=FBCondi1R;end
    if condi==2; BP_L=BPCondi2L; FB_L=FBCondi2L; BP_R=BPCondi2R; FB_R=FBCondi2R;end
    if condi==3; BP_L=BPCondi3L; FB_L=FBCondi3L; BP_R=BPCondi3R; FB_R=FBCondi3R;end
    if condi==4; BP_L=BPCondi4L; FB_L=FBCondi4L; BP_R=BPCondi4R; FB_R=FBCondi4R;end
    % compute
    subplot(4,4,condi); %  (BP_L + FB_L)
        BP_inter1=Calinterval(BP_L')./2; BP_inter2=Calinterval(FB_L')./2;
        BPfrac_time1=linspace(1,5,length(BP_inter1)); BPfrac_time2=linspace(1,5,length(BP_inter2)); 
        plot(BPfrac_time1,BP_inter1,'r');hold on; plot(BPfrac_time2,BP_inter2,'b');
        xlabel('fractional time (min)');ylabel('interval (ms)');
        legend('BP-L','FB-L');
        title({['Condi ' conditionNames{condi}],'L recording'},'Color',condicolors(condi,:));
        ylim([500 1000]);
    subplot(4,4,4+condi); %  (FP_R + BP_R)
        BP_inter1=Calinterval(BP_R')./2; BP_inter2=Calinterval(FB_R')./2;
        BPfrac_time1=linspace(1,5,length(BP_inter1)); BPfrac_time2=linspace(1,5,length(BP_inter2)); 
        plot(BPfrac_time1,BP_inter1,'b');hold on; plot(BPfrac_time2,BP_inter2,'r');
        xlabel('fractional time (min)');ylabel('interval (ms)');
        legend('BP-R','FB-R');
        title('R recording','Color',condicolors(condi,:));
        ylim([500 1000]);
    subplot(4,4,4*2+condi); %   smooth[(BP_L + FB_L)]
        BP_inter1=smoothing(Calinterval(BP_L')./2,win); BP_inter2=smoothing(Calinterval(FB_L')./2,win);
        BPfrac_time1=linspace(1,5,length(BP_inter1)); BPfrac_time2=linspace(1,5,length(BP_inter2)); 
        plot(BPfrac_time1,BP_inter1,'r');hold on; plot(BPfrac_time2,BP_inter2,'b');
        xlabel('fractional time (min)');ylabel('interval (ms)');
        legend('BP-L','FB-L');
        title('L recording smoothed','Color',condicolors(condi,:));
        ylim([500 1000]);
    subplot(4,4,4*3+condi); % smooth[(FP_R + BP_R)]
        BP_inter1=smoothing(Calinterval(BP_R')./2,win); BP_inter2=smoothing(Calinterval(FB_R')./2,win);
        BPfrac_time1=linspace(1,5,length(BP_inter1)); BPfrac_time2=linspace(1,5,length(BP_inter2)); 
        plot(BPfrac_time1,BP_inter1,'b');hold on; plot(BPfrac_time2,BP_inter2,'r');
        xlabel('fractional time (min)');ylabel('interval (ms)');
        legend('BP-R','FB-R');
        title('R recording smoothed','Color',condicolors(condi,:));
        ylim([500 1000]);
end
suptitle(['Fractional button press intervals smoothed win ' num2str(win) ' -- subject ' num2str(seed) '-PLOT-1-1']);

% cd /ssd/zhibin/1overf/20220518_2P/Segmented_data/Plots
figureName=['BP-FB-frac-smooth-win-' num2str(win) '----' num2str(seed) '-PLOT-1-1'];
% save the figure
saveas(gcf,figureName,'jpg');

%% PLOT 1-2:  BP intervals smoothed trend
% [(BP_L + FB_L) + (FP_R + BP_R)]+smooth[(BP_L + FB_L) + (FP_R + BP_R)]
% 6 conditions -subplots(2,6,condi) 
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing

figure('units','normalized','outerposition',[0 0 1 1]);
for condi=1:4
    if condi==1; BP_L=BPCondi1L; FB_L=FBCondi1L; BP_R=BPCondi1R; FB_R=FBCondi1R; end
    if condi==2; BP_L=BPCondi2L; FB_L=FBCondi2L; BP_R=BPCondi2R; FB_R=FBCondi2R; end
    if condi==3; BP_L=BPCondi3L; FB_L=FBCondi3L; BP_R=BPCondi3R; FB_R=FBCondi3R; end
    if condi==4; BP_L=BPCondi4L; FB_L=FBCondi4L; BP_R=BPCondi4R; FB_R=FBCondi4R; end
    % compute
    subplot(4,4,condi); 
        win=1;
        BP_inter1=smoothing(Calinterval(BP_L')./2,win); BP_inter2=smoothing(Calinterval(BP_R')./2,win);
        BPfrac_time1=linspace(1,5,length(BP_inter1)); BPfrac_time2=linspace(1,5,length(BP_inter2)); 
        plot(BPfrac_time1,BP_inter1,'r');hold on; plot(BPfrac_time2,BP_inter2,'b');
        xlabel('fractional time (min)');ylabel('interval (ms)');
        legend('BP-L','BP-R');
        title({['Condi ' conditionNames{condi}],['BP smooth win ' num2str(win)]},'Color',condicolors(condi,:));
        ylim([400 1000]);
    subplot(4,4,4+condi); 
        win=20;
        BP_inter1=smoothing(Calinterval(BP_L')./2,win); BP_inter2=smoothing(Calinterval(BP_R')./2,win);
        BPfrac_time1=linspace(1,5,length(BP_inter1)); BPfrac_time2=linspace(1,5,length(BP_inter2)); 
        plot(BPfrac_time1,BP_inter1,'r');hold on; plot(BPfrac_time2,BP_inter2,'b');
        xlabel('fractional time (min)');ylabel('interval (ms)');
        legend('BP-L','BP-R');
        title(['BP smooth win ' num2str(win)],'Color',condicolors(condi,:));
        ylim([400 1000]);
    subplot(4,4,4*2+condi); 
        win=50;
        BP_inter1=smoothing(Calinterval(BP_L')./2,win); BP_inter2=smoothing(Calinterval(BP_R')./2,win);
        BPfrac_time1=linspace(1,5,length(BP_inter1)); BPfrac_time2=linspace(1,5,length(BP_inter2)); 
        plot(BPfrac_time1,BP_inter1,'r');hold on; plot(BPfrac_time2,BP_inter2,'b');
        xlabel('fractional time (min)');ylabel('interval (ms)');
        legend('BP-L','BP-R');
        title(['BP smooth win ' num2str(win)],'Color',condicolors(condi,:));
        ylim([400 1000]);
    subplot(4,4,4*3+condi); 
        win=100;
        BP_inter1=smoothing(Calinterval(BP_L')./2,win); BP_inter2=smoothing(Calinterval(BP_R')./2,win);
        BPfrac_time1=linspace(1,5,length(BP_inter1)); BPfrac_time2=linspace(1,5,length(BP_inter2)); 
        plot(BPfrac_time1,BP_inter1,'r');hold on; plot(BPfrac_time2,BP_inter2,'b');
        xlabel('fractional time (min)');ylabel('interval (ms)');
        legend('BP-L','BP-R');
        title(['BP smooth win ' num2str(win)],'Color',condicolors(condi,:));
        ylim([400 1000]);
end
suptitle(['Fractional BP intervals smoothed -- subject ' num2str(seed) '-PLOT-1-2']);

% cd /ssd/zhibin/1overf/20220518_2P/Segmented_data/Plots
figureName=['BP-frac-smooth----' num2str(seed) '-PLOT-1-2'];
% save the figure
saveas(gcf,figureName,'jpg');

%% PLOT 1-3:  BP intervals smoothed trend corr
% [(BP_L + FB_L) + (FP_R + BP_R)]+smooth[(BP_L + FB_L) + (FP_R + BP_R)]
% 6 conditions -subplots(2,6,condi) 
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing
% select number of conditions in the experiment
Totalcondi=4;
Totalcondi=5;
Totalcondi=6;

figure('units','normalized','outerposition',[0 0 1 1]);
for condi=1:Totalcondi
    if condi==1; BP_L=BPCondi1L; FB_L=FBCondi1L; BP_R=BPCondi1R; FB_R=FBCondi1R; end
    if condi==2; BP_L=BPCondi2L; FB_L=FBCondi2L; BP_R=BPCondi2R; FB_R=FBCondi2R; end
    if condi==3; BP_L=BPCondi3L; FB_L=FBCondi3L; BP_R=BPCondi3R; FB_R=FBCondi3R; end
    if condi==4; BP_L=BPCondi4L; FB_L=FBCondi4L; BP_R=BPCondi4R; FB_R=FBCondi4R; end
    if condi==5; BP_L=BPCondi5L; FB_L=FBCondi5L; BP_R=BPCondi5R;FB_R=FBCondi5R; end
    if condi==6; BP_L=BPCondi6L; FB_L=FBCondi6L; BP_R=BPCondi6R;FB_R=FBCondi6R; end
    % compute
    subplot(2,Totalcondi,condi); 
        win=13;
        BP_inter1=smoothing(Calinterval(BP_L')./2,win); BP_inter2=smoothing(Calinterval(BP_R')./2,win);
        Maxlength=max([length(BP_inter1) length(BP_inter2)]);
        resam_BP_inter1=resample(BP_inter1',Maxlength,length(BP_inter1));
        resam_BP_inter2=resample(BP_inter2',Maxlength,length(BP_inter2));
        % yyaxis left;ylabel('interval (ms)','color','m'); set(gca,'ycolor','m');ylim([500 1000])
        plot(resam_BP_inter1,'r'); hold on; plot(resam_BP_inter2,'b'); hold off;
        ylabel('interval (ms)','color','m');% xlabel('taps (resampled)');
        title({['Condi ' conditionNames{condi}],['BP smooth win ' num2str(win)]},'Color',condicolors(condi,:));
        % ylim([400 850]);
        legend({'BP-L','BP-R'});
    subplot(2,Totalcondi,Totalcondi+condi); 
        % calculate correlation coeffiences win the same smoothing window
        corrSeries=[];
        for i=1:(Maxlength-win+1)
            corrSeries(i)=corr(resam_BP_inter1(i:i+win-1),resam_BP_inter2(i:i+win-1));
        end
        % plot on top in green
        % yyaxis right; ylabel('corr coef','color',darkgreen);set(gca,'ycolor',darkgreen); ylim([-1 1]);
        plot(smoothing(corrSeries,win),'k'); ylabel('corr coef'); xlabel('taps (resampled)');
        title(['corr coef smooth win ' num2str(win)],'Color',condicolors(condi,:));
        ylim([-1 1]);
        yline(0,'color',deepyellow);
end
suptitle(['Fractional BP intervals smoothed win' num2str(win) 'corr -- subject' num2str(seed) '-PLOT-1-3']);

% cd /ssd/zhibin/1overf/20220518_2P/Segmented_data/Plots
figureName=['BP-frac-smooth-corr----' num2str(seed) '-PLOT-1-3'];
% save the figure
saveas(gcf,figureName,'fig');

%% PLOT 1-3-1:  BP intervals smoothed trend corr
% [(BP_L + FB_L) + (FP_R + BP_R)]+smooth[(BP_L + FB_L) + (FP_R + BP_R)]
% 6 conditions -subplots(2,6,condi) 
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing
% select number of conditions in the experiment
Totalcondi=4;
Totalcondi=5;
Totalcondi=6;

figure('units','normalized','outerposition',[0 0 1 0.8]);
for condi=1:Totalcondi
    if condi==1; BP_L=BPCondi1L; FB_L=FBCondi1L; BP_R=BPCondi1R; FB_R=FBCondi1R; end
    if condi==2; BP_L=BPCondi2L; FB_L=FBCondi2L; BP_R=BPCondi2R; FB_R=FBCondi2R; end
    if condi==3; BP_L=BPCondi3L; FB_L=FBCondi3L; BP_R=BPCondi3R; FB_R=FBCondi3R; end
    if condi==4; BP_L=BPCondi4L; FB_L=FBCondi4L; BP_R=BPCondi4R; FB_R=FBCondi4R; end
    if condi==5; BP_L=BPCondi5L; FB_L=FBCondi5L; BP_R=BPCondi5R;FB_R=FBCondi5R;end
    if condi==6; BP_L=BPCondi6L; FB_L=FBCondi6L; BP_R=BPCondi6R;FB_R=FBCondi6R;end
    % compute
    subplot(2,Totalcondi,condi); 
        win=13;
        BP_inter1=smoothing(Calinterval(BP_L')./2,win); BP_inter2=smoothing(Calinterval(BP_R')./2,win);
        Maxlength=max([length(BP_inter1) length(BP_inter2)]);
        resam_BP_inter1=resample(BP_inter1',Maxlength,length(BP_inter1));
        resam_BP_inter2=resample(BP_inter2',Maxlength,length(BP_inter2));
        % yyaxis left;ylabel('interval (ms)','color','m'); set(gca,'ycolor','m');ylim([500 1000])
        lh=plot(resam_BP_inter1,'r'); hold on; plot(resam_BP_inter2,'b'); hold off;
        ylabel('interval (ms)','color','m');% xlabel('taps (resampled)');
        title({['Condi ' conditionNames{condi}],['BP smooth win ' num2str(win)]},'Color',condicolors(condi,:));
        ylim([200 800]);% legend({'BP-L','BP-R'});
        hold on;
        for i=[5 10 15 20 25 30 35 40 45]
            xline(13*i,'Color',condicolors(condi,:));
        end
        hold off;
    subplot(2,Totalcondi,Totalcondi+condi); 
        % calculate correlation coeffiences win the same smoothing window
        corrSeries=[];
        for i=1:(Maxlength-win+1)
            corrSeries(i)=corr(resam_BP_inter1(i:i+win-1),resam_BP_inter2(i:i+win-1));
        end
        % plot on top in green
        % yyaxis right; ylabel('corr coef','color',darkgreen);set(gca,'ycolor',darkgreen); ylim([-1 1]);
        SmoocorrSeries=smoothing(corrSeries,win);
        posInd=find(SmoocorrSeries>0);negInd=find(SmoocorrSeries<=0);
        plot(posInd,SmoocorrSeries(posInd),'g.'); hold on;
        plot(negInd,SmoocorrSeries(negInd),'k.');
        ylabel('corr coef'); xlabel('taps (resampled)');
        title(['corr coef smooth win ' num2str(win)],'Color',condicolors(condi,:));
        ylim([-1 1]);
        for i=[5 10 15 20 25 30 35 40 45]
            xline(13*i,'Color',condicolors(condi,:));
        end
        hold off;
end
suptitle(['Fractional BP intervals smoothed win' num2str(win) 'corr -- subject' num2str(seed) '-PLOT-1-3-1']);

% cd /ssd/zhibin/1overf/20220518_2P/Segmented_data/Plots
figureName=['BP-frac-smooth-corr----' num2str(seed) '-PLOT-1-3-1'];
% save the figure
saveas(gcf,figureName,'fig');

%% PLOT 1-4:  BP intervals smoothed trend boolean corr
% [(BP_L + FB_L) + (FP_R + BP_R)]+smooth[(BP_L + FB_L) + (FP_R + BP_R)]
% 6 conditions -subplots(2,6,condi) 
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing
% select number of conditions in the experiment
Totalcondi=4;
Totalcondi=5;
Totalcondi=6;

figure('units','normalized','outerposition',[0 0 1 0.8]);
for condi=1:Totalcondi
    if condi==1; BP_L=BPCondi1L; FB_L=FBCondi1L; BP_R=BPCondi1R; FB_R=FBCondi1R; end
    if condi==2; BP_L=BPCondi2L; FB_L=FBCondi2L; BP_R=BPCondi2R; FB_R=FBCondi2R; end
    if condi==3; BP_L=BPCondi3L; FB_L=FBCondi3L; BP_R=BPCondi3R; FB_R=FBCondi3R; end
    if condi==4; BP_L=BPCondi4L; FB_L=FBCondi4L; BP_R=BPCondi4R; FB_R=FBCondi4R; end
    if condi==5; BP_L=BPCondi5L; FB_L=FBCondi5L; BP_R=BPCondi5R;FB_R=FBCondi5R;end
    if condi==6; BP_L=BPCondi6L; FB_L=FBCondi6L; BP_R=BPCondi6R;FB_R=FBCondi6R;end
    % compute
    subplot(2,Totalcondi,condi); 
        win=13;
        BP_inter1=smoothing(Calinterval(BP_L')./2,win); BP_inter2=smoothing(Calinterval(BP_R')./2,win);
        Maxlength=max([length(BP_inter1) length(BP_inter2)]);
        resam_BP_inter1=resample(BP_inter1',Maxlength,length(BP_inter1));
        resam_BP_inter2=resample(BP_inter2',Maxlength,length(BP_inter2));
        % yyaxis left;ylabel('interval (ms)','color','m'); set(gca,'ycolor','m');ylim([500 1000])
        plot(resam_BP_inter1,'r'); hold on; plot(resam_BP_inter2,'b'); hold off;
        ylabel('interval (ms)','color','m');% xlabel('taps (resampled)');
        title({['Condi ' conditionNames{condi}],['BP smooth win ' num2str(win)]},'Color',condicolors(condi,:));
        % ylim([250 800]);
        legend({'BP-L','BP-R'});
    subplot(2,Totalcondi,Totalcondi+condi); 
        % calculate correlation coeffiences win the same smoothing window
        wins=10;
        winsample=sr*wins;
        Minwin=min([floor(length(BP_L)/winsample) floor(length(BP_R)/winsample)]);
        corrSeries=[];
        for i=1:Minwin
            corrSeries(i)=corr(BP_L((winsample*(i-1)+1):(winsample*i),:),BP_R((winsample*(i-1)+1):(winsample*i),:));
        end
        % plot on top in green
        % yyaxis right; ylabel('corr coef','color',darkgreen);set(gca,'ycolor',darkgreen); ylim([-1 1]);
        plot(corrSeries,'k'); ylabel('corr coef'); xlabel('taps (resampled)');
        title(['corr coef win ' num2str(wins) ' s'],'Color',condicolors(condi,:));
end
suptitle(['Fractional BP smoothed win' num2str(win) 's-corr01 -- subject' num2str(seed) '-PLOT-1-4']);

% cd /ssd/zhibin/1overf/20220518_2P/Segmented_data/Plots
figureName=['BP-frac-smooth-corr01----' num2str(seed) '-PLOT-1-4'];
% save the figure
saveas(gcf,figureName,'fig');

%% PLOT 2:  BP interval corr on interval
% [ autocorr(BP_L) + xcorr(BP_L vs FB_L) + xcorr(BP_R vs FB_R) + autocorr(BP_R) ]
% 6 conditions -subplots(4,6,condi) 
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing

figure('units','normalized','outerposition',[0 0 1 1]);
for condi=1:4
    if condi==1; BP_L=BPCondi1L; FB_L=FBCondi1L; BP_R=BPCondi1R; FB_R=FBCondi1R; end
    if condi==2; BP_L=BPCondi2L; FB_L=FBCondi2L; BP_R=BPCondi2R; FB_R=FBCondi2R; end
    if condi==3; BP_L=BPCondi3L; FB_L=FBCondi3L; BP_R=BPCondi3R; FB_R=FBCondi3R; end
    if condi==4; BP_L=BPCondi4L; FB_L=FBCondi4L; BP_R=BPCondi4R; FB_R=FBCondi4R; end
    % compute
    subplot(4,4,condi); %  autocorr(BP_L)
        autocorr(Calinterval(BP_L'),430); xlabel('lags');ylabel('autocorr'); 
        title({['Condi ' conditionNames{condi}],'autocorr BP-L'},'Color',condicolors(condi,:));
    subplot(4,4,4+condi); %  xcorr(BP_L vs FB_L)
        time_series1=Calinterval(BP_L');% plot(time_series1);
        time_series2=Calinterval(FB_L');% plot(time_series2);
        TimeLength=min([length(time_series1) length(time_series2)]);
        time_series1=time_series1-mean(time_series1);% plot(time_series1);
        time_series2=time_series2-mean(time_series2);% plot(time_series2);
        [r,lags]=xcorr(time_series1(1:TimeLength)', time_series2(1:TimeLength)', 15,'coeff');
        plot(lags,r);xlabel('lags');ylabel('Xcorr');xline(0,'--','Color',deepyellow);
        title('xcorr BP-L & FB-L','Color',condicolors(condi,:));
        ylim([-0.2 0.5]);
    subplot(4,4,4*2+condi); %   xcorr(BP_R vs FB_R) 
        time_series1=Calinterval(BP_R');% plot(time_series1);
        time_series2=Calinterval(FB_R');% plot(time_series2);
        TimeLength=min([length(time_series1) length(time_series2)]);
        time_series1=time_series1-mean(time_series1);% plot(time_series1);
        time_series2=time_series2-mean(time_series2);% plot(time_series2);
        [r,lags]=xcorr(time_series1(1:TimeLength), time_series2(1:TimeLength), 15,'normalized');
        plot(lags,r);xlabel('lags');ylabel('Xcorr');xline(0,'--','Color',deepyellow);
        title('xcorr BP-R & FB-R','Color',condicolors(condi,:));
        ylim([-0.2 0.5]);
    subplot(4,4,4*3+condi); %  autocorr(BP_R)
        autocorr(Calinterval(BP_R'),430); xlabel('lags');ylabel('autocorr'); 
        title('autocorr BP-R','Color',condicolors(condi,:));
end
suptitle(['PLOT 2: auto/x corr tapping intervals' ' -- subject ' num2str(seed) '-PLOT-2']);

% cd /ssd/zhibin/1overf/20220518_2P/Segmented_data/Plots    
figureName=['BPcorr--' num2str(seed) '-PLOT-2'];
    % save the figure
    saveas(gcf,figureName,'jpg');
    saveas(gcf,figureName,'fig');

%% PLOT 3:  BP interval corr smooth (better!) (short-range / local statistics; weak anticipation) on interval
% [ autocorr(BP_L) + xcorr(BP_L vs FB_L) + xcorr(BP_R vs FB_R) + autocorr(BP_R) ]
% 6 conditions -subplots(4,6,condi) 
% select number of conditions in the experiment
Totalcondi=4;
Totalcondi=5;
Totalcondi=6;

figure('units','normalized','outerposition',[0 0 1 1]);
win=13; % smooth with a sliding window size
for condi=1:Totalcondi
    if condi==1; BP_L=BPCondi1L; FB_L=FBCondi1L; BP_R=BPCondi1R;FB_R=FBCondi1R;end
    if condi==2; BP_L=BPCondi2L; FB_L=FBCondi2L; BP_R=BPCondi2R;FB_R=FBCondi2R;end
    if condi==3; BP_L=BPCondi3L; FB_L=FBCondi3L; BP_R=BPCondi3R;FB_R=FBCondi3R;end
    if condi==4; BP_L=BPCondi4L; FB_L=FBCondi4L; BP_R=BPCondi4R;FB_R=FBCondi4R;end 
    if condi==5; BP_L=BPCondi5L; FB_L=FBCondi5L; BP_R=BPCondi5R;FB_R=FBCondi5R;end
    if condi==6; BP_L=BPCondi6L; FB_L=FBCondi6L; BP_R=BPCondi6R;FB_R=FBCondi6R;end
    % compute
    subplot(4,Totalcondi,condi); %  autocorr(BP_L)
        autocorr(smoothing(Calinterval(BP_L'),win),430-win); xlabel('lags');ylabel('autocorr'); 
        title({['Condi ' conditionNames{condi}], 'autocorr BP-L'},'Color',condicolors(condi,:));
    subplot(4,Totalcondi,Totalcondi+condi); %  xcorr(BP_L vs FB_L)
        time_series1=smoothing(Calinterval(BP_L'),win);% plot(time_series1);
        time_series2=smoothing(Calinterval(FB_L'),win);% plot(time_series2);
        TimeLength=min([length(time_series1) length(time_series2)]);
        time_series1=time_series1-mean(time_series1);% plot(time_series1);
        time_series2=time_series2-mean(time_series2);% plot(time_series2);
%         [r,lags]=xcorr(time_series1(1:TimeLength), time_series2(1:TimeLength), 430-win,'normalized');
        [r,lags]=xcorr(time_series1(1:TimeLength), time_series2(1:TimeLength), 50,'normalized');
        plot(lags,r);xlabel('lags');ylabel('Xcorr');xline(0,'--','Color',deepyellow);
        title('xcorr BP-L & FB-L','Color',condicolors(condi,:));
        ylim([-0.4 1]);
    subplot(4,Totalcondi,Totalcondi*2+condi); %   xcorr(BP_R vs FB_R) 
        time_series1=smoothing(Calinterval(BP_R'),win);% plot(time_series1);
        time_series2=smoothing(Calinterval(FB_R'),win);% plot(time_series2);
        TimeLength=min([length(time_series1) length(time_series2)]);
        time_series1=time_series1-mean(time_series1);% plot(time_series1);
        time_series2=time_series2-mean(time_series2);% plot(time_series2);
%         [r,lags]=xcorr(time_series1(1:TimeLength), time_series2(1:TimeLength), 430-win,'normalized');
        [r,lags]=xcorr(time_series1(1:TimeLength), time_series2(1:TimeLength), 50,'normalized');
        plot(lags,r);xlabel('lags');ylabel('Xcorr');xline(0,'--','Color',deepyellow);
        title('xcorr BP-R & FB-R','Color',condicolors(condi,:));
        ylim([-0.4 1]);
    subplot(4,Totalcondi,Totalcondi*3+condi); %  autocorr(BP_R)
        autocorr(smoothing(Calinterval(BP_R'),win),430-win); xlabel('lags');ylabel('autocorr'); 
        title('autocorr BP-R','Color',condicolors(condi,:));
end
suptitle(['BP interval corr smooth win ' num2str(win) ' subject ' num2str(seed) '-PLOT-3']);
        
figureName=['BPcorr_smooth--' num2str(seed) '-PLOT-3'];
    % save the figure
    saveas(gcf,figureName,'fig');

%% PLOT 4:  BP corr boolean (low correlation) (short-range / local statistics; weak anticipation) on interval
% [ autocorr(BP_L) + xcorr(BP_L vs FB_L) + xcorr(BP_R vs FB_R) + autocorr(BP_R) ]
% 6 conditions -subplots(4,6,condi) 

figure('units','normalized','outerposition',[0 0 1 1]);
for condi=1:4
    if condi==1; BP_L=BPCondi1L; FB_L=FBCondi1L; BP_R=BPCondi1R;FB_R=FBCondi1R;end
    if condi==2; BP_L=BPCondi2L; FB_L=FBCondi2L; BP_R=BPCondi2R;FB_R=FBCondi2R;end
    if condi==3; BP_L=BPCondi3L; FB_L=FBCondi3L; BP_R=BPCondi3R;FB_R=FBCondi3R;end
    if condi==4; BP_L=BPCondi4L; FB_L=FBCondi4L; BP_R=BPCondi4R;FB_R=FBCondi4R;end
    % compute
    subplot(4,4,condi); %  autocorr(BP_L)
        autocorr(BP_L',1500); xlabel('lags');ylabel('autocorr'); ylim([0 0.03]);
        title({['Condi ' conditionNames{condi}], 'autocorr BP-L'},'Color',condicolors(condi,:));
    subplot(4,4,4+condi); %  xcorr(BP_L vs FB_L)
        time_series1=BP_L';% plot(time_series1);
        time_series2=FB_L';% plot(time_series2);
        TimeLength=min([length(time_series1) length(time_series2)]);
        time_series1=time_series1-mean(time_series1);% plot(time_series1);
        time_series2=time_series2-mean(time_series2);% plot(time_series2);
        [r,lags]=xcorr(time_series1(1:TimeLength), time_series2(1:TimeLength), 1500,'normalized');
        plot(lags./2,r);xlabel('time [ms]');ylabel('Xcorr');xline(0,'--','Color',deepyellow);
        title('xcorr BP-L & FB-L','Color',condicolors(condi,:));
        ylim([-5e-3 15e-3]);
    subplot(4,4,4*2+condi); %   xcorr(BP_R vs FB_R) 
        time_series1=BP_R';% plot(time_series1);
        time_series2=FB_R';% plot(time_series2);
        TimeLength=min([length(time_series1) length(time_series2)]);
        time_series1=time_series1-mean(time_series1);% plot(time_series1);
        time_series2=time_series2-mean(time_series2);% plot(time_series2);
        [r,lags]=xcorr(time_series1(1:TimeLength), time_series2(1:TimeLength), 1500,'normalized');
        plot(lags./2,r);xlabel('time [ms]');ylabel('Xcorr');xline(0,'--','Color',deepyellow);
        title('xcorr BP-R & FB-R','Color',condicolors(condi,:));
        ylim([-5e-3 15e-3]);
    subplot(4,4,4*3+condi); %  autocorr(BP_R)
        autocorr(BP_R',1500); xlabel('lags');ylabel('autocorr'); ylim([0 0.03]);
        title('autocorr BP-R','Color',condicolors(condi,:));
end
suptitle(['BPcorr-boolean ' 'subject ' num2str(seed)]);
        
figureName=['BPcorr_boolean -- ' num2str(seed)];
    % save the figure
    saveas(gcf,figureName,'jpg');

%% PLOT 5:  BP interval spectra DFA  
% Pspectra(BP_L) + Pspectra(BP_R) + DFA(BP_L) + DFA(BP_R)
% 6 conditions -subplots(4,6,i)
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/DFA

figure('units','normalized','outerposition',[0 0 1 1]);
for condi=1:4
    if condi==1; BP_L=BPCondi1L; BP_R=BPCondi1R; end
    if condi==2; BP_L=BPCondi2L; BP_R=BPCondi2R; end
    if condi==3; BP_L=BPCondi3L; BP_R=BPCondi3R; end
    if condi==4; BP_L=BPCondi4L; BP_R=BPCondi4R; end
    subplot(4,4,condi); % Pspectra(BP_L)
        y=Calinterval(BP_L');
        if (condi == 1) | (condi == 2) |(condi == 3) | (condi == 4);
            Fs=2;% assuming a sampling frequency of 2 Hz
        else (condi == 5) | (condi == 6) ;
            Fs=3;% assuming a sampling frequency of 3 Hz
        end
        [freqs,fcoef] = oneoverf(y,Fs);
        xlabel('Log(f)');ylabel('Log(power)');
        title({['Condi ' conditionNames{condi}], 'Spectrum BP-L'},'Color',condicolors(condi,:));
        ylim([-2 8]);xlim([-3 0]);
        
    subplot(4,4,4+condi);  % Pspectra(BP_R) 
        y=Calinterval(BP_R');
        if (condi == 1) | (condi == 2) |(condi == 3) | (condi == 4);
            Fs=2;% assuming a sampling frequency of 2 Hz
        else (condi == 5) | (condi == 6) ;
            Fs=3;% assuming a sampling frequency of 3 Hz
        end
        [freqs,fcoef] = oneoverf(y,Fs);
        xlabel('Log(f)');ylabel('Log(power)');
        title('Spectrum BP-R','Color',condicolors(condi,:));  
        ylim([-2 8]);xlim([-3 0]);
        
    subplot(4,4,4*2+condi);  % DFA(BP_L)
        [D,Alpha1]=DFA_main(Calinterval(BP_L'));
        title('DFA BP-L','Color',condicolors(condi,:));
        ylim([0.5 4]);xlim([0.5 3]);
        
    subplot(4,4,4*3+condi);  % DFA(BP_R)
        [D,Alpha1]=DFA_main(Calinterval(BP_R'));
        title('DFA BP-R','Color',condicolors(condi,:));
        ylim([0.5 4]);xlim([0.5 3]);

end
suptitle(['BP intervals ' ' -- subject ' num2str(seed)]);

figureName=['BP-spectra-DFA -- ' num2str(seed)];
    % save the figure
    saveas(gcf,figureName,'jpg');
   
%% PLOT 5-1:  BP interval spectra DFA   
% Pspectra(BP_L) + Pspectra(BP_R) + DFA(BP_L) + DFA(BP_R)
% 6 conditions -subplots(4,6,i)
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/DFA

figure('units','normalized','outerposition',[0 0 1 1]);
for condi=1:4
    if condi==1; BP_L=BPCondi1L; BP_R=BPCondi1R; end
    if condi==2; BP_L=BPCondi2L; BP_R=BPCondi2R; end
    if condi==3; BP_L=BPCondi3L; BP_R=BPCondi3R; end
    if condi==4; BP_L=BPCondi4L; BP_R=BPCondi4R; end
    subplot(4,4,condi); % Pspectra(BP_L)
        y=Calinterval(BP_L');
        if (condi == 1) | (condi == 2) |(condi == 3) | (condi == 4);
            Fs=2;% assuming a sampling frequency of 2 Hz
        else (condi == 5) | (condi == 6) ;
            Fs=3;% assuming a sampling frequency of 3 Hz
        end
        [freqs,fcoef] = oneoverf(y,Fs);
        xlabel('Log(f)');ylabel('Log(power)');
        title({['Condi ' conditionNames{condi}], 'Spectrum BP-L'},'Color',condicolors(condi,:));
        ylim([-2 8]);xlim([-3 0]);
        
    subplot(4,4,4+condi);  % Pspectra(BP_R) 
        y=Calinterval(BP_R');
        if (condi == 1) | (condi == 2) |(condi == 3) | (condi == 4);
            Fs=2;% assuming a sampling frequency of 2 Hz
        else (condi == 5) | (condi == 6) ;
            Fs=3;% assuming a sampling frequency of 3 Hz
        end
        [freqs,fcoef] = oneoverf(y,Fs);
        xlabel('Log(f)');ylabel('Log(power)');
        title('Spectrum BP-R','Color',condicolors(condi,:));  
        ylim([-2 8]);xlim([-3 0]);
        
    subplot(4,4,4*2+condi);  % DFA(BP_L)
        [D,Alpha1]=DFA_main(Calinterval(BP_L'));
        title('DFA BP-L','Color',condicolors(condi,:));
        ylim([0.5 4]);xlim([0.5 3]);
        
    subplot(4,4,4*3+condi);  % DFA(BP_R)
        [D,Alpha1]=DFA_main(Calinterval(BP_R'));
        title('DFA BP-R','Color',condicolors(condi,:));
        ylim([0.5 4]);xlim([0.5 3]);

end
suptitle(['BP-interval-spectra-DFA-first10-' ' -- subject ' num2str(seed) '-PLOT-5-1']);

figureName=['BP-interval-spectra-DFA-first10- ' num2str(seed) '-PLOT-5-1'];
    % save the figure
    saveas(gcf,figureName,'jpg');
    saveas(gcf,figureName,'fig');
 
%% PLOT 6:  BP interval spectra-DFA smoothed 
% Pspectra(BP_L) + Pspectra(BP_R) + DFA(BP_L) + DFA(BP_R)
% 6 conditions -subplots(4,6,i)
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing

figure('units','normalized','outerposition',[0 0 1 1]);
win=20; % smooth with a sliding window size
for condi=1:4
    if condi==1; BP_L=BPCondi1L; BP_R=BPCondi1R; end
    if condi==2; BP_L=BPCondi2L; BP_R=BPCondi2R; end
    if condi==3; BP_L=BPCondi3L; BP_R=BPCondi3R; end
    if condi==4; BP_L=BPCondi4L; BP_R=BPCondi4R; end
    subplot(4,4,condi); % Pspectra(BP_L)
        y=smoothing(Calinterval(BP_L'),win);
        if (condi == 1) | (condi == 2) |(condi == 3) | (condi == 4);
            Fs=2;% assuming a sampling frequency of 2 Hz
        else (condi == 5) | (condi == 6) ;
            Fs=3;% assuming a sampling frequency of 3 Hz
        end
        [freqs,fcoef] = oneoverf(y,Fs);
        xlabel('Log(f)');ylabel('Log(power)');
        title({['Condi ' conditionNames{condi}], 'Spectrum BP_L'},'Color',condicolors(condi,:));
        ylim([-4 8]);xlim([-3 0]);
        
    subplot(4,4,4+condi);  % Pspectra(BP_R) 
        y=smoothing(Calinterval(BP_R'),win);
        if (condi == 1) | (condi == 2) |(condi == 3) | (condi == 4);
            Fs=2;% assuming a sampling frequency of 2 Hz
        else (condi == 5) | (condi == 6) ;
            Fs=3;% assuming a sampling frequency of 3 Hz
        end
        [freqs,fcoef] = oneoverf(y,Fs);
        xlabel('Log(f)');ylabel('Log(power)');
        title('Spectrum BP_R', 'Color',condicolors(condi,:));    
        ylim([-4 8]);xlim([-3 0]);
        
    subplot(4,4,4*2+condi);  % DFA(BP_L)
        [D,Alpha1]=DFA_main(smoothing(Calinterval(BP_L'),win));
        title('DFA BP-L','Color',condicolors(condi,:));
        ylim([0 4]);xlim([0.5 3]);
        
    subplot(4,4,4*3+condi);  % DFA(BP_R)
        [D,Alpha1]=DFA_main(smoothing(Calinterval(BP_R'),win));
        title('DFA BP-R','Color',condicolors(condi,:));
        ylim([0 4]);xlim([0.5 3]);

end
suptitle(['BP-intervals-smoothed win ' num2str(win) ' subject ' num2str(seed)]);

figureName=['BP-spectra-DFA_smooth-- ' num2str(seed)];
    % save the figure
    saveas(gcf,figureName,'jpg');

%% PLOT 6-1:  BP interval spectra-DFA smoothed 
% Pspectra(BP_L) + Pspectra(BP_R) + DFA(BP_L) + DFA(BP_R)
% 6 conditions -subplots(4,6,i)
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing
% select number of conditions in the experiment
Totalcondi=4;
Totalcondi=5;
Totalcondi=6;

figure('units','normalized','outerposition',[0 0 1 1]);
win=13; % smooth with a sliding window size
for condi=1:Totalcondi
    if condi==1; BP_L=BPCondi1L; BP_R=BPCondi1R; end
    if condi==2; BP_L=BPCondi2L; BP_R=BPCondi2R; end
    if condi==3; BP_L=BPCondi3L; BP_R=BPCondi3R; end
    if condi==4; BP_L=BPCondi4L; BP_R=BPCondi4R; end
    if condi==5; BP_L=BPCondi5L; BP_R=BPCondi5R; end
    if condi==6; BP_L=BPCondi6L; BP_R=BPCondi6R; end
    subplot(4,Totalcondi,condi); % Pspectra(BP_L)
        y=smoothing(Calinterval(BP_L'),win);
        if (condi == 1) | (condi == 2) |(condi == 3) | (condi == 4);
            Fs=1.3;% assuming a sampling frequency of 2 Hz
        else (condi == 5) | (condi == 6) ;
            Fs=3;% assuming a sampling frequency of 3 Hz
        end
        [freqs,fcoef] = oneoverf(y,Fs);
        xlabel('Log(f)');ylabel('Log(power)');
        title({['Condi ' conditionNames{condi}], 'Spectrum BP_L'},'Color',condicolors(condi,:));
        ylim([-3 3]);xlim([-2.5 0]);
        
    subplot(4,Totalcondi,Totalcondi+condi);  % Pspectra(BP_R) 
        y=smoothing(Calinterval(BP_R'),win);
        if (condi == 1) | (condi == 2) |(condi == 3) | (condi == 4);
            Fs=1.3;% assuming a sampling frequency of 2 Hz
        else (condi == 5) | (condi == 6) ;
            Fs=3;% assuming a sampling frequency of 3 Hz
        end
        [freqs,fcoef] = oneoverf(y,Fs);
        xlabel('Log(f)');ylabel('Log(power)');
        title('Spectrum BP_R', 'Color',condicolors(condi,:));    
        ylim([-3 3]);xlim([-2.5 0]);
        
    subplot(4,Totalcondi,Totalcondi*2+condi);  % DFA(BP_L)
        [D,Alpha1]=DFA_main(smoothing(Calinterval(BP_L'),win));
        title('DFA BP-L','Color',condicolors(condi,:));
        ylim([0 4]);xlim([0 4]);
        
    subplot(4,Totalcondi,Totalcondi*3+condi);  % DFA(BP_R)
        [D,Alpha1]=DFA_main(smoothing(Calinterval(BP_R'),win));
        title('DFA BP-R','Color',condicolors(condi,:));
        ylim([0 4]);xlim([0 4]);

end
suptitle(['BP-intervals-smoothed win ' num2str(win) ' subject ' num2str(seed) '-PLOT-6-1']);

figureName=['BP-spectra-DFA_smooth-- ' num2str(seed) '-PLOT-6-1'];
    % save the figure
    saveas(gcf,figureName,'fig');
    
%% PLOT 7:  (give up) BP boolean spectra DFA  (long-range / global statistics; strong anticipation) (very slow, basically white noise, no long range dependence)
% Pspectra(BP_L) + Pspectra(BP_R) + DFA(BP_L) + DFA(BP_R)
% 6 conditions -subplots(4,6,i)
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing
tic
figure('units','normalized','outerposition',[0 0 1 1]);
for condi=1:4
    if condi==1; BP_L=BPCondi1L; BP_R=BPCondi1R; end
    if condi==2; BP_L=BPCondi2L; BP_R=BPCondi2R; end
    if condi==3; BP_L=BPCondi3L; BP_R=BPCondi3R; end
    if condi==4; BP_L=BPCondi4L; BP_R=BPCondi4R; end
    subplot(4,4,condi); % Pspectra(BP_L)
        y=BP_L';
        if (condi == 1) | (condi == 2) |(condi == 3) | (condi == 4);
            Fs=1.3;% assuming a sampling frequency of 2 Hz
        end
        [freqs,fcoef] = oneoverf(y,Fs);
        xlabel('Log(f)');ylabel('Log(power)');
        title({['Condi ' conditionNames{condi}],'Spectrum BP_L'});
        
    subplot(4,4,4+condi);  % Pspectra(BP_R) 
        y=BP_R';
        if (condi == 1) | (condi == 2) |(condi == 3) | (condi == 4);
            Fs=1.3;% assuming a sampling frequency of 2 Hz
        end
        [freqs,fcoef] = oneoverf(y,Fs);
        xlabel('Log(f)');ylabel('Log(power)');
        title({'Spectrum BP_R', ['Condi ' num2str(condi)]});    
        
    subplot(4,4,4*2+condi);  % DFA(BP_L)
        [D,Alpha1]=DFA_main(BP_L');
        title({'DFA BP_L', ['Condi ' num2str(condi)]});
        
    subplot(4,4,4*3+condi);  % DFA(BP_R)
        [D,Alpha1]=DFA_main(BP_R');
        title({'DFA BP_R', ['Condi ' num2str(condi)]});

end
suptitle({'BP-boolean-spectra-DFA -- ',['subject ' num2str(seed)]});
toc

figureName=['BP_boolean-spectra-DFA -- ' num2str(seed)];
    % save the figure
    saveas(gcf,figureName,'jpg');

%% PLOT 8:  P2 Error smooth (synchronization)
% Error + autocorr (Error) + Pspectra(Error)+ DFA(Error)
% uncoupled condition (bidiretional error) -subplots(4,10,1-2)
% unidrectional conditions -subplots(4,10,3-4)
% bidirectional condition (bidiretional error) -subplots(4,10,5-10)
addpath D:\360MoveData\Users\alienware\Documents\GitHub\Motor_cordination\1_over_f\data_analysis\2P_testing
addpath D:\360MoveData\Users\alienware\Documents\GitHub\Motor_cordination\1_over_f\data_analysis\DFA

figure('units','normalized','outerposition',[0 0 1 1]);
win=13; % smooth the error with sliding window size

% select condition to plot
condiSelect=[1 2 3 4 5 6 7 8];
Totalcondi=length(condiSelect);

for c=1:Totalcondi; % 1:Totalcondi*2
    condi=condiSelect(c);
    if condi==1; BP_L=BPCondi1L; FB_L=FBCondi1L; time_series1=find(FB_L==1); time_series2=find(BP_L==1); subcondi=1; end % R's Error base on L (uncouple)
    if condi==2; BP_R=BPCondi1R; FB_R=FBCondi1R; time_series1=find(FB_R==1); time_series2=find(BP_R==1); subcondi=1; end % L's Error base on R (uncouple)
    if condi==3; BP_L=BPCondi2L; FB_L=FBCondi2L; time_series1=find(FB_L==1); time_series2=find(BP_L==1); subcondi=2; end % R's Error base on L (unidirectional)
    if condi==4; BP_R=BPCondi2R; FB_R=FBCondi2R; time_series1=find(FB_R==1); time_series2=find(BP_R==1); subcondi=2; end % L's Error base on R (unidirectional)
    if condi==5; BP_L=BPCondi3L; FB_L=FBCondi3L; time_series1=find(FB_L==1); time_series2=find(BP_L==1); subcondi=3; end % R's Error base on L (unidirectional)
    if condi==6; BP_R=BPCondi3R; FB_R=FBCondi3R; time_series1=find(FB_R==1); time_series2=find(BP_R==1); subcondi=3; end % L's Error base on R (unidirectional)
    if condi==7; BP_L=BPCondi4L; FB_L=FBCondi4L; time_series1=find(FB_L==1); time_series2=find(BP_L==1); subcondi=4; end % R's Error base on L (2Hz bidirection)
    if condi==8; BP_R=BPCondi4R; FB_R=FBCondi4R; time_series1=find(FB_R==1); time_series2=find(BP_R==1); subcondi=4; end % L's Error base on R (2Hz bidirection)
    if condi==9; BP_L=BPCondi5L; FB_L=FBCondi5L; time_series1=find(FB_L==1); time_series2=find(BP_L==1); subcondi=5; end % R's Error base on L (3Hz bidirection)
    if condi==10; BP_R=BPCondi5R; FB_R=FBCondi5R; time_series1=find(FB_R==1); time_series2=find(BP_R==1); subcondi=5; end % L's Error base on R (3Hz bidirection)
    if condi==11; BP_L=BPCondi6L; FB_L=FBCondi6L; time_series1=find(FB_L==1); time_series2=find(BP_L==1); subcondi=6; end % R's Error base on L (>3Hz bidirection)
    if condi==12; BP_R=BPCondi6R; FB_R=FBCondi6R; time_series1=find(FB_R==1); time_series2=find(BP_R==1); subcondi=6; end % L's Error base on R (>3Hz bidirection)% Compute mid-timepoints for time_series1
    TimeLength=length(time_series1)-1;
    Error=[];
    % find the time difference with the closest botton press 
    for i=1:TimeLength % i=StartStim2 % syncopation has one press less than condition 1 and 3
        [minValue,closetIndex]=min(abs(time_series2-time_series1(i))); % closetIndex in time_series2
        Error(i)=time_series2(closetIndex)-time_series1(i);
    end
    % smooth the error with sliding window
    [ErrorSmooth] = smoothing(Error,win); % figure;plot(ErrorSmooth);
%     ErrorSmooth=Error;

    subplot(4,Totalcondi,c); % Error
        plot(ErrorSmooth/2000,'r.');
        xlabel('taps');ylabel('timing error (s)');
        title({['Condi ' conditionNames{subcondi}],'error'},'Color',condicolors(subcondi,:));
        % ylim([-0.5 0.6]);
        
    subplot(4,Totalcondi,Totalcondi+c);  % autocorr (Error) 
        autocorr(ErrorSmooth/2000,400); xlabel('lags');ylabel('autocorr'); 
        title('error','Color',condicolors(subcondi,:));
        % ylim([-0.2 0.8]);
        
    subplot(4,Totalcondi,Totalcondi*2+c);  % Pspectra(Error)
        y=ErrorSmooth/2000;
        if (condi == 1) | (condi == 2) |(condi == 3) | (condi == 4) | (condi == 5) | (condi == 6) | (condi == 7) | (condi == 8);
            Fs=1.3;% assuming a sampling frequency of 2 Hz
        else
            Fs=3;
        end
        [freqs,fcoef] = oneoverf(y,Fs);
        xlabel('Log(f)');ylabel('Log(power)');
        title('Spectrum Error','Color',condicolors(subcondi,:));
        % ylim([-6 -2]);
        
    subplot(4,Totalcondi,Totalcondi*3+c);  %  DFA(Error) 
        [D,Alpha1]=DFA_main(ErrorSmooth/2000);
        title('DFA Error', 'Color',condicolors(subcondi,:));
        % ylim([-2 0.5]);

end
suptitle(['error-smooth ' 'win ' num2str(win) ' subject ' num2str(seed) '-PLOT-8-synchronization']);

figureName=['Error-smooth-win-' num2str(win) '--' num2str(seed) '-PLOT-8-synchronization'];
    % save the figure
    saveas(gcf,figureName,'fig');

%% PLOT 8-1:  P2 Error-Pacer smooth (syncopation)
% Error + autocorr (Error) + Pspectra(Error)+ DFA(Error)
% uncoupled condition (bidiretional error) -subplots(4,10,1-2)
% unidrectional conditions -subplots(4,10,3-4)
% bidirectional condition (bidiretional error) -subplots(4,10,5-10)
addpath D:\360MoveData\Users\alienware\Documents\GitHub\Motor_cordination\1_over_f\data_analysis\2P_testing
addpath D:\360MoveData\Users\alienware\Documents\GitHub\Motor_cordination\1_over_f\data_analysis\DFA

figure('units','normalized','outerposition',[0 0 1 1]);
win=13; % smooth the error with sliding window size

% select condition to plot
condiSelect=[1 2 3 4 5 6 7 8];
Totalcondi=length(condiSelect);

for c=1:Totalcondi; % 1:Totalcondi*2
    condi=condiSelect(c);
%     if condi==1; BP_L=BPCondi1L; BP_R=BPCondi1R; time_series1=find(BP_L==1); time_series2=find(BP_R==1); subcondi=1; end % R's Error base on L (uncouple)
%     if condi==2; BP_L=BPCondi1L; BP_R=BPCondi1R; time_series1=find(BP_R==1); time_series2=find(BP_L==1); subcondi=1; end % L's Error base on R (uncouple)
%     if condi==3; BP_L=BPCondi2L; BP_R=BPCondi2R; time_series1=find(BP_L==1); time_series2=find(BP_R==1); subcondi=2; end % R's Error base on L (unidirectional)
%     if condi==4; BP_L=BPCondi2L; BP_R=BPCondi2R; time_series1=find(BP_R==1); time_series2=find(BP_L==1); subcondi=2; end % L's Error base on R (unidirectional)
%     if condi==5; BP_L=BPCondi3L; BP_R=BPCondi3R; time_series1=find(BP_L==1); time_series2=find(BP_R==1); subcondi=3; end % R's Error base on L (unidirectional)
%     if condi==6; BP_L=BPCondi3L; BP_R=BPCondi3R; time_series1=find(BP_R==1); time_series2=find(BP_L==1); subcondi=3; end % L's Error base on R (unidirectional)
%     if condi==7; BP_L=BPCondi4L; BP_R=BPCondi4R; time_series1=find(BP_L==1); time_series2=find(BP_R==1); subcondi=4; end % R's Error base on L (2Hz bidirection)
%     if condi==8; BP_L=BPCondi4L; BP_R=BPCondi4R; time_series1=find(BP_R==1); time_series2=find(BP_L==1); subcondi=4; end % L's Error base on R (2Hz bidirection)
    if condi==1; BP_L=BPCondi1L; FB_L=FBCondi1L; time_series1=find(FB_L==1); time_series2=find(BP_L==1); subcondi=1; Pacer=find(PacerCondi1L==1); end % R's Error base on L (uncouple)
    if condi==2; BP_R=BPCondi1R; FB_R=FBCondi1R; time_series1=find(FB_R==1); time_series2=find(BP_R==1); subcondi=1; Pacer=find(PacerCondi1R==1); end % L's Error base on R (uncouple)
    if condi==3; BP_L=BPCondi2L; FB_L=FBCondi2L; time_series1=find(FB_L==1); time_series2=find(BP_L==1); subcondi=2; Pacer=find(PacerCondi2L==1); end % R's Error base on L (unidirectional)
    if condi==4; BP_R=BPCondi2R; FB_R=FBCondi2R; time_series1=find(FB_R==1); time_series2=find(BP_R==1); subcondi=2; Pacer=find(PacerCondi2R==1); end % L's Error base on R (unidirectional)
    if condi==5; BP_L=BPCondi3L; FB_L=FBCondi3L; time_series1=find(FB_L==1); time_series2=find(BP_L==1); subcondi=3; Pacer=find(PacerCondi3L==1); end % R's Error base on L (unidirectional)
    if condi==6; BP_R=BPCondi3R; FB_R=FBCondi3R; time_series1=find(FB_R==1); time_series2=find(BP_R==1); subcondi=3; Pacer=find(PacerCondi3R==1); end % L's Error base on R (unidirectional)
    if condi==7; BP_L=BPCondi4L; FB_L=FBCondi4L; time_series1=find(FB_L==1); time_series2=find(BP_L==1); subcondi=4; Pacer=find(PacerCondi4L==1); end % R's Error base on L (2Hz bidirection)
    if condi==8; BP_R=BPCondi4R; FB_R=FBCondi4R; time_series1=find(FB_R==1); time_series2=find(BP_R==1); subcondi=4; Pacer=find(PacerCondi4R==1); end % L's Error base on R (2Hz bidirection)
    % Compute mid-timepoints from Pacer
    if condi==7 | condi==8
        % Compute mid-timepoints for time_series1
        TimeLength=length(time_series1)-1;
        mid_series1=[];
        for i=1:TimeLength
            mid_series1(i)=(time_series1(i)+time_series1(i+1))/2;
        end
    else
        TimeLength=length(Pacer)-1;
        mid_series1=[];
        for i=1:TimeLength
            mid_series1(i)=(Pacer(i)+Pacer(i+1))/2;
        end
    end
    Error=[];
    % find the time difference with the closest botton press 
    for i=1:TimeLength % i=StartStim2 % syncopation has one press less than condition 1 and 3
        [minValue,closetIndex]=min(abs(time_series2-mid_series1(i))); % closetIndex in time_series2
        Error(i)=time_series2(closetIndex)-mid_series1(i);
    end
    % smooth the error with sliding window
    [ErrorSmooth] = smoothing(Error,win); % figure;plot(ErrorSmooth);
%     ErrorSmooth=Error;

    subplot(4,Totalcondi,c); % Error
        plot(ErrorSmooth/2000,'r.');
        xlabel('taps');ylabel('timing error (s)');
        title({['Condi ' conditionNames{subcondi}],'error'},'Color',condicolors(subcondi,:));
        % ylim([-0.5 0.6]);
        
    subplot(4,Totalcondi,Totalcondi+c);  % autocorr (Error) 
        autocorr(ErrorSmooth/2000,350); xlabel('lags');ylabel('autocorr'); 
        title('error','Color',condicolors(subcondi,:));
        % ylim([-0.2 0.8]);
        
    subplot(4,Totalcondi,Totalcondi*2+c);  % Pspectra(Error)
        y=ErrorSmooth/2000;
        if (condi == 1) | (condi == 2) |(condi == 3) | (condi == 4) | (condi == 5) | (condi == 6) | (condi == 7) | (condi == 8);
            Fs=1.3;% assuming a sampling frequency of 2 Hz
        end
        [freqs,fcoef] = oneoverf(y,Fs);
        xlabel('Log(f)');ylabel('Log(power)');
        title('Spectrum Error','Color',condicolors(subcondi,:));
        % ylim([-5 -3]);
        
    subplot(4,Totalcondi,Totalcondi*3+c);  %  DFA(Error) 
        [D,Alpha1]=DFA_main(ErrorSmooth/2000);
        title('DFA Error', 'Color',condicolors(subcondi,:));
        % ylim([-2 0.5]);

end
suptitle(['error-pacer-smooth ' 'win ' num2str(win) ' subject ' num2str(seed) '-PLOT-8-1-syncopation']);

figureName=['Error-pacer-smooth-win' num2str(win) '--subject--' num2str(seed) '-PLOT-8-1-syncopation'];
% save the figure
saveas(gcf,figureName,'fig');
    
%% PLOT 9:  P2 Error smooth normalized (short & long -range statistics)
% Error + autocorr (Error) + Pspectra(Error)+ DFA(Error)
% uncoupled condition (bidiretional error) -subplots(4,10,1-2)
% unidrectional conditions -subplots(4,10,3-4)
% bidirectional condition (bidiretional error) -subplots(4,10,5-10)
addpath D:\360MoveData\Users\alienware\Documents\GitHub\Motor_cordination\1_over_f\data_analysis\2P_testing
addpath D:\360MoveData\Users\alienware\Documents\GitHub\Motor_cordination\1_over_f\data_analysis\DFA

figure('units','normalized','outerposition',[0 0 1 1]);
win=5; % smooth the error with sliding window size
for condi=1:8
    if condi==1; BP_L=BPCondi1L; FB_L=FBCondi1L; time_series1=find(BP_L==1); time_series2=find(FB_L==1); subcondi=1; end % R's Error base on L (uncouple)
    if condi==2; BP_R=BPCondi1R; FB_R=FBCondi1R; time_series1=find(BP_R==1); time_series2=find(FB_R==1); subcondi=1; end % L's Error base on R (uncouple)
    if condi==3; BP_L=BPCondi2L; FB_L=FBCondi2L; time_series1=find(BP_L==1); time_series2=find(FB_L==1); subcondi=2; end % R's Error base on L (unidirectional)
    if condi==4; BP_R=BPCondi2R; FB_R=FBCondi2R; time_series1=find(BP_R==1); time_series2=find(FB_R==1); subcondi=2; end % L's Error base on R (unidirectional)
    if condi==5; BP_L=BPCondi3L; FB_L=FBCondi3L; time_series1=find(BP_L==1); time_series2=find(FB_L==1); subcondi=3; end % R's Error base on L (unidirectional)
    if condi==6; BP_R=BPCondi3R; FB_R=FBCondi3R; time_series1=find(BP_R==1); time_series2=find(FB_R==1); subcondi=3; end % L's Error base on R (unidirectional)
    if condi==7; BP_L=BPCondi4L; FB_L=FBCondi4L; time_series1=find(BP_L==1); time_series2=find(FB_L==1); subcondi=4; end % R's Error base on L (2Hz bidirection)
    if condi==8; BP_R=BPCondi4R; FB_R=FBCondi4R; time_series1=find(BP_R==1); time_series2=find(FB_R==1); subcondi=4; end % L's Error base on R (2Hz bidirection)
    % Compute mid-timepoints for time_series1
    TimeLength=length(time_series1)-1;
    mid_series1=[];
    for i=1:TimeLength
        mid_series1(i)=(time_series1(i)+time_series1(i+1))/2;
    end
    Error=[];
    % find the time difference with the closest botton press 
    for i=1:TimeLength % i=StartStim2 % syncopation has one press less than condition 1 and 3
        [minValue,closetIndex]=min(abs(time_series2-mid_series1(i))); % closetIndex in time_series2
        Error(i)=time_series2(closetIndex)-mid_series1(i);
    end
    % smooth the error with sliding window
    [ErrorSmooth] = smoothing(Error,win); % figure;plot(ErrorSmooth);
%     ErrorSmooth=Error;
    % normalize
    ErrorNorm=normalize(ErrorSmooth,'center','median');% figure;plot(ErrorSmooth);
    
    subplot(4,8,condi); % Error
        plot(ErrorNorm,'r.');
        xlabel('taps');ylabel('error');
        title({['Condi ' conditionNames{subcondi}],'error'},'Color',condicolors(subcondi,:));
        % ylim([-0.5 0.6]);
        
    subplot(4,8,8+condi);  % autocorr (Error) 
        autocorr(ErrorNorm,400); xlabel('lags');ylabel('autocorr'); 
        title('error','Color',condicolors(subcondi,:));
        % ylim([-0.2 0.8]);
        
    subplot(4,8,8*2+condi);  % Pspectra(Error)
        y=ErrorNorm;
        if (condi == 1) | (condi == 2) |(condi == 3) | (condi == 4) | (condi == 5) | (condi == 6) | (condi == 7) | (condi == 8);
            Fs=1.3;% assuming a sampling frequency of 2 Hz
        end
        [freqs,fcoef] = oneoverf(y,Fs);
        xlabel('Log(f)');ylabel('Log(power)');
        title('Spectrum Error','Color',condicolors(subcondi,:));
        % ylim([-10 0]);
        
    subplot(4,8,8*3+condi);  %  DFA(Error) 
        [D,Alpha1]=DFA_main(ErrorNorm);
        title('DFA Error', 'Color',condicolors(subcondi,:));
        % ylim([-2 0.5]);

end
suptitle(['error-normalized ' 'win ' num2str(win) ' subject ' num2str(seed)]);

figureName=['Error-normalized -- ' num2str(seed)];
    % save the figure
    saveas(gcf,figureName,'jpg');
    saveas(gcf,figureName,'fig');
    close all;
    
%% PLOT 10:  EMG corr (short-range / local statistics; weak anticipation)
% [ autocorr(EMG_L) + xcorr (EMG_L vs EMG_R) + autocorr(EMG_R) ]
% 6 conditions - subplots(3,6,i)
% select number of conditions in the experiment
Totalcondi=4;
Totalcondi=5;
Totalcondi=6;

figure('units','normalized','outerposition',[0 0 1 1]);
for condi=1:Totalcondi
    if condi==1; EMG_L=EMGCondi1L; EMG_R=EMGCondi1R; end % R's Error base on L (uncouple)
    if condi==2; EMG_L=EMGCondi2L; EMG_R=EMGCondi2R; end % L's Error base on R (uncouple)
    if condi==3; EMG_L=EMGCondi3L; EMG_R=EMGCondi3R; end % R's Error base on L (unidirectional)
    if condi==4; EMG_L=EMGCondi4L; EMG_R=EMGCondi4R; end % L's Error base on R (unidirectional)
    if condi==5; EMG_L=EMGCondi5L; EMG_R=EMGCondi5R; end % R's Error base on L (2Hz bidirection)
    if condi==6; EMG_L=EMGCondi6L; EMG_R=EMGCondi6R; end % L's Error base on R (2Hz bidirection)
    
    subplot(3,Totalcondi,condi); % autocorr(EMG_L)
        lags=1000; % round(length(EMG_L))-1;
        autocorr(EMG_L',lags); xlabel('lags of samples');ylabel('autocorr'); 
        title({['Condi ' conditionNames{condi}],'EMG-L'},'Color',condicolors(condi,:));
        % ylim([-0.5 1]);
        % xlim([0 2000]);
        
    subplot(3,Totalcondi,Totalcondi+condi); % xcorr (EMG_L vs EMG_R) 
        time_series1=EMG_L';% plot(time_series1);
        time_series2=EMG_R';% plot(time_series2);
        time_series1=time_series1-mean(time_series1);
        time_series2=time_series2-mean(time_series2);
        TimeLength=1000; % min([length(time_series1) length(time_series2)]);
        [r,lags]=xcorr(time_series1(1:TimeLength), time_series2(1:TimeLength), TimeLength-1,'normalized');
        plot(lags./2,r);xlabel('lag of time [ms]');ylabel('Xcorr');
        title('EMG-L & EMG-R','Color',condicolors(condi,:));
        xline(0,'--','Color',deepyellow);
        ylim([0 1]);% xlim([-5e5 5e5]);
        % xlim([-2000 2000]);
        
    subplot(3,Totalcondi,Totalcondi*2+condi); %  autocorr(EMG_R) 
        lags=1000; % round(length(EMG_R))-1;
        autocorr(EMG_R',lags); xlabel('lags of samples');ylabel('autocorr'); 
        title('EMG-R','Color',condicolors(condi,:));
        % ylim([-0.5 1]);
        % xlim([0 2000]);
end
suptitle(['EMG corr ' '-- subject ' num2str(seed)]);

figureName=['EMGcorr' num2str(seed)];
    % save the figure
    saveas(gcf,figureName,'jpg');

%% PLOT 11:  EMG spectra (long-range / global statistics; strong anticipation)
% [ Pspectra(EMG_L) + Xspectra (EMG_L vs EMG_R) + Pspectra(EMG_R) ]
% 6 conditions - subplots(3,6,i)
% select number of conditions in the experiment
Totalcondi=4;
Totalcondi=5;

maxfreq=10;sr=2000;
figure('units','normalized','outerposition',[0 0 1 1]);
for condi=1:Totalcondi
    if condi==1; EMG_L=EMGCondi1L; EMG_R=EMGCondi1R; end % R's Error base on L (uncouple)
    if condi==2; EMG_L=EMGCondi2L; EMG_R=EMGCondi2R; end % L's Error base on R (uncouple)
    if condi==3; EMG_L=EMGCondi3L; EMG_R=EMGCondi3R; end % R's Error base on L (unidirectional)
    if condi==4; EMG_L=EMGCondi4L; EMG_R=EMGCondi4R; end % L's Error base on R (unidirectional)
    if condi==5; EMG_L=EMGCondi5L; EMG_R=EMGCondi5R; end % R's Error base on L (2Hz bidirection)
    if condi==6; EMG_L=EMGCondi6L; EMG_R=EMGCondi6R; end % L's Error base on R (2Hz bidirection)
    
    subplot(3,Totalcondi,condi); % Pspectra(EMG_L) 
        [fcoef1,fcoef2,cprod, freqs] = spectra(EMG_L,EMG_L,maxfreq,sr);
        plot(freqs,abs(cprod));xlabel('freqs');ylabel('power');
        title({['Condi ' conditionNames{condi}], 'EMG-L'},'Color',condicolors(condi,:));
        xticks(linspace(0,10,5));xticklabels({'0','2.5','5','7.5','10'});
%         ylim([0 3e-3]);
        
    subplot(3,Totalcondi,Totalcondi+condi); % Xspectra (EMG_L vs EMG_R) 
        [fcoef1,fcoef2,cprod, freqs] = spectra(EMG_L,EMG_R,maxfreq,sr);
        plot(freqs,abs(cprod));xlabel('freqs');ylabel('xspectra');
        title('EMG-L & EMG-R','Color',condicolors(condi,:));
        xticks(linspace(0,10,5));xticklabels({'0','2.5','5','7.5','10'});
%         ylim([0 6]);
        
    subplot(3,Totalcondi,Totalcondi*2+condi); %  Pspectra (EMG_R) 
        [fcoef1,fcoef2,cprod, freqs] = spectra(EMG_R,EMG_R,maxfreq,sr);
        plot(freqs,abs(cprod));xlabel('freqs');ylabel('power');
        title('EMG-R','Color',condicolors(condi,:));
        xticks(linspace(0,10,5));xticklabels({'0','2.5','5','7.5','10'});
%         ylim([0 2e-3]);
end
suptitle(['EMG-spectra -- ' 'subject ' num2str(seed)]);

figureName=['EMG-spectra' num2str(seed)];
    % save the figure
    saveas(gcf,figureName,'jpg');

%% PLOT 12:  EEG auto/xcorr single chan (weak anticipation, short-range / local statistics)
% [ autocorr(EEG_L) + xcorr(EEG_L vs EEG_R) + autocorr(EEG_R) ]
% each of the 6 conditions (subplots(3,6,i)) for each of the 32 channels;
% auto save the 32 figures in a folder named - ShortEEG
cd /ssd/zhibin/1overf/20220515_2P/Segmented_data/Plots/EEGcorr_chans
cd /ssd/zhibin/1overf/20220517_2P/Segmented_data/Plots/ShortEEG
cd /ssd/zhibin/1overf/20220518_2P/Segmented_data/Plots/EEGcorr_chans

cd /ssd/zhibin/1overf/20220610_2P/Segmented_data/Plots/EEGcorr_chans

for chan = 1:32
    figure('units','normalized','outerposition',[0 0 1 1]);
    for condi=1:4
        if condi==1; EEG_L=EEGCondi1L; EEG_R=EEGCondi1R; end % R's Error base on L (uncouple)
        if condi==2; EEG_L=EEGCondi2L; EEG_R=EEGCondi2R; end % L's Error base on R (uncouple)
        if condi==3; EEG_L=EEGCondi3L; EEG_R=EEGCondi3R; end % R's Error base on L (unidirectional)
        if condi==4; EEG_L=EEGCondi4L; EEG_R=EEGCondi4R; end % L's Error base on R (unidirectional)

        subplot(3,4,condi); % autocorr(EEG_L) 
            lags=1000;% round(length(EEG_L(:,chan)))-1;
            autocorr(EEG_L(:,chan)',lags); xlabel('lags of samples');ylabel('autocorr'); 
            title({['Condi ' conditionNames{condi}],'EEG-L'},'Color',condicolors(condi,:));
            ylim([-0.5 1]);
            % xlim([0 2000]);

        subplot(3,4,4+condi); % xcorr(EEG_L vs EEG_R)  
            time_series1=EEG_L(:,chan)';% plot(time_series1);
            time_series2=EEG_R(:,chan)';% plot(time_series2);
            TimeLength=1000; % min([length(time_series1) length(time_series2)]);
            [r,lags]=xcorr(time_series1(1:TimeLength), time_series2(1:TimeLength), TimeLength-1,'normalized');
            plot(lags./2,r);xlabel('lag of time [ms]');ylabel('Xcorr');
            title('EEG-L & EEG-R','Color',condicolors(condi,:));
            ylim([-1 1]);
            % xlim([-2000 2000]);

        subplot(3,4,4*2+condi); % autocorr(EEG_R)
            lags=1000; % round(length(EEG_R(:,chan)))-1;
            autocorr(EEG_R(:,chan)',lags); xlabel('lags of samples');ylabel('autocorr'); 
            title('EEG-R', 'Color',condicolors(condi,:));
            ylim([-0.5 1]);
            % xlim([0 2000]);
    end
    suptitle(['EEGcorr channel ' labels{chan} ' subject ' num2str(seed)]);
    figureName=['EEGcorr-chan' labels{chan}];
    % save the figure
    saveas(gcf,figureName,'jpg');
    close all;
end

%% PLOT 13:  EEG xcorr all chans (weak anticipation, short-range / local statistics)
% [ autocorr(EEG_L) + xcorr(EEG_L vs EEG_R) + autocorr(EEG_R) ]
% each of the 6 conditions (subplots(3,6,i)) for each of the 32 channels;
% auto save the 32 figures in a folder named - ShortEEG
cd /ssd/zhibin/1overf/20220515_2P/Segmented_data/Plots/
cd /ssd/zhibin/1overf/20220517_2P/Segmented_data/Plots/
cd /ssd/zhibin/1overf/20220518_2P/Segmented_data/Plots/

figure('units','normalized','outerposition',[0 0 1 1]); 
for condi=1:4
    if condi==1; EEG_L=EEGCondi1L; EEG_R=EEGCondi1R; end % R's Error base on L (uncouple)
    if condi==2; EEG_L=EEGCondi2L; EEG_R=EEGCondi2R; end % L's Error base on R (uncouple)
    if condi==3; EEG_L=EEGCondi3L; EEG_R=EEGCondi3R; end % R's Error base on L (unidirectional)
    if condi==4; EEG_L=EEGCondi4L; EEG_R=EEGCondi4R; end % L's Error base on R (unidirectional)

    subplot(3,4,condi); % autocorr(EEG_L) 
        r=[];
        for chan = 1:32;
            lags=1000; % round(length(EEG_L(:,chan)))-1;
            [r(:,chan),lags]=autocorr(EEG_L(:,chan),lags); 
        end
        plotx(lags,r);xlabel('lags of samples');ylabel('autocorr'); yline(0,'r');
        title({['Condi ' conditionNames{condi}],'EEG-L'},'Color',condicolors(condi,:));
        % ylim([-0.5 1]);
        % xlim([0 2000]);
        

    subplot(3,4,4+condi); % xcorr(EEG_L vs EEG_R)  
        r=[];
        for chan = 1:32;
            time_series1=EEG_L(:,chan)';% plot(time_series1);
            time_series2=EEG_R(:,chan)';% plot(time_series2);
            TimeLength=1000; % min([length(time_series1) length(time_series2)]);
            [r(:,chan), lags]=xcorr(time_series1(1:TimeLength), time_series2(1:TimeLength), TimeLength-1,'normalized');
        end
        plotx(lags./2,r);xlabel('lag of time [ms]');ylabel('Xcorr');
        title('EEG-L & EEG-R','Color',condicolors(condi,:));
        % ylim([-0.3 0.3]);
        % xlim([-2000 2000]);

    subplot(3,4,4*2+condi); % autocorr(EEG_R)
        r=[];
        for chan = 1:32;
            lags=1000; % round(length(EEG_R(:,chan)))-1;
            [r(:,chan),lags]= autocorr(EEG_R(:,chan)',lags); 
        end
        plotx(lags,r);xlabel('lags of samples');ylabel('autocorr'); yline(0,'r');
        title('EEG-R', 'Color',condicolors(condi,:));
        % ylim([-0.5 1]);
        % xlim([0 2000]);
end
suptitle(['EEGcorr-all-chans -- ' ' subject ' num2str(seed)]);
legend(labels);

figureName=['EEGcorr-all-chan' num2str(seed)];
% save the figure
saveas(gcf,figureName,'jpg');
saveas(gcf,figureName,'fig'); % save as fig to examine plotx
close all;

%% PLOT 14:  EEG Power spectra Xspectra - single chan (strong anticipation, long-range / global statistics)
% [ Pspectra(EEG_L) + Xspectra(EEG_L vs EEG_R) + Pspectra(EEG_R) ]
% each of the 6 conditions (subplots(3,6,i)) for each of the 32 channels;
% auto save the 32 figures in a folder named - LongEEG
cd /ssd/zhibin/1overf/20220515_2P/Segmented_data/EEGspectra_chans
cd /ssd/zhibin/1overf/20220517_2P/Segmented_data/LongEEG
cd /ssd/zhibin/1overf/20220518_2P/Segmented_data/EEGspectra_chans

maxfreq=50;sr=2000;
for chan = 1:32
    figure('units','normalized','outerposition',[0 0 1 1]);
    for condi=1:4
        if condi==1; EEG_L=EEGCondi1L; EEG_R=EEGCondi1R; end % R's Error base on L (uncouple)
        if condi==2; EEG_L=EEGCondi2L; EEG_R=EEGCondi2R; end % L's Error base on R (uncouple)
        if condi==3; EEG_L=EEGCondi3L; EEG_R=EEGCondi3R; end % R's Error base on L (unidirectional)
        if condi==4; EEG_L=EEGCondi4L; EEG_R=EEGCondi4R; end % L's Error base on R (unidirectional)

        subplot(3,4,condi); % Pspectra(EEG_L) 
            [fcoef1,fcoef2,cprod, freqs] = spectra(EEG_L(:,chan),EEG_L(:,chan),maxfreq,sr);
            plot(freqs,abs(cprod),'Color',condicolors(condi,:));xlabel('freqs');ylabel('power');
            title({['Condi ' conditionNames{condi}],'EEG-L'},'Color',condicolors(condi,:));
    %         ylim([]);

        subplot(3,4,4+condi); % xcorr (EEG_L vs EEG_R) 
            [fcoef1,fcoef2,cprod, freqs] = spectra(EEG_L(:,chan),EEG_R(:,chan),maxfreq,sr);      
            plot(freqs,abs(cprod),'Color',condicolors(condi,:));xlabel('freqs');ylabel('Xspectra');
            title('EEG-L & EEG-R','Color',condicolors(condi,:));
    %         ylim([]);

        subplot(3,4,4*2+condi); %  Pspectra(EEG_R) 
            [fcoef1,fcoef2,cprod, freqs] = spectra(EEG_R(:,chan),EEG_R(:,chan),maxfreq,sr);
            plot(freqs,abs(cprod),'Color',condicolors(condi,:));xlabel('freqs');ylabel('power');
            title('EEG-R','Color',condicolors(condi,:));
    %         ylim([]);
    end
    suptitle(['EEG-spectra ' labels{chan} ' subject ' num2str(seed)]);
    figureName=['EEG-spectra-chan' labels{chan}];
    % save the figure
    saveas(gcf,figureName,'jpg');
    close all;
end

%% PLOT 14-1:  EEG Spectra - single chan - smooth
% [ Pspectra(EEG_L) + Xspectra(EEG_L vs EEG_R) + Pspectra(EEG_R) ]
% each of the 6 conditions (subplots(3,6,i)) for each of the 32 channels;
% auto save the 32 figures in a folder 
cd /ssd/zhibin/1overf/20220518_2P/Segmented_data/EEGspectra_chans_coh_PLOT_14_1
% select number of conditions in the experiment
Totalcondi=4;
Totalcondi=5;
Totalcondi=6;

maxfreq=50;sr=2000;win=10; % win in seconds (for cutting into epochs)
% fselect=[1 2 4 8 16];
% % select freq of [1 2 4 8 16] 
% fselect=[1 2 4 8 16];
% df=1/20; indf=fselect/df;

for chan = 1:32
    figure('units','normalized','outerposition',[0 0 1 1]);
    for condi=1:Totalcondi
        if condi==1; EEG_L=EEGCondi1L; EEG_R=EEGCondi1R; end 
        if condi==2; EEG_L=EEGCondi2L; EEG_R=EEGCondi2R; end 
        if condi==3; EEG_L=EEGCondi3L; EEG_R=EEGCondi3R; end
        if condi==4; EEG_L=EEGCondi4L; EEG_R=EEGCondi4R; end 
        if condi==5; EEG_L=EEGCondi5L; EEG_R=EEGCondi5R; end % R's Error base on L (3Hz bidirection)
        if condi==6; EEG_L=EEGCondi6L; EEG_R=EEGCondi6R; end 
        
        [P1,P2,Px,coh,phase,freqs]  = spectra(EEG_L(:,chan),EEG_R(:,chan),maxfreq,sr,win);
        subplot(2,Totalcondi,condi); % Pspectra(EEG_L) 
            plot(freqs,P1,'r');hold on;plot(freqs,P2,'b');plot(freqs,abs(Px),'color',darkgreen);
            xlabel('freqs');ylabel('Power');legend('L','R','Xspectra');
            title({['Condi ' conditionNames{condi}],'spectra'},'Color',condicolors(condi,:));
            ylim([0 0.02]);xlim([0 50]);

        subplot(2,Totalcondi,Totalcondi+condi); % Xspetra (EEG_L vs EEG_R)   
            plot(freqs,coh,'Color',condicolors(condi,:));xlabel('freqs');ylabel('coh');
            title('coh: EEG-L & EEG-R','Color',condicolors(condi,:));
            ylim([0 0.25]);xlim([0 50]);

    end
    suptitle(['EEG-spectra-coh ' labels{chan} ' subject ' num2str(seed) '-PLOT-14-1']);
    figureName=['EEG-spectra-coh-chan' labels{chan} '-PLOT-14-1'];
    % save the figure
    saveas(gcf,figureName,'fig');
    close all;
end

%% PLOT 14-2:  EEG Spectra - single chan - log log smooth at 10s
% [ Pspectra(EEG_L) + Xspectra(EEG_L vs EEG_R) + Pspectra(EEG_R) ]
% each of the 6 conditions (subplots(3,6,i)) for each of the 32 channels;
% auto save the 32 figures in a folder 
cd /ssd/zhibin/1overf/20220518_2P/Segmented_data/EEGspectra_chans_coh_PLOT_14_2
% select number of conditions in the experiment
Totalcondi=4;
% Totalcondi=5;
% Totalcondi=6;

maxfreq=50;sr=2000;win=10; % win in seconds (for cutting into epochs)

for chan = 1:32
    figure('units','normalized','outerposition',[0 0 1 1]);
    for condi=1:Totalcondi
        if condi==1; EEG_L=EEGCondi1L; EEG_R=EEGCondi1R; end 
        if condi==2; EEG_L=EEGCondi2L; EEG_R=EEGCondi2R; end 
        if condi==3; EEG_L=EEGCondi3L; EEG_R=EEGCondi3R; end
        if condi==4; EEG_L=EEGCondi4L; EEG_R=EEGCondi4R; end 
        if condi==5; EEG_L=EEGCondi5L; EEG_R=EEGCondi5R; end % R's Error base on L (3Hz bidirection)
        if condi==6; EEG_L=EEGCondi6L; EEG_R=EEGCondi6R; end 
        
        % calculate for each condition
        [P1,P2,Px,coh,phase,freqs]  = spectra(EEG_L(:,chan),EEG_R(:,chan),maxfreq,sr,win);
        
        subplot(1,4,1); % P1
            hold on;
            loglog(smoothing(freqs,win),smoothing(P1,win),'Color',condicolors(condi,:));
            if condi==Totalcondi;
                legend(conditionNames(1:Totalcondi));
                xlabel('log-freqs');ylabel('log-Power');
                ylim([0 0.005]);xlim([1.5 50]);
                title('Spectra EEG-L')
                hold off;
            end
            
        subplot(1,4,2); % P2
            hold on;
            loglog(smoothing(freqs,win),smoothing(P2,win),'Color',condicolors(condi,:));
            if condi==Totalcondi;
                legend(conditionNames(1:Totalcondi));
                xlabel('log-freqs');ylabel('log-Power');
                ylim([0 0.005]);xlim([1.5 50]);
                title('Spectra EEG-R')
                hold off;
            end
            
        subplot(1,4,3); % Px
            hold on;
            loglog(smoothing(freqs,win),smoothing(abs(Px),win),'Color',condicolors(condi,:));
            if condi==Totalcondi;
                legend(conditionNames(1:Totalcondi));
                xlabel('log-freqs');ylabel('log-Power');
                ylim([0 0.001]);xlim([1.5 50]);
                title('XSpectra')
                hold off;
            end
            
        subplot(1,4,4); % Coh
            hold on;
            % plot(freqs,coh,'Color',condicolors(condi,:));imagesc(coh);colorbar;
            plot(freqs,coh,'Color',condicolors(condi,:));
            if condi==Totalcondi;
                legend(conditionNames(1:Totalcondi));
                xlabel('freqs');ylabel('coh');
                ylim([0 0.25]);xlim([1.5 50]);
                title('Coh')
                hold off;
            end

    end
    suptitle(['EEG-loglogspectra-coh-10s' labels{chan} ' subject ' num2str(seed) '-PLOT-14-2']);
    figureName=['EEG-loglogspectra-coh-10s-chan' labels{chan} '-PLOT-14-2'];
    % save the figure
    saveas(gcf,figureName,'fig');
    close all;
end

%% PLOT 14-3:  EEG Spectra - single chan - log log smooth at 10s plot each level with 4 conditions
% [ Pspectra(EEG_L) + Xspectra(EEG_L vs EEG_R) + Pspectra(EEG_R) ]
% each of the 6 conditions (subplots(3,6,i)) for each of the 32 channels;
% auto save the 32 figures in a folder 
cd /ssd/zhibin/1overf/20220518_2P/Segmented_data/EEGspectra_chans_coh_PLOT_14_3
% select number of conditions in the experiment
Totalcondi=4;
% Totalcondi=5;
% Totalcondi=6;

maxfreq=50;sr=2000;win=10; % win in seconds (for cutting into epochs)
HML_labels={'High corr','Medium corr','Low corr'};

for chan = 1:32
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
            
            subplot(3,4,4*(i-1)+1); % P1
                hold on;
                loglog(smoothing(freqs,win),smoothing(P1,win),'Color',condicolors(condi,:));
                if condi==Totalcondi;
                    legend(conditionNames(1:Totalcondi));
                    xlabel('log-freqs');ylabel('log-Power');
                    ylim([0 0.01]);xlim([1.5 50]);
                    title({['Spectra EEG-L -' HML_labels{i}]});
                    hold off;
                end

            subplot(3,4,4*(i-1)+2); % P2
                hold on;
                loglog(smoothing(freqs,win),smoothing(P2,win),'Color',condicolors(condi,:));
                if condi==Totalcondi;
                    legend(conditionNames(1:Totalcondi));
                    xlabel('log-freqs');ylabel('log-Power');
                    ylim([0 0.01]);xlim([1.5 50]);
                    title({['Spectra EEG-R -' HML_labels{i}]});
                    hold off;
                end

            subplot(3,4,4*(i-1)+3); % Px
                hold on;
                loglog(smoothing(freqs,win),smoothing(abs(Px),win),'Color',condicolors(condi,:));
                if condi==Totalcondi;
                    legend(conditionNames(1:Totalcondi));
                    xlabel('log-freqs');ylabel('log-Power');
                    ylim([0 0.002]);xlim([1.5 50]);
                    title({['XSpectra -' HML_labels{i}]});
                    hold off;
                end

            subplot(3,4,4*(i-1)+4); % Coh
                hold on;
                % plot(freqs,coh,'Color',condicolors(condi,:));imagesc(coh);colorbar;
                plot(freqs,coh,'Color',condicolors(condi,:));
                if condi==Totalcondi;
                    legend(conditionNames(1:Totalcondi),'location','northeast');
                    xlabel('freqs');ylabel('coh');
                    ylim([0 1]);xlim([1.5 50]);
                    title({['Coh -' HML_labels{i}]});
                    hold off;
                end
        end

    end
    suptitle(['EEG-loglogspectra-coh-10s' labels{chan} ' subject ' num2str(seed) '-PLOT-14-3']);
    figureName=['EEG-loglogspectra-coh-10s-chan' labels{chan} '-PLOT-14-3'];
    % save the figure
    saveas(gcf,figureName,'fig');
    close all;
end

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

for chan = 1:32
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
                    ylim([0 0.01]);xlim([1.5 50]);
                    title({['Condi ' conditionNames{condi}],['Spectra EEG-R']},'Color',condicolors(condi,:));
                    hold off;
                end

            subplot(4,4,4+condi); % P2
                hold on;
                loglog(smoothing(freqs,win),smoothing(P2,win),'Color',HNLcolors(i,:));
                if i==3;
                    legend(HML_labels(1:3));
                    xlabel('log-freqs');ylabel('log-Power');
                    ylim([0 0.01]);xlim([1.5 50]);
                    title({['Spectra EEG-R']},'Color',condicolors(condi,:));
                    hold off;
                end

            subplot(4,4,4*2+condi); % Px
                hold on;
                loglog(smoothing(freqs,win),smoothing(abs(Px),win),'Color',HNLcolors(i,:));
                if i==3;
                    legend(HML_labels(1:3));
                    xlabel('log-freqs');ylabel('log-Power');
                    ylim([0 0.002]);xlim([1.5 50]);
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
                    ylim([0 1]);xlim([1.5 50]);
                    title({['Coh']},'Color',condicolors(condi,:));
                    hold off;
                end
        end

    end
    suptitle(['EEG-loglogspectra-coh-10s' labels{chan} ' subject ' num2str(seed) '-PLOT-14-4']);
    figureName=['EEG-loglogspectra-coh-10s-chan' labels{chan} '-PLOT-14-4'];
    % save the figure
    saveas(gcf,figureName,'fig');
    close all;
end


%% PLOT 15:  Xcorr EMG vs BP/FB(boolean) (short-range):(equal length) 
% [ xcorr(EMG_L vs BP_L) + xcorr(EMG_L vs FB_L) + xcorr(EMG_R vs BP_R) + xcorr(EMG_R vs FB_R) ]
% 6 conditions -subplots(4,6,condi) 
% select number of conditions in the experiment
Totalcondi=4;
Totalcondi=5;

figure('units','normalized','outerposition',[0 0 1 1]);
for condi=1:Totalcondi
    if condi==1; EMG_L=EMGCondi1L; BP_L=BPCondi1L; FB_L=FBCondi1L; EMG_R=EMGCondi1R; BP_R=BPCondi1R; FB_R=FBCondi1R; end
    if condi==2; EMG_L=EMGCondi2L; BP_L=BPCondi2L; FB_L=FBCondi2L; EMG_R=EMGCondi2R; BP_R=BPCondi2R; FB_R=FBCondi2R; end
    if condi==3; EMG_L=EMGCondi3L; BP_L=BPCondi3L; FB_L=FBCondi3L; EMG_R=EMGCondi3R; BP_R=BPCondi3R; FB_R=FBCondi3R; end
    if condi==4; EMG_L=EMGCondi4L; BP_L=BPCondi4L; FB_L=FBCondi4L; EMG_R=EMGCondi4R; BP_R=BPCondi4R; FB_R=FBCondi4R; end
    if condi==5; EMG_L=EMGCondi5L; BP_L=BPCondi5L; FB_L=FBCondi5L; EMG_R=EMGCondi5R; BP_R=BPCondi5R; FB_R=FBCondi5R; end
    if condi==6; EMG_L=EMGCondi6L; BP_L=BPCondi6L; FB_L=FBCondi6L; EMG_R=EMGCondi6R; BP_R=BPCondi6R; FB_R=FBCondi6R; end
    % compute
    subplot(4,Totalcondi,condi); %  xcorr(EMG_L vs BP_L)
        [r,lags]=xcorr(EMG_L', BP_L', 1000,'normalized');
        plot(lags./2,r);xlabel('time [ms]');ylabel('Xcorr');xline(0,'Color',deepyellow);
        title({['Condi ' conditionNames{condi}],'EMG-L & BP-L'},'Color',condicolors(condi,:));
        % xlim([-2e3 2e3]);ylim([-0.005 0.1]);
    subplot(4,Totalcondi,Totalcondi+condi); %  xcorr(EMG_L vs FB_L)
        [r,lags]=xcorr(EMG_L', FB_L', 1000,'normalized');
        plot(lags./2,r);xlabel('time [ms]');ylabel('Xcorr');xline(0,'Color',deepyellow);
        title('EMG-L & FB-L','Color',condicolors(condi,:));
        % xlim([-2e3 2e3]);ylim([-0.005 0.1]);
    subplot(4,Totalcondi,Totalcondi*2+condi); %   xcorr(EMG_R vs BP_R)
        [r,lags]=xcorr(EMG_R', BP_R', 1000,'normalized');
        plot(lags./2,r);xlabel('time [ms]');ylabel('Xcorr');xline(0,'Color',deepyellow);
        title('EMG-R & BP-R','Color',condicolors(condi,:));
        % xlim([-2e3 2e3]);ylim([-0.005 0.1]);
    subplot(4,Totalcondi,Totalcondi*3+condi); %  xcorr(EMG_R vs FB_R)
        [r,lags]=xcorr(EMG_R', FB_R', 1000, 'normalized');
        plot(lags./2,r);xlabel('time [ms]');ylabel('Xcorr');xline(0,'Color',deepyellow);
        title('EMG-R & FB-R','Color',condicolors(condi,:));
        % xlim([-2e3 2e3]);ylim([-0.005 0.1]);
end
suptitle(['xcorr: EMG-BP/FB ' ' -- subject ' num2str(seed)]);
        
figureName=['xcorrEMG-BPFBboolean' num2str(seed)];
    % save the figure
    saveas(gcf,figureName,'jpg');

%% PLOT 16:  Xcorr EEG vs BP/FB(boolean) - single chan (short-range): (equal length) (low correlation)
% [ xcorr(EMG_L vs BP_L) + xcorr(EMG_L vs FB_L) + xcorr(EMG_R vs BP_R) + xcorr(EMG_R vs FB_R) ]
% 6 conditions -subplots(4,6,condi) 
cd /ssd/zhibin/1overf/20220518_2P/Segmented_data/XcorrEEG_BP_FB
% select number of conditions in the experiment
Totalcondi=4;
Totalcondi=5;

for chan=1:32
    figure('units','normalized','outerposition',[0 0 1 1]);
    for condi=1:Totalcondi
        if condi==1; EEG_L=EEGCondi1L; BP_L=BPCondi1L; FB_L=FBCondi1L; EEG_R=EEGCondi1R; BP_R=BPCondi1R; FB_R=FBCondi1R; end
        if condi==2; EEG_L=EEGCondi2L; BP_L=BPCondi2L; FB_L=FBCondi2L; EEG_R=EEGCondi2R; BP_R=BPCondi2R; FB_R=FBCondi2R; end
        if condi==3; EEG_L=EEGCondi3L; BP_L=BPCondi3L; FB_L=FBCondi3L; EEG_R=EEGCondi3R; BP_R=BPCondi3R; FB_R=FBCondi3R; end
        if condi==4; EEG_L=EEGCondi4L; BP_L=BPCondi4L; FB_L=FBCondi4L; EEG_R=EEGCondi4R; BP_R=BPCondi4R; FB_R=FBCondi4R; end
        if condi==5; EEG_L=EEGCondi5L; BP_L=BPCondi5L; FB_L=FBCondi5L; EEG_R=EEGCondi5R; BP_R=BPCondi5R; FB_R=FBCondi5R; end
        if condi==6; EEG_L=EEGCondi6L; BP_L=BPCondi6L; FB_L=FBCondi6L; EEG_R=EEGCondi6R; BP_R=BPCondi6R; FB_R=FBCondi6R; end
        % compute
        subplot(4,Totalcondi,condi); %  xcorr(EEG_L vs BP_L)
            [r,lags]=xcorr(EEG_L(:,chan)', BP_L', 1000,'normalized');
            plot(lags./2,r);xlabel('time [ms]');ylabel('Xcorr');xline(0,'Color',deepyellow);
            title({['Condi ' conditionNames{condi}],'EEG-L & BP-L'},'Color',condicolors(condi,:));
            % xlim([-2e3 2e3]);ylim([-0.005 0.1]);
        subplot(4,Totalcondi,Totalcondi+condi); %  xcorr(EEG_L vs FB_L)
            [r,lags]=xcorr(EEG_L(:,chan)', FB_L', 1000,'normalized');
            plot(lags./2,r);xlabel('time [ms]');ylabel('Xcorr');xline(0,'Color',deepyellow);
            title('EEG-L & FB-L', 'Color',condicolors(condi,:));
            % xlim([-2e3 2e3]);ylim([-0.005 0.1]);
        subplot(4,Totalcondi,Totalcondi*2+condi); %   xcorr(EEG_R vs BP_R)
            [r,lags]=xcorr(EEG_R(:,chan)', BP_R', 1000,'normalized');
            plot(lags./2,r);xlabel('time [ms]');ylabel('Xcorr');xline(0,'Color',deepyellow);
            title('EEG-R & BP-R','Color',condicolors(condi,:));
            % xlim([-2e3 2e3]);ylim([-0.005 0.1]);
        subplot(4,Totalcondi,Totalcondi*3+condi); %  xcorr(EEG_R vs FB_R)
            [r,lags]=xcorr(EEG_R(:,chan)', FB_R', 1000,'normalized');
            plot(lags./2,r);xlabel('time [ms]');ylabel('Xcorr');xline(0,'Color',deepyellow);
            title('EEG-R & FB-R','Color',condicolors(condi,:));
            % xlim([-2e3 2e3]);ylim([-0.005 0.1]);
    end
    suptitle(['xcorr: EEG-BP/FB EEG channel ' labels{chan} '   subject ' num2str(seed)]);
        
figureName=['xcorrEEG-BPFBboolean_EEG_chan' labels{chan} '--subject--' num2str(seed)];
    % save the figure
    saveas(gcf,figureName,'jpg');
end

close all;

%% PLOT 17:  Xcorr EEG vs BP/FB(boolean) - all chans (short-range): (equal length) (low correlation)
% [ xcorr(EMG_L vs BP_L) + xcorr(EMG_L vs FB_L) + xcorr(EMG_R vs BP_R) + xcorr(EMG_R vs FB_R) ]
% 6 conditions -subplots(4,6,condi) 
cd /ssd/zhibin/1overf/20220518_2P/Segmented_data/Plots

figure('units','normalized','outerposition',[0 0 1 1]);
for condi=1:4
    if condi==1; EEG_L=EEGCondi1L; BP_L=BPCondi1L; FB_L=FBCondi1L; EEG_R=EEGCondi1R; BP_R=BPCondi1R; FB_R=FBCondi1R; end
    if condi==2; EEG_L=EEGCondi2L; BP_L=BPCondi2L; FB_L=FBCondi2L; EEG_R=EEGCondi2R; BP_R=BPCondi2R; FB_R=FBCondi2R; end
    if condi==3; EEG_L=EEGCondi3L; BP_L=BPCondi3L; FB_L=FBCondi3L; EEG_R=EEGCondi3R; BP_R=BPCondi3R; FB_R=FBCondi3R; end
    if condi==4; EEG_L=EEGCondi4L; BP_L=BPCondi4L; FB_L=FBCondi4L; EEG_R=EEGCondi4R; BP_R=BPCondi4R; FB_R=FBCondi4R; end
    % compute
    subplot(4,4,condi); %  xcorr(EEG_L vs BP_L)
        r=[];
        for chan=1:32
            [r(:,chan),lags]=xcorr(EEG_L(:,chan), BP_L, 1000,'normalized');
        end
        plotx(lags./2,r);xlabel('time [ms]');ylabel('Xcorr');xline(0,'Color',deepyellow);
        title({['Condi ' conditionNames{condi}],'EEG-L & BP-L'},'Color',condicolors(condi,:));
        % xlim([-2e3 2e3]);ylim([-0.005 0.1]);
    subplot(4,4,4+condi); %  xcorr(EEG_L vs FB_L)
        r=[];
        for chan=1:32
            [r(:,chan),lags]=xcorr(EEG_L(:,chan), FB_L, 1000,'normalized');
        end
        plotx(lags./2,r);xlabel('time [ms]');ylabel('Xcorr');xline(0,'Color',deepyellow);
        title('EEG-L & FB-L', 'Color',condicolors(condi,:));
        % xlim([-2e3 2e3]);ylim([-0.005 0.1]);
    subplot(4,4,4*2+condi); %   xcorr(EEG_R vs BP_R)
        r=[];
        for chan=1:32
            [r(:,chan),lags]=xcorr(EEG_R(:,chan), BP_R, 1000,'normalized');
        end
        plotx(lags./2,r);xlabel('time [ms]');ylabel('Xcorr');xline(0,'Color',deepyellow);
        title('EEG-R & BP-R','Color',condicolors(condi,:));
        % xlim([-2e3 2e3]);ylim([-0.005 0.1]);
    subplot(4,4,4*3+condi); %  xcorr(EEG_R vs FB_R)
        r=[];
        for chan=1:32
            [r(:,chan),lags]=xcorr(EEG_R(:,chan), FB_R, 1000,'normalized');
        end
        plotx(lags./2,r);xlabel('time [ms]');ylabel('Xcorr');xline(0,'Color',deepyellow);
        title('EEG-R & FB-R','Color',condicolors(condi,:));
        % xlim([-2e3 2e3]);ylim([-0.005 0.1]);
end
suptitle(['xcorr: EEG-BP/FB EEG all channels -- subject ' num2str(seed)]);
legend(labels);

figureName=['xcorrEEG-BPFBboolean_EEG_all-chans--subject--' num2str(seed)];
% save the figure
saveas(gcf,figureName,'jpg');
saveas(gcf,figureName,'fig');
close all;

%% PLOT 18:  Xcorr EEG vs EMG (equal length) - single chan (short-range): (low correlation) 
% [ xcorr(EEG_L vs EMG_L) + xcorr(EEG_L vs EMG_R) xcorr(EEG_R vs EMG_L) xcorr(EEG_R vs EMG_R) ]
% 6 conditions -subplots(2,6,condi) 
cd /ssd/zhibin/1overf/20220518_2P/Segmented_data/XcorrEEG_EMG
% select number of conditions in the experiment
Totalcondi=4;
Totalcondi=5;

for chan=1:32
    figure('units','normalized','outerposition',[0 0 1 1]);
    for condi=1:Totalcondi
        if condi==1; EEG_L=EEGCondi1L; EMG_L=EMGCondi1L; EEG_R=EEGCondi1R; EMG_R=EMGCondi1R; end
        if condi==2; EEG_L=EEGCondi2L; EMG_L=EMGCondi2L; EEG_R=EEGCondi2R; EMG_R=EMGCondi2R; end
        if condi==3; EEG_L=EEGCondi3L; EMG_L=EMGCondi3L; EEG_R=EEGCondi3R; EMG_R=EMGCondi3R; end
        if condi==4; EEG_L=EEGCondi4L; EMG_L=EMGCondi4L; EEG_R=EEGCondi4R; EMG_R=EMGCondi4R; end
        if condi==5; EEG_L=EEGCondi5L; EMG_L=EMGCondi5L; EEG_R=EEGCondi5R; EMG_R=EMGCondi5R; end
        % compute
        subplot(4,Totalcondi,condi); %  xcorr(EEG_L vs EMG_L)
            [r,lags]=xcorr(EEG_L(:,chan)', EMG_L', 1000,'normalized');
            plot(lags./2,r);xlabel('time [ms]');ylabel('Xcorr');xline(0,'Color',deepyellow);
            title({['Condi ' conditionNames{condi}], 'EEG-L & EMG-L'},'Color',condicolors(condi,:));
            % xlim([-2e3 2e3]);ylim([-0.1 0.1]);
        subplot(4,Totalcondi,Totalcondi+condi); %  xcorr(EEG_L vs EMG_R)
            minLength=1000; % min([size(EEG_L,1) length(EMG_R)]);
            [r,lags]=xcorr(EEG_L(1:minLength,chan)', EMG_R(1:minLength)', minLength,'normalized');
            plot(lags./2,r);xlabel('time [ms]');ylabel('Xcorr');xline(0,'Color',deepyellow);
            title('EEG-L & EMG-R','Color',condicolors(condi,:));
            % xlim([-2e3 2e3]);ylim([-0.1 0.1]);
        subplot(4,Totalcondi,Totalcondi*2+condi); %   xcorr(EEG_R vs EMG_L)
            minLength=1000; % min([size(EEG_L,1) length(EMG_R)]);
            [r,lags]=xcorr(EEG_R(1:minLength,chan)', EMG_L(1:minLength)', minLength,'normalized');
            plot(lags./2,r);xlabel('time [ms]');ylabel('Xcorr');xline(0,'Color',deepyellow);
            title('EEG-R & EMG-L','Color',condicolors(condi,:));
            % xlim([-2e3 2e3]);ylim([-0.1 0.1]);
        subplot(4,Totalcondi,Totalcondi*3+condi); %  xcorr(EEG_R vs FB_R)
            [r,lags]=xcorr(EEG_R(:,chan)', EMG_R', 1000,'normalized');
            plot(lags./2,r);xlabel('time [ms]');ylabel('Xcorr');xline(0,'Color',deepyellow);
            title('EEG-R & EMG-R', 'Color',condicolors(condi,:));
            % xlim([-2e3 2e3]);ylim([-0.1 0.1]);
    end
    suptitle(['xcorr-EEG-EMG  EEG channel ' labels{chan} '   subject ' num2str(seed)]);
        
figureName=['xcorr-EEG-EMG--EEG--chan' labels{chan} '--subject--' num2str(seed)];
    % save the figure
    saveas(gcf,figureName,'jpg');
end

close all;

0%% Photocell
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing

figure('units','normalized','outerposition',[0 0 1 1]);
win=20;
for condi=1:4
    if condi==1; BP_L=BPCondi1L; FB_L=FBCondi1L; locs_L=FBlocsCondi1L; BP_R=BPCondi1R; FB_R=FBCondi1R; locs_R=FBlocsCondi1R; end
    if condi==2; BP_L=BPCondi2L; FB_L=FBCondi2L; locs_L=FBlocsCondi2L; BP_R=BPCondi2R; FB_R=FBCondi2R; locs_R=FBlocsCondi2R; end
    if condi==3; BP_L=BPCondi3L; FB_L=FBCondi3L; locs_L=FBlocsCondi3L; BP_R=BPCondi3R; FB_R=FBCondi3R; locs_R=FBlocsCondi3R; end
    if condi==4; BP_L=BPCondi4L; FB_L=FBCondi4L; locs_L=FBlocsCondi4L; BP_R=BPCondi4R; FB_R=FBCondi4R; locs_R=FBlocsCondi4R; end
    % compute
    subplot(6,4,condi); %  hist(BP_L)
        histogram(Calinterval(BP_L')./2,600,'FaceColor',condicolors(condi,:),'EdgeColor',condicolors(condi,:));
        xlabel('BP-intervals (ms)');ylabel('number of BP intervals');xlim([0 1200]);
        legend(['total ' num2str(length(Calinterval(BP_L')))]);
        title({['Condi ' conditionNames{condi}],' BP-R interval'},'Color',condicolors(condi,:));
    subplot(6,4,4+condi); %  hist(FB_L)
        histogram(Calinterval(FB_L')./2,600,'FaceColor',condicolors(condi,:),'EdgeColor',condicolors(condi,:));
        xlabel('FB-intervals (ms)');ylabel('number of FB intervals');xlim([0 1200]);
        legend(['total ' num2str(length(Calinterval(FB_L')))]);
        title('FB-L','Color',condicolors(condi,:));
    subplot(6,4,4*2+condi); %  hist(FB_L)
        histogram(Calinterval(locs_L)./2,600,'FaceColor',condicolors(condi,:),'EdgeColor',condicolors(condi,:));
        xlabel('locs-intervals (ms)');ylabel('number of FB intervals');xlim([0 1200]);
        legend(['total ' num2str(length(Calinterval(locs_L)))]);
        title('locs-L','Color',condicolors(condi,:));
    subplot(6,4,4*3+condi); %  hist(FB_L)
        histogram(Calinterval(locs_R)./2,600,'FaceColor',condicolors(condi,:),'EdgeColor',condicolors(condi,:));
        xlabel('locs-intervals (ms)');ylabel('number of FB intervals');xlim([0 1200]);
        legend(['total ' num2str(length(Calinterval(locs_R)))]);
        title('locs-L','Color',condicolors(condi,:));
    subplot(6,4,4*4+condi); %   hist(FP_R)
        histogram(Calinterval(FB_R')./2,600,'FaceColor',condicolors(condi,:),'EdgeColor',condicolors(condi,:));
        xlabel('FB-intervals (ms)');ylabel('number of FB intervals');xlim([0 1200]);
        legend(['total ' num2str(length(Calinterval(FB_R')))]);
        title('FB-R','Color',condicolors(condi,:));
    subplot(6,4,4*5+condi); %  hist(BP_R)
        histogram(Calinterval(BP_R')./2,600,'FaceColor',condicolors(condi,:),'EdgeColor',condicolors(condi,:));
        xlabel('BP-intervals (ms)');ylabel('number of BP intervals');xlim([0 1200]);
        legend(['total ' num2str(length(Calinterval(BP_R')))]);
        title('BP-R','Color',condicolors(condi,:));
end
suptitle(['BP-FB-locs-intervals ' ' -- subject ' num2str(seed)]);

figureName=['BP-FB-locs-hist -- ' num2str(seed)];
% save the figure
saveas(gcf,figureName,'jpg');

figure('units','normalized','outerposition',[0 0 1 1]);
win=20;
for condi=1:4
    if condi==1; BP_L=BPCondi1L; FB_L=FBCondi1L; locs_L=FBlocsCondi1L; BP_R=BPCondi1R; FB_R=FBCondi1R; locs_R=FBlocsCondi1R; end
    if condi==2; BP_L=BPCondi2L; FB_L=FBCondi2L; locs_L=FBlocsCondi2L; BP_R=BPCondi2R; FB_R=FBCondi2R; locs_R=FBlocsCondi2R; end
    if condi==3; BP_L=BPCondi3L; FB_L=FBCondi3L; locs_L=FBlocsCondi3L; BP_R=BPCondi3R; FB_R=FBCondi3R; locs_R=FBlocsCondi3R; end
    if condi==4; BP_L=BPCondi4L; FB_L=FBCondi4L; locs_L=FBlocsCondi4L; BP_R=BPCondi4R; FB_R=FBCondi4R; locs_R=FBlocsCondi4R; end
    % compute
    subplot(4,4,condi); %  (BP_L + FB_L)
        BP_inter1=Calinterval(locs_L)./2; BP_inter2=Calinterval(FB_L')./2;
        BPfrac_time1=linspace(1,5,length(BP_inter1)); BPfrac_time2=linspace(1,5,length(BP_inter2)); 
        plot(BPfrac_time1,BP_inter1,'b');hold on; plot(BPfrac_time2,BP_inter2,'r');
        xlabel('fractional time (min)');ylabel('interval (ms)');
        legend('locs-L','FB-L');
        title('L recording','Color',condicolors(condi,:));
    subplot(4,4,4+condi); %  (FP_R + BP_R)
        BP_inter1=Calinterval(locs_R)./2; BP_inter2=Calinterval(FB_R')./2;
        BPfrac_time1=linspace(1,5,length(BP_inter1)); BPfrac_time2=linspace(1,5,length(BP_inter2)); 
        plot(BPfrac_time1,BP_inter1,'b');hold on; plot(BPfrac_time2,BP_inter2,'r');
        xlabel('fractional time (min)');ylabel('interval (ms)');
        legend('locs-R','FB-R');
        title('R recording','Color',condicolors(condi,:));
    subplot(4,4,4*2+condi); %  
        BP_inter1=smoothing(Calinterval(locs_L)./2,win); BP_inter2=smoothing(Calinterval(FB_L')./2,win);
        BPfrac_time1=linspace(1,5,length(BP_inter1)); BPfrac_time2=linspace(1,5,length(BP_inter2)); 
        plot(BPfrac_time1,BP_inter1,'b');hold on; plot(BPfrac_time2,BP_inter2,'r');
        xlabel('fractional time (min)');ylabel('interval (ms)');
        legend('locs-L','FB-L');
        title('L recording smoothed','Color',condicolors(condi,:));
    subplot(4,4,4*3+condi); % 
        BP_inter1=smoothing(Calinterval(locs_R)./2,win); BP_inter2=smoothing(Calinterval(FB_R')./2,win);
        BPfrac_time1=linspace(1,5,length(BP_inter1)); BPfrac_time2=linspace(1,5,length(BP_inter2)); 
        plot(BPfrac_time1,BP_inter1,'b');hold on; plot(BPfrac_time2,BP_inter2,'r');
        xlabel('fractional time (min)');ylabel('interval (ms)');
        legend('locs-R','FB-R');
        title('R recording smoothed','Color',condicolors(condi,:));
end
suptitle(['Fractional button press intervals smoothed win ' num2str(win) ' -- subject ' num2str(seed)]);

