%% Plots - syncopation - color scheme
red   = [1 0 0];
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

% % test color
% showcolor=deepyellow;
% imagesc(cat(3,showcolor(1),showcolor(2),showcolor(3)));
%%  conditionNames
conditionNames;

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
suptitle(['Button press intervals ' ' -- subject ' num2str(seed)]);

% cd /ssd/zhibin/1overf/20220518_2P/Segmented_data/Plots
figureName=['BP-hist -- ' num2str(seed)];
% save the figure
saveas(gcf,figureName,'jpg');

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
    subplot(4,4,4+condi); %  (FP_R + BP_R)
        BP_inter1=Calinterval(BP_R')./2; BP_inter2=Calinterval(FB_R')./2;
        BPfrac_time1=linspace(1,5,length(BP_inter1)); BPfrac_time2=linspace(1,5,length(BP_inter2)); 
        plot(BPfrac_time1,BP_inter1,'b');hold on; plot(BPfrac_time2,BP_inter2,'r');
        xlabel('fractional time (min)');ylabel('interval (ms)');
        legend('BP-R','FB-R');
        title('R recording','Color',condicolors(condi,:));
    subplot(4,4,4*2+condi); %   smooth[(BP_L + FB_L)]
        BP_inter1=smoothing(Calinterval(BP_L')./2,win); BP_inter2=smoothing(Calinterval(FB_L')./2,win);
        BPfrac_time1=linspace(1,5,length(BP_inter1)); BPfrac_time2=linspace(1,5,length(BP_inter2)); 
        plot(BPfrac_time1,BP_inter1,'r');hold on; plot(BPfrac_time2,BP_inter2,'b');
        xlabel('fractional time (min)');ylabel('interval (ms)');
        legend('BP-L','FB-L');
        title('L recording smoothed','Color',condicolors(condi,:));
    subplot(4,4,4*3+condi); % smooth[(FP_R + BP_R)]
        BP_inter1=smoothing(Calinterval(BP_R')./2,win); BP_inter2=smoothing(Calinterval(FB_R')./2,win);
        BPfrac_time1=linspace(1,5,length(BP_inter1)); BPfrac_time2=linspace(1,5,length(BP_inter2)); 
        plot(BPfrac_time1,BP_inter1,'b');hold on; plot(BPfrac_time2,BP_inter2,'r');
        xlabel('fractional time (min)');ylabel('interval (ms)');
        legend('BP-R','FB-R');
        title('R recording smoothed','Color',condicolors(condi,:));
end
suptitle(['Fractional button press intervals smoothed win ' num2str(win) ' -- subject ' num2str(seed)]);

% cd /ssd/zhibin/1overf/20220518_2P/Segmented_data/Plots
figureName=['BP-FB-frac-smooth-win-' num2str(win) '----' num2str(seed)];
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
    subplot(4,4,4+condi); 
        win=20;
        BP_inter1=smoothing(Calinterval(BP_L')./2,win); BP_inter2=smoothing(Calinterval(BP_R')./2,win);
        BPfrac_time1=linspace(1,5,length(BP_inter1)); BPfrac_time2=linspace(1,5,length(BP_inter2)); 
        plot(BPfrac_time1,BP_inter1,'r');hold on; plot(BPfrac_time2,BP_inter2,'b');
        xlabel('fractional time (min)');ylabel('interval (ms)');
        legend('BP-L','BP-R');
        title(['BP smooth win ' num2str(win)],'Color',condicolors(condi,:));
    subplot(4,4,4*2+condi); 
        win=50;
        BP_inter1=smoothing(Calinterval(BP_L')./2,win); BP_inter2=smoothing(Calinterval(BP_R')./2,win);
        BPfrac_time1=linspace(1,5,length(BP_inter1)); BPfrac_time2=linspace(1,5,length(BP_inter2)); 
        plot(BPfrac_time1,BP_inter1,'r');hold on; plot(BPfrac_time2,BP_inter2,'b');
        xlabel('fractional time (min)');ylabel('interval (ms)');
        legend('BP-L','BP-R');
        title(['BP smooth win ' num2str(win)],'Color',condicolors(condi,:));
    subplot(4,4,4*3+condi); 
        win=100;
        BP_inter1=smoothing(Calinterval(BP_L')./2,win); BP_inter2=smoothing(Calinterval(BP_R')./2,win);
        BPfrac_time1=linspace(1,5,length(BP_inter1)); BPfrac_time2=linspace(1,5,length(BP_inter2)); 
        plot(BPfrac_time1,BP_inter1,'r');hold on; plot(BPfrac_time2,BP_inter2,'b');
        xlabel('fractional time (min)');ylabel('interval (ms)');
        legend('BP-L','BP-R');
        title(['BP smooth win ' num2str(win)],'Color',condicolors(condi,:));
end
suptitle(['Fractional BP intervals smoothed -- subject ' num2str(seed)]);

% cd /ssd/zhibin/1overf/20220518_2P/Segmented_data/Plots
figureName=['BP-frac-smooth----' num2str(seed)];
% save the figure
saveas(gcf,figureName,'jpg');


%% PLOT 2:  BP interval corr (short-range / local statistics; weak anticipation) on interval
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
        [r,lags]=xcorr(time_series1(1:TimeLength)', time_series2(1:TimeLength)', 430,'coeff');
        plot(lags,r);xlabel('lags');ylabel('Xcorr');xline(0,'--','Color',deepyellow);
        title('xcorr BP-L & FB-L','Color',condicolors(condi,:));
        ylim([-0.2 0.5]);
    subplot(4,4,4*2+condi); %   xcorr(BP_R vs FB_R) 
        time_series1=Calinterval(BP_R');% plot(time_series1);
        time_series2=Calinterval(FB_R');% plot(time_series2);
        TimeLength=min([length(time_series1) length(time_series2)]);
        time_series1=time_series1-mean(time_series1);% plot(time_series1);
        time_series2=time_series2-mean(time_series2);% plot(time_series2);
        [r,lags]=xcorr(time_series1(1:TimeLength), time_series2(1:TimeLength), 430,'normalized');
        plot(lags,r);xlabel('lags');ylabel('Xcorr');xline(0,'--','Color',deepyellow);
        title('xcorr BP-R & FB-R','Color',condicolors(condi,:));
        ylim([-0.2 0.5]);
    subplot(4,4,4*3+condi); %  autocorr(BP_R)
        autocorr(Calinterval(BP_R'),430); xlabel('lags');ylabel('autocorr'); 
        title('autocorr BP-R','Color',condicolors(condi,:));
end
suptitle(['auto/x corr tapping intervals' ' -- subject ' num2str(seed)]);

% cd /ssd/zhibin/1overf/20220518_2P/Segmented_data/Plots    
figureName=['BPcorr -- ' num2str(seed)];
    % save the figure
    saveas(gcf,figureName,'jpg');


%% PLOT 3:  BP interval corr smooth (better!) (short-range / local statistics; weak anticipation) on interval
% [ autocorr(BP_L) + xcorr(BP_L vs FB_L) + xcorr(BP_R vs FB_R) + autocorr(BP_R) ]
% 6 conditions -subplots(4,6,condi) 

figure('units','normalized','outerposition',[0 0 1 1]);
win=20; % smooth with a sliding window size
for condi=1:4
    if condi==1; BP_L=BPCondi1L; FB_L=FBCondi1L; BP_R=BPCondi1R;FB_R=FBCondi1R;end
    if condi==2; BP_L=BPCondi2L; FB_L=FBCondi2L; BP_R=BPCondi2R;FB_R=FBCondi2R;end
    if condi==3; BP_L=BPCondi3L; FB_L=FBCondi3L; BP_R=BPCondi3R;FB_R=FBCondi3R;end
    if condi==4; BP_L=BPCondi4L; FB_L=FBCondi4L; BP_R=BPCondi4R;FB_R=FBCondi4R;end 
    % compute
    subplot(4,4,condi); %  autocorr(BP_L)
        autocorr(smoothing(Calinterval(BP_L'),win),430-win); xlabel('lags');ylabel('autocorr'); 
        title({['Condi ' conditionNames{condi}], 'autocorr BP-L'},'Color',condicolors(condi,:));
    subplot(4,4,4+condi); %  xcorr(BP_L vs FB_L)
        time_series1=smoothing(Calinterval(BP_L'),win);% plot(time_series1);
        time_series2=smoothing(Calinterval(FB_L'),win);% plot(time_series2);
        TimeLength=min([length(time_series1) length(time_series2)]);
        time_series1=time_series1-mean(time_series1);% plot(time_series1);
        time_series2=time_series2-mean(time_series2);% plot(time_series2);
        [r,lags]=xcorr(time_series1(1:TimeLength), time_series2(1:TimeLength), 430-win,'normalized');
        plot(lags,r);xlabel('lags');ylabel('Xcorr');xline(0,'--','Color',deepyellow);
        title('xcorr BP-L & FB-L','Color',condicolors(condi,:));
        ylim([-0.5 1]);
    subplot(4,4,4*2+condi); %   xcorr(BP_R vs FB_R) 
        time_series1=smoothing(Calinterval(BP_R'),win);% plot(time_series1);
        time_series2=smoothing(Calinterval(FB_R'),win);% plot(time_series2);
        TimeLength=min([length(time_series1) length(time_series2)]);
        time_series1=time_series1-mean(time_series1);% plot(time_series1);
        time_series2=time_series2-mean(time_series2);% plot(time_series2);
        [r,lags]=xcorr(time_series1(1:TimeLength), time_series2(1:TimeLength), 430-win,'normalized');
        plot(lags,r);xlabel('lags');ylabel('Xcorr');xline(0,'--','Color',deepyellow);
        title('xcorr BP-R & FB-R','Color',condicolors(condi,:));
        ylim([-0.5 1]);
    subplot(4,4,4*3+condi); %  autocorr(BP_R)
        autocorr(smoothing(Calinterval(BP_R'),win),430-win); xlabel('lags');ylabel('autocorr'); 
        title('autocorr BP-R','Color',condicolors(condi,:));
end
suptitle(['BP interval corr smooth win ' num2str(win) ' subject ' num2str(seed)]);
        
figureName=['BPcorr_smooth--' num2str(seed)];
    % save the figure
    saveas(gcf,figureName,'jpg');


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

    %% PLOT 5:  BP interval spectra DFA (long-range / global statistics; strong anticipation)  
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


%% PLOT 6:  BP spectra-DFA smoothed (long-range / global statistics; strong anticipation) on interval
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

%% PLOT 8:  P2 Error smooth (short & long -range statistics)
% Error + autocorr (Error) + Pspectra(Error)+ DFA(Error)
% uncoupled condition (bidiretional error) -subplots(4,10,1-2)
% unidrectional conditions -subplots(4,10,3-4)
% bidirectional condition (bidiretional error) -subplots(4,10,5-10)
addpath D:\360MoveData\Users\alienware\Documents\GitHub\Motor_cordination\1_over_f\data_analysis\2P_testing
addpath D:\360MoveData\Users\alienware\Documents\GitHub\Motor_cordination\1_over_f\data_analysis\DFA

figure('units','normalized','outerposition',[0 0 1 1]);
win=5; % smooth the error with sliding window size
for condi=1:8
%     if condi==1; BP_L=BPCondi1L; BP_R=BPCondi1R; time_series1=find(BP_L==1); time_series2=find(BP_R==1); subcondi=1; end % R's Error base on L (uncouple)
%     if condi==2; BP_L=BPCondi1L; BP_R=BPCondi1R; time_series1=find(BP_R==1); time_series2=find(BP_L==1); subcondi=1; end % L's Error base on R (uncouple)
%     if condi==3; BP_L=BPCondi2L; BP_R=BPCondi2R; time_series1=find(BP_L==1); time_series2=find(BP_R==1); subcondi=2; end % R's Error base on L (unidirectional)
%     if condi==4; BP_L=BPCondi2L; BP_R=BPCondi2R; time_series1=find(BP_R==1); time_series2=find(BP_L==1); subcondi=2; end % L's Error base on R (unidirectional)
%     if condi==5; BP_L=BPCondi3L; BP_R=BPCondi3R; time_series1=find(BP_L==1); time_series2=find(BP_R==1); subcondi=3; end % R's Error base on L (unidirectional)
%     if condi==6; BP_L=BPCondi3L; BP_R=BPCondi3R; time_series1=find(BP_R==1); time_series2=find(BP_L==1); subcondi=3; end % L's Error base on R (unidirectional)
%     if condi==7; BP_L=BPCondi4L; BP_R=BPCondi4R; time_series1=find(BP_L==1); time_series2=find(BP_R==1); subcondi=4; end % R's Error base on L (2Hz bidirection)
%     if condi==8; BP_L=BPCondi4L; BP_R=BPCondi4R; time_series1=find(BP_R==1); time_series2=find(BP_L==1); subcondi=4; end % L's Error base on R (2Hz bidirection)
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

    subplot(4,8,condi); % Error
        plot(ErrorSmooth/2000,'r.');
        xlabel('taps');ylabel('timing error (s)');
        title({['Condi ' conditionNames{subcondi}],'error'},'Color',condicolors(subcondi,:));
        ylim([-0.5 0.6]);
        
    subplot(4,8,8+condi);  % autocorr (Error) 
        autocorr(ErrorSmooth/2000,400); xlabel('lags');ylabel('autocorr'); 
        title('error','Color',condicolors(subcondi,:));
        ylim([-0.2 0.8]);
        
    subplot(4,8,8*2+condi);  % Pspectra(Error)
        y=ErrorSmooth/2000;
        if (condi == 1) | (condi == 2) |(condi == 3) | (condi == 4) | (condi == 5) | (condi == 6) | (condi == 7) | (condi == 8);
            Fs=1.3;% assuming a sampling frequency of 2 Hz
        end
        [freqs,fcoef] = oneoverf(y,Fs);
        xlabel('Log(f)');ylabel('Log(power)');
        title('Spectrum Error','Color',condicolors(subcondi,:));
        ylim([-10 0]);
        
    subplot(4,8,8*3+condi);  %  DFA(Error) 
        [D,Alpha1]=DFA_main(ErrorSmooth/2000);
        title('DFA Error', 'Color',condicolors(subcondi,:));
        ylim([-2 0.5]);

end
suptitle(['error-smooth ' 'win ' num2str(win) ' subject ' num2str(seed)]);

figureName=['Error-smooth -- ' num2str(seed)];
    % save the figure
    saveas(gcf,figureName,'jpg');
    saveas(gcf,figureName,'fig');
    close all;

    
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
    
%%









%%








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

maxfreq=25;sr=2000;
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


%%





















%% PLOT 16:  Xcorr EEG vs BP/FB(boolean) - single chan (short-range): (equal length) (low correlation)
% [ xcorr(EMG_L vs BP_L) + xcorr(EMG_L vs FB_L) + xcorr(EMG_R vs BP_R) + xcorr(EMG_R vs FB_R) ]
% 6 conditions -subplots(4,6,condi) 
cd /ssd/zhibin/1overf/20220518_2P/Segmented_data/XcorrEEG_BP_FB

for chan=1:32
    figure('units','normalized','outerposition',[0 0 1 1]);
    for condi=1:4
        if condi==1; EEG_L=EEGCondi1L; BP_L=BPCondi1L; FB_L=FBCondi1L; EEG_R=EEGCondi1R; BP_R=BPCondi1R; FB_R=FBCondi1R; end
        if condi==2; EEG_L=EEGCondi2L; BP_L=BPCondi2L; FB_L=FBCondi2L; EEG_R=EEGCondi2R; BP_R=BPCondi2R; FB_R=FBCondi2R; end
        if condi==3; EEG_L=EEGCondi3L; BP_L=BPCondi3L; FB_L=FBCondi3L; EEG_R=EEGCondi3R; BP_R=BPCondi3R; FB_R=FBCondi3R; end
        if condi==4; EEG_L=EEGCondi4L; BP_L=BPCondi4L; FB_L=FBCondi4L; EEG_R=EEGCondi4R; BP_R=BPCondi4R; FB_R=FBCondi4R; end
        % compute
        subplot(4,4,condi); %  xcorr(EEG_L vs BP_L)
            [r,lags]=xcorr(EEG_L(:,chan)', BP_L', 1000,'normalized');
            plot(lags./2,r);xlabel('time [ms]');ylabel('Xcorr');xline(0,'Color',deepyellow);
            title({['Condi ' conditionNames{condi}],'EEG-L & BP-L'},'Color',condicolors(condi,:));
            % xlim([-2e3 2e3]);ylim([-0.005 0.1]);
        subplot(4,4,4+condi); %  xcorr(EEG_L vs FB_L)
            [r,lags]=xcorr(EEG_L(:,chan)', FB_L', 1000,'normalized');
            plot(lags./2,r);xlabel('time [ms]');ylabel('Xcorr');xline(0,'Color',deepyellow);
            title('EEG-L & FB-L', 'Color',condicolors(condi,:));
            % xlim([-2e3 2e3]);ylim([-0.005 0.1]);
        subplot(4,4,4*2+condi); %   xcorr(EEG_R vs BP_R)
            [r,lags]=xcorr(EEG_R(:,chan)', BP_R', 1000,'normalized');
            plot(lags./2,r);xlabel('time [ms]');ylabel('Xcorr');xline(0,'Color',deepyellow);
            title('EEG-R & BP-R','Color',condicolors(condi,:));
            % xlim([-2e3 2e3]);ylim([-0.005 0.1]);
        subplot(4,4,4*3+condi); %  xcorr(EEG_R vs FB_R)
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


%%
















%% Photocell
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

