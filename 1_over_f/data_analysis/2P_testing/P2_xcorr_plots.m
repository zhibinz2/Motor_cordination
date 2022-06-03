%% load data
clear;close all;
cd /ssd/zhibin/1overf/20220515_2P
cd /ssd/zhibin/1overf/20220517_2P
cd /ssd/zhibin/1overf/20220518_2P
cd G:\1overf\20220515_2P\Segmented_data\

%% color scheme
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
% test color
showcolor=deepyellow;
imagesc(cat(3,showcolor(1),showcolor(2),showcolor(3)));

%% P2 behaviral (short-range / local statistics; weak anticipation) on interval
% [ autocorr(BP_L) + xcorr(BP_L vs FB_L) + xcorr(BP_R vs FB_R) + autocorr(BP_R) ]
% 6 conditions -subplots(4,6,condi) 

figure('units','normalized','outerposition',[0 0 1 1]);
for condi=1:6
    if condi==1; BP_L=BPCondi1L; FB_L=FBCondi1L; BP_R=BPCondi1R;FB_R=FBCondi1R;end
    if condi==2; BP_L=BPCondi2L; FB_L=FBCondi2L; BP_R=BPCondi2R;FB_R=FBCondi2R;end
    if condi==3; BP_L=BPCondi3L; FB_L=FBCondi3L; BP_R=BPCondi3R;FB_R=FBCondi3R;end
    if condi==4; BP_L=BPCondi4L; FB_L=FBCondi4L; BP_R=BPCondi4R;FB_R=FBCondi4R;end
    if condi==5; BP_L=BPCondi5L; FB_L=FBCondi5L; BP_R=BPCondi5R;FB_R=FBCondi5R;end
    if condi==6; BP_L=BPCondi6L; FB_L=FBCondi6L; BP_R=BPCondi6R;FB_R=FBCondi6R;end
    % compute
    subplot(4,6,condi); %  autocorr(BP_L)
        autocorr(Calinterval(BP_L'),430); xlabel('lags');ylabel('autocorr'); 
        title({['Condi ' num2str(condi)],'autocorr BP-L'},'Color',condicolors(condi,:));
    subplot(4,6,6+condi); %  xcorr(BP_L vs FB_L)
        time_series1=Calinterval(BP_L');% plot(time_series1);
        time_series2=Calinterval(FB_L');% plot(time_series2);
        TimeLength=min([length(time_series1) length(time_series2)]);
        [r,lags]=xcorr(time_series1(1:TimeLength), time_series2(1:TimeLength), 1000,'normalized');
        plot(lags./2,r);xlabel('time [ms]');ylabel('Xcorr');xline(0,'--','Color',deepyellow);
        title('xcorr BP-L & FB-L','Color',condicolors(condi,:));
    subplot(4,6,6*2+condi); %   xcorr(BP_R vs FB_R) 
        time_series1=Calinterval(BP_R');% plot(time_series1);
        time_series2=Calinterval(FB_R');% plot(time_series2);
        TimeLength=min([length(time_series1) length(time_series2)]);
        [r,lags]=xcorr(time_series1(1:TimeLength), time_series2(1:TimeLength), 1000,'normalized');
        plot(lags./2,r);xlabel('time [ms]');ylabel('Xcorr');xline(0,'--','Color',deepyellow);
        title('xcorr BP-R & FB-R','Color',condicolors(condi,:));
    subplot(4,6,6*3+condi); %  autocorr(BP_R)
        autocorr(Calinterval(BP_R'),430); xlabel('lags');ylabel('autocorr'); 
        title('autocorr BP-R','Color',condicolors(condi,:));
end
suptitle(['short-range statistics: tapping intervals' ' subject ' num2str(seed)]);
        
figureName=['ShortBP' num2str(seed)];
    % save the figure
    saveas(gcf,figureName,'jpg');
    
%% P2 behaviral smooth (better!) (short-range / local statistics; weak anticipation) on interval
% [ autocorr(BP_L) + xcorr(BP_L vs FB_L) + xcorr(BP_R vs FB_R) + autocorr(BP_R) ]
% 6 conditions -subplots(4,6,condi) 

win=20; % smooth with a sliding window size

figure('units','normalized','outerposition',[0 0 1 1]);
for condi=1:6
    if condi==1; BP_L=BPCondi1L; FB_L=FBCondi1L; BP_R=BPCondi1R;FB_R=FBCondi1R;end
    if condi==2; BP_L=BPCondi2L; FB_L=FBCondi2L; BP_R=BPCondi2R;FB_R=FBCondi2R;end
    if condi==3; BP_L=BPCondi3L; FB_L=FBCondi3L; BP_R=BPCondi3R;FB_R=FBCondi3R;end
    if condi==4; BP_L=BPCondi4L; FB_L=FBCondi4L; BP_R=BPCondi4R;FB_R=FBCondi4R;end
    if condi==5; BP_L=BPCondi5L; FB_L=FBCondi5L; BP_R=BPCondi5R;FB_R=FBCondi5R;end
    if condi==6; BP_L=BPCondi6L; FB_L=FBCondi6L; BP_R=BPCondi6R;FB_R=FBCondi6R;end
    % compute
    subplot(4,6,condi); %  autocorr(BP_L)
        autocorr(smoothError(Calinterval(BP_L'),win),430-win); xlabel('lags');ylabel('autocorr'); 
        title({['Condi ' num2str(condi)], 'autocorr BP-L'},'Color',condicolors(condi,:));
    subplot(4,6,6+condi); %  xcorr(BP_L vs FB_L)
        time_series1=smoothError(Calinterval(BP_L'),win);% plot(time_series1);
        time_series2=smoothError(Calinterval(FB_L'),win);% plot(time_series2);
        TimeLength=min([length(time_series1) length(time_series2)]);
        [r,lags]=xcorr(time_series1(1:TimeLength), time_series2(1:TimeLength), 1000,'normalized');
        plot(lags./2,r);xlabel('time [ms]');ylabel('Xcorr');xline(0,'--','Color',deepyellow);
        title('xcorr BP-L & FB-L','Color',condicolors(condi,:));
    subplot(4,6,6*2+condi); %   xcorr(BP_R vs FB_R) 
        time_series1=smoothError(Calinterval(BP_R'),win);% plot(time_series1);
        time_series2=smoothError(Calinterval(FB_R'),win);% plot(time_series2);
        TimeLength=min([length(time_series1) length(time_series2)]);
        [r,lags]=xcorr(time_series1(1:TimeLength), time_series2(1:TimeLength), 1000,'normalized');
        plot(lags./2,r);xlabel('time [ms]');ylabel('Xcorr');xline(0,'--','Color',deepyellow);
        title('xcorr BP-R & FB-R','Color',condicolors(condi,:));
    subplot(4,6,6*3+condi); %  autocorr(BP_R)
        autocorr(smoothError(Calinterval(BP_R'),win),430-win); xlabel('lags');ylabel('autocorr'); 
        title('autocorr BP-R','Color',condicolors(condi,:));
end
suptitle(['short-range statistics: tapping intervals smoothed ' 'subject ' num2str(seed)]);
        
figureName=['ShortBP_smooth' num2str(seed)];
    % save the figure
    saveas(gcf,figureName,'jpg');
    
%% P2 behaviral boolean (short-range / local statistics; weak anticipation) on interval
% [ autocorr(BP_L) + xcorr(BP_L vs FB_L) + xcorr(BP_R vs FB_R) + autocorr(BP_R) ]
% 6 conditions -subplots(4,6,condi) 

figure('units','normalized','outerposition',[0 0 1 1]);
for condi=1:6
    if condi==1; BP_L=BPCondi1L; FB_L=FBCondi1L; BP_R=BPCondi1R;FB_R=FBCondi1R;end
    if condi==2; BP_L=BPCondi2L; FB_L=FBCondi2L; BP_R=BPCondi2R;FB_R=FBCondi2R;end
    if condi==3; BP_L=BPCondi3L; FB_L=FBCondi3L; BP_R=BPCondi3R;FB_R=FBCondi3R;end
    if condi==4; BP_L=BPCondi4L; FB_L=FBCondi4L; BP_R=BPCondi4R;FB_R=FBCondi4R;end
    if condi==5; BP_L=BPCondi5L; FB_L=FBCondi5L; BP_R=BPCondi5R;FB_R=FBCondi5R;end
    if condi==6; BP_L=BPCondi6L; FB_L=FBCondi6L; BP_R=BPCondi6R;FB_R=FBCondi6R;end
    % compute
    subplot(4,6,condi); %  autocorr(BP_L)
        autocorr(BP_L',430); xlabel('lags');ylabel('autocorr'); 
        title({['Condi ' num2str(condi)], 'autocorr BP-L'},'Color',condicolors(condi,:));
    subplot(4,6,6+condi); %  xcorr(BP_L vs FB_L)
        time_series1=BP_L';% plot(time_series1);
        time_series2=FB_L';% plot(time_series2);
        TimeLength=min([length(time_series1) length(time_series2)]);
        [r,lags]=xcorr(time_series1(1:TimeLength), time_series2(1:TimeLength), 1000,'normalized');
        plot(lags./2,r);xlabel('time [ms]');ylabel('Xcorr');xline(0,'--','Color',deepyellow);
        title('xcorr BP-L & FB-L','Color',condicolors(condi,:));
    subplot(4,6,6*2+condi); %   xcorr(BP_R vs FB_R) 
        time_series1=BP_R';% plot(time_series1);
        time_series2=FB_R';% plot(time_series2);
        TimeLength=min([length(time_series1) length(time_series2)]);
        [r,lags]=xcorr(time_series1(1:TimeLength), time_series2(1:TimeLength), 1000,'normalized');
        plot(lags./2,r);xlabel('time [ms]');ylabel('Xcorr');xline(0,'--','Color',deepyellow);
        title('xcorr BP-R & FB-R','Color',condicolors(condi,:));
    subplot(4,6,6*3+condi); %  autocorr(BP_R)
        autocorr(BP_R',430); xlabel('lags');ylabel('autocorr'); 
        title('autocorr BP-R','Color',condicolors(condi,:));
end
suptitle(['short-range statistics: tapping boolean ' 'subject ' num2str(seed)]);
        
figureName=['ShortBP_boolean' num2str(seed)];
    % save the figure
    saveas(gcf,figureName,'jpg');
    
%% P2 behaviral (long-range / global statistics; strong anticipation) on interval
% Pspectra(BP_L) + Pspectra(BP_R) + DFA(BP_L) + DFA(BP_R)
% 6 conditions -subplots(4,6,i)
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing

figure('units','normalized','outerposition',[0 0 1 1]);
for condi=1:6
    if condi==1; BP_L=BPCondi1L; BP_R=BPCondi1R; end
    if condi==2; BP_L=BPCondi2L; BP_R=BPCondi2R; end
    if condi==3; BP_L=BPCondi3L; BP_R=BPCondi3R; end
    if condi==4; BP_L=BPCondi4L; BP_R=BPCondi4R; end
    if condi==5; BP_L=BPCondi5L; BP_R=BPCondi5R; end
    if condi==6; BP_L=BPCondi6L; BP_R=BPCondi6R; end
    subplot(4,6,condi); % Pspectra(BP_L)
        y=Calinterval(BP_L');
        if (condi == 1) | (condi == 2) |(condi == 3) | (condi == 4);
            Fs=2;% assuming a sampling frequency of 2 Hz
        else (condi == 5) | (condi == 6) ;
            Fs=3;% assuming a sampling frequency of 3 Hz
        end
        [freqs,fcoef] = oneoverf(y,Fs);
        xlabel('Log(f)');ylabel('Log(power)');
        title({['Condi ' num2str(condi)], 'Spectrum BP_L'},'Color',condicolors(condi,:));
        ylim([-2 8]);xlim([-3 0]);
        
    subplot(4,6,6+condi);  % Pspectra(BP_R) 
        y=Calinterval(BP_R');
        if (condi == 1) | (condi == 2) |(condi == 3) | (condi == 4);
            Fs=2;% assuming a sampling frequency of 2 Hz
        else (condi == 5) | (condi == 6) ;
            Fs=3;% assuming a sampling frequency of 3 Hz
        end
        [freqs,fcoef] = oneoverf(y,Fs);
        xlabel('Log(f)');ylabel('Log(power)');
        title('Spectrum BP_R','Color',condicolors(condi,:));  
        ylim([-2 8]);xlim([-3 0]);
        
    subplot(4,6,6*2+condi);  % DFA(BP_L)
        [D,Alpha1]=DFA_main(Calinterval(BP_L'));
        title('DFA BP-L','Color',condicolors(condi,:));
        ylim([0.5 4]);xlim([0.5 3]);
        
    subplot(4,6,6*3+condi);  % DFA(BP_R)
        [D,Alpha1]=DFA_main(Calinterval(BP_R'));
        title('DFA BP-R','Color',condicolors(condi,:));
        ylim([0.5 4]);xlim([0.5 3]);

end
suptitle(['long-range statistics: tapping intervals ' 'subject ' num2str(seed)]);

figureName=['LongBP' num2str(seed)];
    % save the figure
    saveas(gcf,figureName,'jpg');

    
%% P2 behaviral smoothed (long-range / global statistics; strong anticipation) on interval
% Pspectra(BP_L) + Pspectra(BP_R) + DFA(BP_L) + DFA(BP_R)
% 6 conditions -subplots(4,6,i)
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing

win=20; % smooth with a sliding window size

figure('units','normalized','outerposition',[0 0 1 1]);
for condi=1:6
    if condi==1; BP_L=BPCondi1L; BP_R=BPCondi1R; end
    if condi==2; BP_L=BPCondi2L; BP_R=BPCondi2R; end
    if condi==3; BP_L=BPCondi3L; BP_R=BPCondi3R; end
    if condi==4; BP_L=BPCondi4L; BP_R=BPCondi4R; end
    if condi==5; BP_L=BPCondi5L; BP_R=BPCondi5R; end
    if condi==6; BP_L=BPCondi6L; BP_R=BPCondi6R; end
    subplot(4,6,condi); % Pspectra(BP_L)
        y=smoothError(Calinterval(BP_L'),win);
        if (condi == 1) | (condi == 2) |(condi == 3) | (condi == 4);
            Fs=2;% assuming a sampling frequency of 2 Hz
        else (condi == 5) | (condi == 6) ;
            Fs=3;% assuming a sampling frequency of 3 Hz
        end
        [freqs,fcoef] = oneoverf(y,Fs);
        xlabel('Log(f)');ylabel('Log(power)');
        title({['Condi ' num2str(condi)], 'Spectrum BP_L'},'Color',condicolors(condi,:));
        ylim([-4 8]);xlim([-3 0]);
        
    subplot(4,6,6+condi);  % Pspectra(BP_R) 
        y=smoothError(Calinterval(BP_R'),win);
        if (condi == 1) | (condi == 2) |(condi == 3) | (condi == 4);
            Fs=2;% assuming a sampling frequency of 2 Hz
        else (condi == 5) | (condi == 6) ;
            Fs=3;% assuming a sampling frequency of 3 Hz
        end
        [freqs,fcoef] = oneoverf(y,Fs);
        xlabel('Log(f)');ylabel('Log(power)');
        title('Spectrum BP_R', 'Color',condicolors(condi,:));    
        
    subplot(4,6,6*2+condi);  % DFA(BP_L)
        [D,Alpha1]=DFA_main(smoothError(Calinterval(BP_L'),win));
        title('DFA BP_L','Color',condicolors(condi,:));
        ylim([0 4]);xlim([0.5 3]);
        
    subplot(4,6,6*3+condi);  % DFA(BP_R)
        [D,Alpha1]=DFA_main(smoothError(Calinterval(BP_R'),win));
        title('DFA BP_R','Color',condicolors(condi,:));
        ylim([0 4]);xlim([0.5 3]);

end
suptitle(['long-range statistics: tapping intervals smoothed ' 'subject ' num2str(seed)]);

figureName=['LongBP_smooth' num2str(seed)];
    % save the figure
    saveas(gcf,figureName,'jpg');

%% P2 behaviral boolean (long-range / global statistics; strong anticipation) on interval (skip, take time too long to finish)
% Pspectra(BP_L) + Pspectra(BP_R) + DFA(BP_L) + DFA(BP_R)
% 6 conditions -subplots(4,6,i)
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing
tic
figure('units','normalized','outerposition',[0 0 1 1]);
for condi=1:6
    if condi==1; BP_L=BPCondi1L; BP_R=BPCondi1R; end
    if condi==2; BP_L=BPCondi2L; BP_R=BPCondi2R; end
    if condi==3; BP_L=BPCondi3L; BP_R=BPCondi3R; end
    if condi==4; BP_L=BPCondi4L; BP_R=BPCondi4R; end
    if condi==5; BP_L=BPCondi5L; BP_R=BPCondi5R; end
    if condi==6; BP_L=BPCondi6L; BP_R=BPCondi6R; end
    subplot(4,6,condi); % Pspectra(BP_L)
        y=BP_L';
        if (condi == 1) | (condi == 2) |(condi == 3) | (condi == 4);
            Fs=2;% assuming a sampling frequency of 2 Hz
        else (condi == 5) | (condi == 6) ;
            Fs=3;% assuming a sampling frequency of 3 Hz
        end
        [freqs,fcoef] = oneoverf(y,Fs);
        xlabel('Log(f)');ylabel('Log(power)');
        title({'Spectrum BP_L', ['Condi ' num2str(condi)]});
        
    subplot(4,6,6+condi);  % Pspectra(BP_R) 
        y=BP_R';
        if (condi == 1) | (condi == 2) |(condi == 3) | (condi == 4);
            Fs=2;% assuming a sampling frequency of 2 Hz
        else (condi == 5) | (condi == 6) ;
            Fs=3;% assuming a sampling frequency of 3 Hz
        end
        [freqs,fcoef] = oneoverf(y,Fs);
        xlabel('Log(f)');ylabel('Log(power)');
        title({'Spectrum BP_R', ['Condi ' num2str(condi)]});    
        
    subplot(4,6,6*2+condi);  % DFA(BP_L)
        [D,Alpha1]=DFA_main(BP_L');
        title({'DFA BP_L', ['Condi ' num2str(condi)]});
        
    subplot(4,6,6*3+condi);  % DFA(BP_R)
        [D,Alpha1]=DFA_main(BP_R');
        title({'DFA BP_R', ['Condi ' num2str(condi)]});

end
suptitle({'long-range statistics: tapping intervals',['subject ' num2str(seed)]});

figureName=['LongBP_boolean' num2str(seed)];
    % save the figure
    saveas(gcf,figureName,'jpg');
toc
%% P2 Error (short & long -range statistics)
% Error + autocorr (Error) + Pspectra(Error)+ DFA(Error)
% uncoupled condition (bidiretional error) -subplots(4,10,1-2)
% unidrectional conditions -subplots(4,10,3-4)
% bidirectional condition (bidiretional error) -subplots(4,10,5-10)
addpath D:\360MoveData\Users\alienware\Documents\GitHub\Motor_cordination\1_over_f\data_analysis\2P_testing
addpath D:\360MoveData\Users\alienware\Documents\GitHub\Motor_cordination\1_over_f\data_analysis\DFA

win=5; % smooth the error with sliding window size

figure('units','normalized','outerposition',[0 0 1 1]);
for condi=1:12
%     if condi==1; BP_L=BPCondi1L; BP_R=BPCondi1R; time_series1=find(BP_L==1); time_series2=find(BP_R==1); subcondi=1; end % R's Error base on L (uncouple)
%     if condi==2; BP_L=BPCondi1L; BP_R=BPCondi1R; time_series1=find(BP_R==1); time_series2=find(BP_L==1); subcondi=1; end % L's Error base on R (uncouple)
%     if condi==3; BP_L=BPCondi2L; BP_R=BPCondi2R; time_series1=find(BP_L==1); time_series2=find(BP_R==1); subcondi=2; end % R's Error base on L (unidirectional)
%     if condi==4; BP_L=BPCondi2L; BP_R=BPCondi2R; time_series1=find(BP_R==1); time_series2=find(BP_L==1); subcondi=2; end % L's Error base on R (unidirectional)
%     if condi==5; BP_L=BPCondi3L; BP_R=BPCondi3R; time_series1=find(BP_L==1); time_series2=find(BP_R==1); subcondi=3; end % R's Error base on L (unidirectional)
%     if condi==6; BP_L=BPCondi3L; BP_R=BPCondi3R; time_series1=find(BP_R==1); time_series2=find(BP_L==1); subcondi=3; end % L's Error base on R (unidirectional)
%     if condi==7; BP_L=BPCondi4L; BP_R=BPCondi4R; time_series1=find(BP_L==1); time_series2=find(BP_R==1); subcondi=4; end % R's Error base on L (2Hz bidirection)
%     if condi==8; BP_L=BPCondi4L; BP_R=BPCondi4R; time_series1=find(BP_R==1); time_series2=find(BP_L==1); subcondi=4; end % L's Error base on R (2Hz bidirection)
%     if condi==9; BP_L=BPCondi5L; BP_R=BPCondi5R; time_series1=find(BP_L==1); time_series2=find(BP_R==1); subcondi=5; end % R's Error base on L (3Hz bidirection)
%     if condi==10; BP_L=BPCondi5L; BP_R=BPCondi5R; time_series1=find(BP_R==1); time_series2=find(BP_L==1); subcondi=5; end % L's Error base on R (3Hz bidirection)
%     if condi==11; BP_L=BPCondi6L; BP_R=BPCondi6R; time_series1=find(BP_L==1); time_series2=find(BP_R==1); subcondi=6; end % R's Error base on L (>3Hz bidirection)
%     if condi==12; BP_L=BPCondi6L; BP_R=BPCondi6R; time_series1=find(BP_R==1); time_series2=find(BP_L==1); subcondi=6; end % L's Error base on R (>3Hz bidirection)
    if condi==1; BP_L=BPCondi1L; FB_L=FBCondi1L; time_series1=find(BP_L==1); time_series2=find(FB_L==1); subcondi=1; end % R's Error base on L (uncouple)
    if condi==2; BP_R=BPCondi1R; FB_R=FBCondi1R; time_series1=find(BP_R==1); time_series2=find(FB_R==1); subcondi=1; end % L's Error base on R (uncouple)
    if condi==3; BP_L=BPCondi2L; FB_L=FBCondi2L; time_series1=find(BP_L==1); time_series2=find(FB_L==1); subcondi=2; end % R's Error base on L (unidirectional)
    if condi==4; BP_R=BPCondi2R; FB_R=FBCondi2R; time_series1=find(BP_R==1); time_series2=find(FB_R==1); subcondi=2; end % L's Error base on R (unidirectional)
    if condi==5; BP_L=BPCondi3L; FB_L=FBCondi3L; time_series1=find(BP_L==1); time_series2=find(FB_L==1); subcondi=3; end % R's Error base on L (unidirectional)
    if condi==6; BP_R=BPCondi3R; FB_R=FBCondi3R; time_series1=find(BP_R==1); time_series2=find(FB_R==1); subcondi=3; end % L's Error base on R (unidirectional)
    if condi==7; BP_L=BPCondi4L; FB_L=FBCondi4L; time_series1=find(BP_L==1); time_series2=find(FB_L==1); subcondi=4; end % R's Error base on L (2Hz bidirection)
    if condi==8; BP_R=BPCondi4R; FB_R=FBCondi4R; time_series1=find(BP_R==1); time_series2=find(FB_R==1); subcondi=4; end % L's Error base on R (2Hz bidirection)
    if condi==9; BP_L=BPCondi5L; FB_L=FBCondi5L; time_series1=find(BP_L==1); time_series2=find(FB_L==1); subcondi=5; end % R's Error base on L (3Hz bidirection)
    if condi==10; BP_R=BPCondi5R; FB_R=FBCondi5R; time_series1=find(BP_R==1); time_series2=find(FB_R==1); subcondi=5; end % L's Error base on R (3Hz bidirection)
    if condi==11; BP_L=BPCondi6L; FB_L=FBCondi6L; time_series1=find(BP_L==1); time_series2=find(FB_L==1); subcondi=6; end % R's Error base on L (>3Hz bidirection)
    if condi==12; BP_R=BPCondi6R; FB_R=FBCondi6R; time_series1=find(BP_R==1); time_series2=find(FB_R==1); subcondi=6; end % L's Error base on R (>3Hz bidirection)
    % Compute error
    % plot(time_series1); % plot(time_series2);
    % TimeLength=min([length(time_series1) length(time_series2)]);
    % TimeLength=max([length(time_series1) length(time_series2)]);
    TimeLength=length(time_series1); % measure each error in time_series1
    Error=[];
    % find the time difference with the closest botton press 
    for i=1:TimeLength % i=StartStim2 % syncopation has one press less than condition 1 and 3
        [minValue,closetIndex]=min(abs(time_series2-time_series1(i))); % closetIndex in time_series2
        Error(i)=time_series2(closetIndex)-time_series1(i);
    end
    % smooth the error with sliding window
    [ErrorSmooth] = smoothError(Error,win);
%     ErrorSmooth=Error;

    subplot(4,12,condi); % Error
        plot(ErrorSmooth/2000,'r.');
        xlabel('taps');ylabel('timing error (s)');
        title({['Condi ' num2str(subcondi)],'error'},'Color',condicolors(subcondi,:));
        ylim([-0.5 0.6]);
        
    subplot(4,12,12+condi);  % autocorr (Error) 
        autocorr(ErrorSmooth/2000,400); xlabel('lags');ylabel('autocorr'); 
        title('error','Color',condicolors(subcondi,:));
        ylim([-0.2 0.8]);
        
    subplot(4,12,12*2+condi);  % Pspectra(Error)
        y=ErrorSmooth/2000;
        if (condi == 1) | (condi == 2) |(condi == 3) | (condi == 4) | (condi == 5) | (condi == 6) | (condi == 7) | (condi == 8);
            Fs=2;% assuming a sampling frequency of 2 Hz
        else (condi == 9) | (condi == 10) | (condi == 11) | (condi == 12);
            Fs=3;% assuming a sampling frequency of 3 Hz
        end
        [freqs,fcoef] = oneoverf(y,Fs);
        xlabel('Log(f)');ylabel('Log(power)');
        title('Spectrum Error','Color',condicolors(subcondi,:));
        ylim([-10 0]);
        
    subplot(4,12,12*3+condi);  %  DFA(Error) 
        [D,Alpha1]=DFA_main(ErrorSmooth/2000);
        title('DFA Error', 'Color',condicolors(subcondi,:));
        % ylim([-3 1]);

end
suptitle(['short & long-range statistics: error-smooth ' 'win ' num2str(win) ' subject ' num2str(seed)]);

figureName=['Error-smooth' num2str(seed)];
    % save the figure
    saveas(gcf,figureName,'jpg');
    

%% P2 EMG (short-range / local statistics; weak anticipation)
% [ autocorr(EMG_L) + xcorr (EMG_L vs EMG_R) + autocorr(EMG_R) ]
% 6 conditions - subplots(3,6,i)
figure('units','normalized','outerposition',[0 0 1 1]);
for condi=1:6
    if condi==1; EMG_L=EMGCondi1L; EMG_R=EMGCondi1R; end % R's Error base on L (uncouple)
    if condi==2; EMG_L=EMGCondi2L; EMG_R=EMGCondi2R; end % L's Error base on R (uncouple)
    if condi==3; EMG_L=EMGCondi3L; EMG_R=EMGCondi3R; end % R's Error base on L (unidirectional)
    if condi==4; EMG_L=EMGCondi4L; EMG_R=EMGCondi4R; end % L's Error base on R (unidirectional)
    if condi==5; EMG_L=EMGCondi5L; EMG_R=EMGCondi5R; end % R's Error base on L (2Hz bidirection)
    if condi==6; EMG_L=EMGCondi6L; EMG_R=EMGCondi6R; end % L's Error base on R (2Hz bidirection)
    
    subplot(3,6,condi); % autocorr(EMG_L)
        lags=round(length(EMG_L))-1;
        autocorr(EMG_L',lags); xlabel('lags of samples');ylabel('autocorr'); 
        title({['Condi ' num2str(condi)],'EMG-L'},'Color',condicolors(condi,:));
        ylim([-0.5 1]);
        xlim([0 2000]);
        
    subplot(3,6,6+condi); % xcorr (EMG_L vs EMG_R) 
        time_series1=EMG_L';% plot(time_series1);
        time_series2=EMG_R';% plot(time_series2);
        TimeLength=min([length(time_series1) length(time_series2)]);
        [r,lags]=xcorr(time_series1(1:TimeLength), time_series2(1:TimeLength), TimeLength-1,'normalized');
        plot(lags./2,r);xlabel('lag of time [ms]');ylabel('Xcorr');
        title('EMG-L & EMG-R','Color',condicolors(condi,:));
        xline(0,'--','Color',deepyellow);
        ylim([-0.1 0.6]);% xlim([-5e5 5e5]);
        xlim([-2000 2000]);
        
    subplot(3,6,6*2+condi); %  autocorr(EMG_R) 
        lags=round(length(EMG_R))-1;
        autocorr(EMG_R',lags); xlabel('lags of samples');ylabel('autocorr'); 
        title('EMG-R','Color',condicolors(condi,:));
        ylim([-0.5 1]);
        xlim([0 2000]);
end
suptitle(['short-range statistics: EMG ' 'subject ' num2str(seed)]);

figureName=['ShortEMG' num2str(seed)];
    % save the figure
    saveas(gcf,figureName,'jpg');
    
%% P2 EMG (long-range / global statistics; strong anticipation)
% [ Pspectra(EMG_L) + Xspectra (EMG_L vs EMG_R) + Pspectra(EMG_R) ]
% 6 conditions - subplots(3,6,i)
maxfreq=10;sr=2000;
figure('units','normalized','outerposition',[0 0 1 1]);
for condi=1:6
    if condi==1; EMG_L=EMGCondi1L; EMG_R=EMGCondi1R; end % R's Error base on L (uncouple)
    if condi==2; EMG_L=EMGCondi2L; EMG_R=EMGCondi2R; end % L's Error base on R (uncouple)
    if condi==3; EMG_L=EMGCondi3L; EMG_R=EMGCondi3R; end % R's Error base on L (unidirectional)
    if condi==4; EMG_L=EMGCondi4L; EMG_R=EMGCondi4R; end % L's Error base on R (unidirectional)
    if condi==5; EMG_L=EMGCondi5L; EMG_R=EMGCondi5R; end % R's Error base on L (2Hz bidirection)
    if condi==6; EMG_L=EMGCondi6L; EMG_R=EMGCondi6R; end % L's Error base on R (2Hz bidirection)
    
    subplot(3,6,condi); % Pspectra(EMG_L) 
        [fcoef1,fcoef2,cprod, freqs] = spectra(EMG_L,EMG_L,maxfreq,sr);
        plot(freqs,abs(cprod));xlabel('freqs');ylabel('power');
        title({['Condi ' num2str(condi)], 'EMG-L'},'Color',condicolors(condi,:));
        xticks(linspace(0,10,5));xticklabels({'0','2.5','5','7.5','10'});
%         ylim([0 3e-3]);
        
    subplot(3,6,6+condi); % Xspectra (EMG_L vs EMG_R) 
        [fcoef1,fcoef2,cprod, freqs] = spectra(EMG_L,EMG_R,maxfreq,sr);
        plot(freqs,abs(cprod));xlabel('freqs');ylabel('xspectra');
        title('EMG-L & EMG-R','Color',condicolors(condi,:));
        xticks(linspace(0,10,5));xticklabels({'0','2.5','5','7.5','10'});
%         ylim([0 6]);
        
    subplot(3,6,6*2+condi); %  Pspectra (EMG_R) 
        [fcoef1,fcoef2,cprod, freqs] = spectra(EMG_R,EMG_R,maxfreq,sr);
        plot(freqs,abs(cprod));xlabel('freqs');ylabel('power');
        title('EMG-R','Color',condicolors(condi,:));
        xticks(linspace(0,10,5));xticklabels({'0','2.5','5','7.5','10'});
%         ylim([0 2e-3]);
end
suptitle(['long-range statistics: EMG ' 'subject ' num2str(seed)]);

figureName=['LongEMG' num2str(seed)];
    % save the figure
    saveas(gcf,figureName,'jpg');
    
%% EEG auto xcorr (weak anticipation, short-range / local statistics)
% [ autocorr(EEG_L) + xcorr(EEG_L vs EEG_R) + autocorr(EEG_R) ]
% each of the 6 conditions (subplots(3,6,i)) for each of the 32 channels;
% auto save the 32 figures in a folder named - ShortEEG
cd /ssd/zhibin/1overf/20220515_2P/Segmented_data/ShortEEG
cd /ssd/zhibin/1overf/20220517_2P/Segmented_data/ShortEEG
cd /ssd/zhibin/1overf/20220518_2P/Segmented_data/ShortEEG

maxfreq=25;sr=2000;
for chan = 1:32
    figure('units','normalized','outerposition',[0 0 1 1]);
    for condi=1:6
        if condi==1; EEG_L=EEGCondi1L; EEG_R=EEGCondi1R; end % R's Error base on L (uncouple)
        if condi==2; EEG_L=EEGCondi2L; EEG_R=EEGCondi2R; end % L's Error base on R (uncouple)
        if condi==3; EEG_L=EEGCondi3L; EEG_R=EEGCondi3R; end % R's Error base on L (unidirectional)
        if condi==4; EEG_L=EEGCondi4L; EEG_R=EEGCondi4R; end % L's Error base on R (unidirectional)
        if condi==5; EEG_L=EEGCondi5L; EEG_R=EEGCondi5R; end % R's Error base on L (2Hz bidirection)
        if condi==6; EEG_L=EEGCondi6L; EEG_R=EEGCondi6R; end % L's Error base on R (2Hz bidirection)

        subplot(3,6,condi); % autocorr(EEG_L) 
            lags=round(length(EEG_L(:,chan)))-1;
            autocorr(EEG_L(:,chan)',lags); xlabel('lags of samples');ylabel('autocorr'); 
            title({['Condi ' num2str(condi)],'EEG-L'},'Color',condicolors(condi,:));
            ylim([-0.5 1]);
            xlim([0 2000]);

        subplot(3,6,6+condi); % xcorr(EEG_L vs EEG_R)  
            time_series1=EEG_L(:,chan)';% plot(time_series1);
            time_series2=EEG_R(:,chan)';% plot(time_series2);
            TimeLength=min([length(time_series1) length(time_series2)]);
            [r,lags]=xcorr(time_series1(1:TimeLength), time_series2(1:TimeLength), TimeLength-1,'normalized');
            plot(lags./2,r);xlabel('lag of time [ms]');ylabel('Xcorr');
            title('EEG-L & EEG-R','Color',condicolors(condi,:));
            ylim([-0.3 0.3]);
            xlim([-2000 2000]);

        subplot(3,6,6*2+condi); % autocorr(EEG_R)
            lags=round(length(EEG_R(:,chan)))-1;
            autocorr(EEG_R(:,chan)',lags); xlabel('lags of samples');ylabel('autocorr'); 
            title('EEG-R', 'Color',condicolors(condi,:));
            ylim([-0.5 1]);
            xlim([0 2000]);
    end
    suptitle(['short-range statistics: EEG channel ' num2str(chan) ' subject ' num2str(seed)]);
    figureName=['ShortEEGchan' num2str(chan)];
    % save the figure
    saveas(gcf,figureName,'jpg');
    close all;
end


%% EEG Power spectra and Xspectra (strong anticipation, long-range / global statistics)
% [ Pspectra(EEG_L) + Xspectra(EEG_L vs EEG_R) + Pspectra(EEG_R) ]
% each of the 6 conditions (subplots(3,6,i)) for each of the 32 channels;
% auto save the 32 figures in a folder named - LongEEG
cd /ssd/zhibin/1overf/20220515_2P/Segmented_data/LongEEG
cd /ssd/zhibin/1overf/20220517_2P/Segmented_data/LongEEG
cd /ssd/zhibin/1overf/20220518_2P/Segmented_data/LongEEG

maxfreq=25;sr=2000;
for chan = 1:32
    figure('units','normalized','outerposition',[0 0 1 1]);
    for condi=1:6
        if condi==1; EEG_L=EEGCondi1L; EEG_R=EEGCondi1R; end % R's Error base on L (uncouple)
        if condi==2; EEG_L=EEGCondi2L; EEG_R=EEGCondi2R; end % L's Error base on R (uncouple)
        if condi==3; EEG_L=EEGCondi3L; EEG_R=EEGCondi3R; end % R's Error base on L (unidirectional)
        if condi==4; EEG_L=EEGCondi4L; EEG_R=EEGCondi4R; end % L's Error base on R (unidirectional)
        if condi==5; EEG_L=EEGCondi5L; EEG_R=EEGCondi5R; end % R's Error base on L (2Hz bidirection)
        if condi==6; EEG_L=EEGCondi6L; EEG_R=EEGCondi6R; end % L's Error base on R (2Hz bidirection)

        subplot(3,6,condi); % Pspectra(EEG_L) 
            [fcoef1,fcoef2,cprod, freqs] = spectra(EEG_L(:,chan),EEG_L(:,chan),maxfreq,sr);
            plot(freqs,abs(cprod));xlabel('freqs');ylabel('power');
            title({['Condi ' num2str(condi)],'EEG-L'},'Color',condicolors(condi,:));
    %         ylim([]);

        subplot(3,6,6+condi); % xcorr (EEG_L vs EEG_R) 
            [fcoef1,fcoef2,cprod, freqs] = spectra(EEG_L(:,chan),EEG_R(:,chan),maxfreq,sr);
            plot(freqs,abs(cprod));xlabel('freqs');ylabel('Xspectra');
            title('EEG-L & EEG-R','Color',condicolors(condi,:));
    %         ylim([]);

        subplot(3,6,6*2+condi); %  Pspectra(EEG_R) 
            [fcoef1,fcoef2,cprod, freqs] = spectra(EEG_R(:,chan),EEG_R(:,chan),maxfreq,sr);
            plot(freqs,abs(cprod));xlabel('freqs');ylabel('power');
            title('EEG-R','Color',condicolors(condi,:));
    %         ylim([]);
    end
    suptitle(['long-range statistics: EEG channel ' num2str(chan) ' subject ' num2str(seed)]);
    figureName=['LongEEGchan' num2str(chan)];
    % save the figure
    saveas(gcf,figureName,'jpg');
    close all;
end


%% Xcorr (short-range): EMG vs FB(boolean) (equal length)
% [ xcorr(EMG_L vs BP_L) + xcorr(EMG_L vs FB_L) + xcorr(EMG_R vs BP_R) + xcorr(EMG_R vs FB_R) ]
% 6 conditions -subplots(4,6,condi) 

figure('units','normalized','outerposition',[0 0 1 1]);
for condi=1:6
    if condi==1; EMG_L=EMGCondi1L; BP_L=BPCondi1L; FB_L=FBCondi1L; EMG_R=EMGCondi1R; BP_R=BPCondi1R; FB_R=FBCondi1R; end
    if condi==2; EMG_L=EMGCondi2L; BP_L=BPCondi2L; FB_L=FBCondi2L; EMG_R=EMGCondi2R; BP_R=BPCondi2R; FB_R=FBCondi2R; end
    if condi==3; EMG_L=EMGCondi3L; BP_L=BPCondi3L; FB_L=FBCondi3L; EMG_R=EMGCondi3R; BP_R=BPCondi3R; FB_R=FBCondi3R; end
    if condi==4; EMG_L=EMGCondi4L; BP_L=BPCondi4L; FB_L=FBCondi4L; EMG_R=EMGCondi4R; BP_R=BPCondi4R; FB_R=FBCondi4R; end
    if condi==5; EMG_L=EMGCondi5L; BP_L=BPCondi5L; FB_L=FBCondi5L; EMG_R=EMGCondi5R; BP_R=BPCondi5R; FB_R=FBCondi5R; end
    if condi==6; EMG_L=EMGCondi6L; BP_L=BPCondi6L; FB_L=FBCondi6L; EMG_R=EMGCondi6R; BP_R=BPCondi6R; FB_R=FBCondi6R; end
    % compute
    subplot(4,6,condi); %  xcorr(EMG_L vs BP_L)
        [r,lags]=xcorr(EMG_L', BP_L', length(EMG_L),'normalized');
        plot(lags./2,r);xlabel('time [ms]');ylabel('Xcorr');
        title({['Condi ' num2str(condi)],'EMG-L & BP-L'},'Color',condicolors(condi,:));
        xlim([-2e3 2e3]);ylim([-0.005 0.1]);
    subplot(4,6,6+condi); %  xcorr(EMG_L vs FB_L)
        [r,lags]=xcorr(EMG_L', FB_L', length(EMG_L),'normalized');
        plot(lags./2,r);xlabel('time [ms]');ylabel('Xcorr');
        title('EMG-L & FB-L','Color',condicolors(condi,:));
        xlim([-2e3 2e3]);ylim([-0.005 0.1]);
    subplot(4,6,6*2+condi); %   xcorr(EMG_R vs BP_R)
        [r,lags]=xcorr(EMG_R', BP_R', length(EMG_R),'normalized');
        plot(lags./2,r);xlabel('time [ms]');ylabel('Xcorr');
        title('EMG-R & BP-R','Color',condicolors(condi,:));
        xlim([-2e3 2e3]);ylim([-0.005 0.1]);
    subplot(4,6,6*3+condi); %  xcorr(EMG_R vs FB_R)
        [r,lags]=xcorr(EMG_R', FB_R', length(EMG_R),'normalized');
        plot(lags./2,r);xlabel('time [ms]');ylabel('Xcorr');
        title('EMG-R & FB-R','Color',condicolors(condi,:));
        xlim([-2e3 2e3]);ylim([-0.005 0.1]);
end
suptitle(['xcorr: EMG-BP/FB ' 'subject ' num2str(seed)]);
        
figureName=['xcorrEMG-BPFBboolean' num2str(seed)];
    % save the figure
    saveas(gcf,figureName,'jpg');
    
%% Xcorr (short-range): EEG vs FB(boolean) (equal length) (low correlation)
% [ xcorr(EMG_L vs BP_L) + xcorr(EMG_L vs FB_L) + xcorr(EMG_R vs BP_R) + xcorr(EMG_R vs FB_R) ]
% 6 conditions -subplots(4,6,condi) 
cd /ssd/zhibin/1overf/20220518_2P/Segmented_data/XcorrEEG_BP_FB

for chan=1:32
    figure('units','normalized','outerposition',[0 0 1 1]);
    for condi=1:6
        if condi==1; EEG_L=EEGCondi1L; BP_L=BPCondi1L; FB_L=FBCondi1L; EEG_R=EEGCondi1R; BP_R=BPCondi1R; FB_R=FBCondi1R; end
        if condi==2; EEG_L=EEGCondi2L; BP_L=BPCondi2L; FB_L=FBCondi2L; EEG_R=EEGCondi2R; BP_R=BPCondi2R; FB_R=FBCondi2R; end
        if condi==3; EEG_L=EEGCondi3L; BP_L=BPCondi3L; FB_L=FBCondi3L; EEG_R=EEGCondi3R; BP_R=BPCondi3R; FB_R=FBCondi3R; end
        if condi==4; EEG_L=EEGCondi4L; BP_L=BPCondi4L; FB_L=FBCondi4L; EEG_R=EEGCondi4R; BP_R=BPCondi4R; FB_R=FBCondi4R; end
        if condi==5; EEG_L=EEGCondi5L; BP_L=BPCondi5L; FB_L=FBCondi5L; EEG_R=EEGCondi5R; BP_R=BPCondi5R; FB_R=FBCondi5R; end
        if condi==6; EEG_L=EEGCondi6L; BP_L=BPCondi6L; FB_L=FBCondi6L; EEG_R=EEGCondi6R; BP_R=BPCondi6R; FB_R=FBCondi6R; end
        % compute
        subplot(4,6,condi); %  xcorr(EEG_L vs BP_L)
            [r,lags]=xcorr(EEG_L(:,chan)', BP_L', length(EEG_L),'normalized');
            plot(lags./2,r);xlabel('time [ms]');ylabel('Xcorr');
            title({['Condi ' num2str(condi)],'EEG-L & BP-L'},'Color',condicolors(condi,:));
            xlim([-2e3 2e3]);ylim([-0.005 0.1]);
        subplot(4,6,6+condi); %  xcorr(EEG_L vs FB_L)
            [r,lags]=xcorr(EEG_L(:,chan)', FB_L', length(EEG_L),'normalized');
            plot(lags./2,r);xlabel('time [ms]');ylabel('Xcorr');
            title('EEG-L & FB-L', 'Color',condicolors(condi,:));
            xlim([-2e3 2e3]);ylim([-0.005 0.1]);
        subplot(4,6,6*2+condi); %   xcorr(EEG_R vs BP_R)
            [r,lags]=xcorr(EEG_R(:,chan)', BP_R', length(EEG_R),'normalized');
            plot(lags./2,r);xlabel('time [ms]');ylabel('Xcorr');
            title('EEG-R & BP-R','Color',condicolors(condi,:));
            xlim([-2e3 2e3]);ylim([-0.005 0.1]);
        subplot(4,6,6*3+condi); %  xcorr(EEG_R vs FB_R)
            [r,lags]=xcorr(EEG_R(:,chan)', FB_R', length(EEG_R),'normalized');
            plot(lags./2,r);xlabel('time [ms]');ylabel('Xcorr');
            title('EEG-R & FB-R','Color',condicolors(condi,:));
            xlim([-2e3 2e3]);ylim([-0.005 0.1]);
    end
    suptitle(['xcorr: EEG-BP/FB    EEG channel ' num2str(chan) '   subject ' num2str(seed)]);
        
figureName=['xcorrEEG-BPFBboolean EEG channel ' num2str(chan) ' subject ' num2str(seed)];
    % save the figure
    saveas(gcf,figureName,'jpg');
end

close all;

%% Xcorr (short-range): EEG vs EMG (equal length) (low correlation)
% [ xcorr(EEG_L vs EMG_L) + xcorr(EEG_L vs EMG_R) xcorr(EEG_R vs EMG_L) xcorr(EEG_R vs EMG_R) ]
% 6 conditions -subplots(2,6,condi) 
cd /ssd/zhibin/1overf/20220518_2P/Segmented_data/XcorrEEG_EMG

for chan=1:32
    figure('units','normalized','outerposition',[0 0 1 1]);
    for condi=1:6
        if condi==1; EEG_L=EEGCondi1L; EMG_L=EMGCondi1L; EEG_R=EEGCondi1R; EMG_R=EMGCondi1R; end
        if condi==2; EEG_L=EEGCondi2L; EMG_L=EMGCondi2L; EEG_R=EEGCondi2R; EMG_R=EMGCondi2R; end
        if condi==3; EEG_L=EEGCondi3L; EMG_L=EMGCondi3L; EEG_R=EEGCondi3R; EMG_R=EMGCondi3R; end
        if condi==4; EEG_L=EEGCondi4L; EMG_L=EMGCondi4L; EEG_R=EEGCondi4R; EMG_R=EMGCondi4R; end
        if condi==5; EEG_L=EEGCondi5L; EMG_L=EMGCondi5L; EEG_R=EEGCondi5R; EMG_R=EMGCondi5R; end
        if condi==6; EEG_L=EEGCondi6L; EMG_L=EMGCondi6L; EEG_R=EEGCondi6R; EMG_R=EMGCondi6R; end
        % compute
        subplot(4,6,condi); %  xcorr(EEG_L vs EMG_L)
            [r,lags]=xcorr(EEG_L(:,chan)', EMG_L', length(EMG_L),'normalized');
            plot(lags./2,r);xlabel('time [ms]');ylabel('Xcorr');
            title({['Condi ' num2str(condi)], 'EEG-L & EMG-L'},'Color',condicolors(condi,:));
            xlim([-2e3 2e3]);ylim([-0.1 0.1]);
        subplot(4,6,6+condi); %  xcorr(EEG_L vs EMG_R)
            minLength=min([size(EEG_L,1) length(EMG_R)]);
            [r,lags]=xcorr(EEG_L(1:minLength,chan)', EMG_R(1:minLength)', minLength,'normalized');
            plot(lags./2,r);xlabel('time [ms]');ylabel('Xcorr');
            title('EEG-L & EMG-R','Color',condicolors(condi,:));
            xlim([-2e3 2e3]);ylim([-0.1 0.1]);
        subplot(4,6,6*2+condi); %   xcorr(EEG_R vs EMG_L)
            minLength=min([size(EEG_L,1) length(EMG_R)]);
            [r,lags]=xcorr(EEG_R(1:minLength,chan)', EMG_L(1:minLength)', minLength,'normalized');
            plot(lags./2,r);xlabel('time [ms]');ylabel('Xcorr');
            title('EEG-R & EMG-L','Color',condicolors(condi,:));
            xlim([-2e3 2e3]);ylim([-0.1 0.1]);
        subplot(4,6,6*3+condi); %  xcorr(EEG_R vs FB_R)
            [r,lags]=xcorr(EEG_R(:,chan)', EMG_R', length(EMG_R),'normalized');
            plot(lags./2,r);xlabel('time [ms]');ylabel('Xcorr');
            title('EEG-R & EMG-R', 'Color',condicolors(condi,:));
            xlim([-2e3 2e3]);ylim([-0.1 0.1]);
    end
    suptitle(['xcorr: EEG-EMG    EEG channel ' num2str(chan) '   subject ' num2str(seed)]);
        
figureName=['xcorr EEG-EMG  EEG channel ' num2str(chan) ' subject ' num2str(seed)];
    % save the figure
    saveas(gcf,figureName,'jpg');
end

close all;
%%
