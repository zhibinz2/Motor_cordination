%% load data
clear;close all;
cd /ssd/zhibin/1overf/20220515_2P
cd /ssd/zhibin/1overf/20220517_2P
cd /ssd/zhibin/1overf/20220518_2P
%% P2 behaviral (short-range / local statistics; weak anticipation) on interval
% [ autocorr(BP_L) + xcorr(BP_L vs FB_L) + xcorr(BP_R vs FB_R) + autocorr(BP_R) ]
% 6 conditions -subplots(4,6,condi) 

figure;
for condi=1:6
    if condi==1; BP_L=BPCondi1L; FB_L=; BP_R=;FB_R=;end
    if condi==2; BP_L=BPCondi2L; FB_L=; BP_R=;FB_R=;end
    if condi==3; BP_L=BPCondi3L; FB_L=; BP_R=;FB_R=;end
    if condi==4; BP_L=BPCondi4L; FB_L=; BP_R=;FB_R=;end
    if condi==5; BP_L=BPCondi5L; FB_L=; BP_R=;FB_R=;end
    if condi==6; BP_L=BPCondi6L; FB_L=; BP_R=;FB_R=;end
    % compute
    subplot(4,6,condi); %  autocorr(BP_L)
        autocorr(Calinterval(BP_L'),500); xlabel('lags');ylabal('autocorrelation'); 
        title(['autocorr BP_L' '\nCondi ' num2str(condi)]);
    subplot(4,6,6+condi); %  xcorr(BP_L vs FB_L)
        time_series1=Calinterval(BP_L');% plot(time_series1);
        time_series2=Calinterval(FB_L');% plot(time_series2);
        TimeLength=min([length(time_series1) length(time_series2)]);
        [r,lags]=xcorr(time_series1(1:TimeLength), time_series2(1:TimeLength), 1000,'normalized');
        plot(lags./2,r);xlabel('time [ms]');ylabel('cross correlation');
        title(['xcorr BP_L & FB_L' '\nCondi ' num2str(condi)]);
    subplot(4,6,6*2+condi); %   xcorr(BP_R vs FB_R) 
        time_series1=Calinterval(BP_R');% plot(time_series1);
        time_series2=Calinterval(FB_R');% plot(time_series2);
        TimeLength=min([length(time_series1) length(time_series2)]);
        [r,lags]=xcorr(time_series1(1:TimeLength), time_series2(1:TimeLength), 1000,'normalized');
        plot(lags./2,r);xlabel('time [ms]');ylabel('cross correlation');
        title(['xcorr BP_R & FB_R' '\nCondi ' num2str(condi)]);
    subplot(4,6,6*3+condi); %  autocorr(BP_R)
        autocorr(Calinterval(BP_R'),500); xlabel('lags');ylabal('autocorrelation'); 
        title(['autocorr BP_R' '\nCondi ' num2str(condi)]);
    suptitle('short-range statistics: tapping intervals');
end
    
        

%% P2 behaviral (long-range / global statistics; strong anticipation)
% Pspectra(BP_L) + Pspectra(BP_R) + DFA(BP_L) + DFA(BP_R)
% 6 conditions -subplots(4,6,i)
maxfreq=5;sr=2000;
figure;
for condi=1:6
    if condi==1; BP_L=BPCondi1L; BP_R=;end
    if condi==2; BP_L=BPCondi2L; BP_R=;end
    if condi==3; BP_L=BPCondi3L; BP_R=;end
    if condi==4; BP_L=BPCondi4L; BP_R=;end
    if condi==5; BP_L=BPCondi5L; BP_R=;end
    if condi==6; BP_L=BPCondi6L; BP_R=;end
    subplot(4,6,condi); % Pspectra(BP_L)
        [fcoef1,fcoef2,cprod, freqs] = spectra(Calinterval(BP_L'),Calinterval(BP_L'),maxfreq,sr);
        plot(freqs,abs(cprod)); title('abs(cprod)');% cross spectra
        xlabel('freqs');ylabal('power'); 
        title(['Power spectra BP_L' '\nCondi ' num2str(condi)]);
        hold on;
        tbl=table(freqs',abs(cprod'));
        mdl=fitlm(tbl,'linear');
        plotAdded(mdl);
        
    subplot(4,6,6+condi);  % Pspectra(BP_R) 
        [fcoef1,fcoef2,cprod, freqs] = spectra(Calinterval(BP_L'),Calinterval(BP_L'),maxfreq,sr);
        plot(freqs,abs(cprod)); title('abs(cprod)');% cross spectra
        xlabel('freqs');ylabal('power'); 
        title(['Power spectra BP_R' '\nCondi ' num2str(condi)]);
        hold on;
        tbl=table(freqs',abs(cprod'));
        mdl=fitlm(tbl,'linear');
        plotAdded(mdl);
        
    subplot(4,6,6*2+condi);  % DFA(BP_L)
        addpath 
        [D,Alpha1]=DFA_main(Calinterval(BP_L'));
        title('['DFA BP_L' '\nCondi ' num2str(condi)]');
        
    subplot(4,6,6*2+condi);  % DFA(BP_R)
        [D,Alpha1]=DFA_main(Calinterval(BP_R'));
        title('['DFA BP_R' '\nCondi ' num2str(condi)]');
      
    suptitle('long-range statistics: tapping intervals');
end

%% P2 Error (short & long -range statistics)
% Error + autocorr (Error) + Pspectra(Error)+ DFA(Error)
% uncoupled condition (bidiretional error) -subplots(4,10,1-2)
% unidrectional conditions -subplots(4,10,3-4)
% bidirectional condition (bidiretional error) -subplots(4,10,5-10)

for condi=1:10
    if condi==1; BP_L=BPCondi1L; BP_R=; time_series1=find(BP_L==1); time_series2=find(BP_R==1); end
    if condi==2; BP_L=BPCondi2L; BP_R=; time_series1=find(BP_R==1); time_series2=find(BP_L==1); end
    if condi==3; BP_L=BPCondi3L; BP_R=; time_series1=find(BP_L==1); time_series2=find(BP_R==1); end
    if condi==4; BP_L=BPCondi4L; BP_R=;end
    if condi==5; BP_L=BPCondi5L; BP_R=;end
    if condi==6; BP_L=BPCondi6L; BP_R=;end
    if condi==7; BP_L=BPCondi6L; BP_R=;end
    if condi==8; BP_L=BPCondi6L; BP_R=;end
    if condi==9; BP_L=BPCondi6L; BP_R=;end
    if condi==10; BP_L=BPCondi6L; BP_R=;end
    % Compute error
    time_series1=find(BP_L==1);% plot(time_series1);
    time_series2=find(BP_R==1);% plot(time_series2);
    % TimeLength=min([length(time_series1) length(time_series2)]);
    % TimeLength=max([length(time_series1) length(time_series2)]);
    TimeLength=length(time_series1);
    Error=[];
    % find the time difference with the closest botton press 
    for i=1:TimeLength % i=StartStim2 % syncopation has one press less than condition 1 and 3
        [minValue,closetIndex]=min(abs(time_series2-time_series1(i))); % closetIndex in BottonPressTime
        Error(i)=time_series2(closetIndex)-time_series1(i);
    end

    subplot(4,10,condi); % Error base on L (uncoupled)
        plot(Error/2000,'r.');title('error');xlabel('taps');ylabel('timing error (s)');
        title(['error (BP_L ' '\nCondi ' num2str(condi)]);
        hold on;
        tbl=table(freqs',abs(cprod'));
        mdl=fitlm(tbl,'linear');
        plotAdded(mdl);
        
    subplot(4,10,10+condi);  % Error base on R (uncoupled)
        
    subplot(4,10,10*2+condi);  %  Error base on L (unidirection)
        
    subplot(4,10,10*3+condi);  %  Error base on R (unidirection)
        
    suptitle('short & long -range statistics: error');
end


% use a sliding window and do it again
% 6 conditions -subplots(3,6,i)


%% P2 EMG (short-range / local statistics; weak anticipation)
% [ autocorr(EMG_L) + xcorr (EMG_L vs EMG_R) + autocorr(EMG_R) ]
% 6 conditions - subplots(3,6,i)


%% P2 EMG (long-range / global statistics; strong anticipation)
% [ Pspectra(EMG_L) + Xspectra (EMG_L vs EMG_R) + Pspectra(EMG_R) ]
% 6 conditions - subplots(3,6,i)


%% EEG auto xcorr (weak anticipation, short-range / local statistics)
% [ autocorr(EEG_L) + xcorr(EEG_L vs EEG_R) + autocorr(EEG_R) ]
% each of the 6 conditions (subplots(3,6,i)) for each of the 32 channels;
% auto save the 32 figures in a folder named -


%% EEG Power spectra and Xspectra (strong anticipation, long-range / global statistics)
% [ Pspectra(EEG_L) + Xspectra(EEG_L vs EEG_R) + Pspectra(EEG_R) ]
% each of the 6 conditions (subplots(3,6,i)) for each of the 32 channels;
% auto save the 32 figures in a folder named - 





