%% PLOT-7-2 sliding 100 taps for all
cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/ARFIMA/Sim_plots
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing
cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/DFA
cd /ssd/zhibin/1overf/20220610_2P/Segmented_data/1_50Hz_ICAautomized

cd D:\360MoveData\Users\alienware\Documents\GitHub\Motor_cordination\1_over_f\data_analysis\ARFIMA\Sim_plots
addpath D:\360MoveData\Users\alienware\Documents\GitHub\Motor_cordination\1_over_f\data_analysis\2P_testing
cd D:\360MoveData\Users\alienware\Documents\GitHub\Motor_cordination\1_over_f\data_analysis\2P_testing
addpath D:\360MoveData\Users\alienware\Documents\GitHub\Motor_cordination\1_over_f\data_analysis\DFA
cd G:\1overf\20220610_2P
cd G:\1overf\20220610_2P\Segmented_data\1_50Hz_ICAautomized


clear
load('BP20220610.mat')
load('FB20220610.mat')

cd D:\360MoveData\Users\alienware\Documents\GitHub\Motor_cordination\1_over_f\data_analysis\ARFIMA\Sim_plots\20220610

% syncopation expt
conditionNames={'uncoupled' 'L-lead' 'R-lead' 'mutual-1.3Hz'};
conditionNames={'uncoupled' 'L-lead' 'R-lead' 'mutual-2Hz'};

% number of samples for 100 taps
sr=2000;
hz=2; % hz=1.3; % tapping frequency
Inter100=round(100/hz,3)*sr;

for condi=1:4
    
if condi==1; BP_L=BPCondi1L; FB_L=FBCondi1L; BP_R=BPCondi1R;FB_R=FBCondi1R;end
if condi==2; BP_L=BPCondi2L; FB_L=FBCondi2L; BP_R=BPCondi2R;FB_R=FBCondi2R;end
if condi==3; BP_L=BPCondi3L; FB_L=FBCondi3L; BP_R=BPCondi3R;FB_R=FBCondi3R;end
if condi==4; BP_L=BPCondi4L; FB_L=FBCondi4L; BP_R=BPCondi4R;FB_R=FBCondi4R;end

% cut off from the last bottom press
BP_L=cutofflast(BP_L);
BP_R=cutofflast(BP_R);

% total number of taps
y1=Calinterval(BP_L')./2;
y2=Calinterval(BP_R')./2;

% minimal number of sliding window (every 0.2*Inter100)
increment=round(0.2*Inter100);
nrepeat=(min([length(BP_L) length(BP_R)])-Inter100)/increment;
nrepeat=floor(nrepeat);

% initialize some variables
y1=[];y2=[];
lags=50;
acfmean1=zeros(1,lags+1);
length1=[];
boundsmean1=zeros(1,2);
acfmean2=zeros(1,lags+1);
boundsmean2=zeros(1,2);
rmean12=zeros(1,2*lags+1);

for i=1:nrepeat
    
    y1(i).BPint=Calinterval(BP_L((increment*(i-1)+1):((increment*(i-1)+1+Inter100-1)))')./2;
    y2(i).BPint=Calinterval(BP_R((increment*(i-1)+1):((increment*(i-1)+1+Inter100-1)))')./2;
    
    MinLength=min([length(y1(i).BPint) length(y2(i).BPint)]);
    
    y1(i).BPint = y1(i).BPint-mean(y1(i).BPint);
    y2(i).BPint = y2(i).BPint-mean(y2(i).BPint);
    
    [acf1,lags1,bounds1] = autocorr(y1(i).BPint(1:MinLength),lags);
    acfmean1 = acfmean1+ acf1;
    boundsmean1 = boundsmean1 + bounds1;
    [acf2,lags2,bounds2] = autocorr(y2(i).BPint(1:MinLength),lags);
    acfmean2 = acfmean2+ acf2;
    boundsmean2 = boundsmean2 + bounds2;
    [r12,lags12]=xcorr(y1(i).BPint(1:MinLength),y2(i).BPint(1:MinLength),lags,'normalized');
    rmean12 = rmean12 + r12;
    
end;

acfmean1 = acfmean1/nrepeat;
acfmean2 = acfmean2/nrepeat;
rmean12 = rmean12/nrepeat;
boundsmean1=boundsmean1/nrepeat;
boundsmean2=boundsmean2/nrepeat;

% calculate H in the 5 window
Alphas1=[];Alphas2=[];% Hs1Err=[];Hs2Err=[];  [HssErr; std(Hs)/sqrt(N)];

% find the minimal length
minLen1=[];minLen2=[];
for i=1:nrepeat
    minLen1=[minLen1 size(y1(i).BPint,2)];
    minLen2=[minLen2 size(y2(i).BPint,2)];
end
minLen=min([minLen1 minLen2]);  
    
for i=1:nrepeat
    [D,Alpha1,n,F_n,FitValues]=DFA_main(y1(i).BPint(1:minLen));
    Alphas1=[Alphas1 Alpha1];
    % Hs1Err=[Hs1Err std(Alphas1)/sqrt(N)];
    [D,Alpha2,n,F_n,FitValues]=DFA_main(y2(i).BPint(1:minLen));
    Alphas2=[Alphas2 Alpha2];
end

% figure;
% plot(1:5,Alphas1,'r');hold on; plot(1:5,Alphas2,'b')
% xlabel('window');ylabel('H');xlim([0 6]);
% xticks(1:5);xticklabels(string(1:5));
% title(['H values in 5 consecutive window']);
% legend({'Player L','Player R'},'Location','southeast');

% Plotting is afterwards
figure('units','normalized','outerposition',[0 0 0.6 1]);

subplot(4,1,1);
plot(BP_L,'r');hold on;plot(BP_R,'b');
xlabel('time (ms)');
yticks(1);yticklabels('button pressed');
xlim([sr*0 Inter100/3]); ylim([0 1.5]);
legend('Player L','Player R');
title(['Syncopated tapping sequences of 2 players' conditionNames(condi)]);

subplot(4,2,3);
stem(lags1,acfmean1,'r'); xlabel('Lag');ylabel('\rho(k)');
% plot(lags1,acfmean1,'r'); xlabel('Lag');ylabel('\rho(k)');
axis([0 lags -0.3 1]) 
hold on;
h = line(lags1,boundsmean1(1)*ones(length(acfmean1),1));
h1 = line(lags1,boundsmean1(2)*ones(length(acfmean1),1));
set(h,'color',[1 0 1]);
set(h1,'color',[1 0 1]);
title('Autocorr: Player L');

subplot(4,2,5);
stem(lags2,acfmean2,'b'); xlabel('Lag');ylabel('\rho(k)');
% plot(lags2,acfmean2,'b');xlabel('Lag');ylabel('\rho(k)');
axis([0 lags -0.3 1]) 
hold on;
h = line(lags2,boundsmean2(1)*ones(length(acfmean2),1));
h1 = line(lags2,boundsmean2(2)*ones(length(acfmean2),1));
set(h,'color',[1 0 1]);
set(h1,'color',[1 0 1]);
title('Autocorr: Player R');

subplot(4,2,7);
plot(lags12,rmean12,'color',[0 0.5 0]);xlabel('Lag');ylabel('\rho(k)');
yline(0,'color',[1 0.8 0.2]);xline(0,'color',[1 0.8 0.2]);
axis([-lags lags -0.3 1]);
title('Xcorr: Player L & R');

subplot(4,2,[4 6]);
i=floor(nrepeat/2);
[D,Alpha1,n,F_n,FitValues]=DFA_main(y1(i).BPint(1:minLen));
plot(log10(n),log10(F_n),'rx');
hold on; plot(log10(n),FitValues,'r');
[D,Alpha2,n,F_n,FitValues]=DFA_main(y2(i).BPint(1:minLen));
plot(log10(n),log10(F_n),'bx');
hold on; plot(log10(n),FitValues,'b');
legend({'Player L',['H=' num2str(Alpha1)],'Player R',['H=' num2str(Alpha2)],},'Location','southeast');
xlabel('Scale [log10(n)]') % win_lengths 
ylabel('RMS [log10(F-n)]') % RMS values
title(['DFA in window ' num2str(i)]);

subplot(4,2,[8]);
cla;
plot(1:nrepeat,Alphas1,'r');hold on; plot(1:nrepeat,Alphas2,'b')
xlabel('sliding windows');ylabel('H');xlim([0 nrepeat+1]);
xticks(1:nrepeat);xticklabels(string(1:nrepeat));
% title(['H values in 5 consecutive window']);
legend({'Player L','Player R'},'Location','southeast');

% cd /ssd/zhibin/1overf/20220610_2P/Segmented_data/Plots
% figureName=['condi_' num2str(condi)];
% saveas(gcf,figureName,'fig');

end

%% PLOT Spetra and DFA for EEG
cd /ssd/zhibin/1overf/20220610_2P/Segmented_data/1_50Hz_ICAautomized
clear
load('EEG20220610.mat')

% syncopation expt
conditionNames={'uncoupled' 'L-lead' 'R-lead' 'mutual-1.3Hz'};
% conditionNames={'uncoupled' 'L-lead' 'R-lead' 'mutual-2Hz'};


tempEEG = load('/ssd/zhibin/1overf/20220610_2P/Segmented_data/1_50Hz_ICAautomized/EEG20220610.mat'); 
% tempBP = load('/ssd/zhibin/1overf/20220610_2P/Segmented_data/1_50Hz_ICAautomized/BP20220610.mat'); 

data(1,1).EEG = tempEEG.EEGCondi1L;
data(1,2).EEG = tempEEG.EEGCondi1R;
data(2,1).EEG = tempEEG.EEGCondi2L;
data(2,2).EEG = tempEEG.EEGCondi2R;
data(3,1).EEG = tempEEG.EEGCondi3L;
data(3,2).EEG = tempEEG.EEGCondi3R;
data(4,1).EEG = tempEEG.EEGCondi4L;
data(4,2).EEG = tempEEG.EEGCondi4R;

% number of samples for 100 taps
sr=2000;
hz=2; % hz=1.3; % tapping frequency
Inter100=round(100/hz,3)*sr;

for condi=1:4
    if condi==1; EEG_L=EEGCondi1L; EEG_R=EEGCondi1R; end 
    if condi==2; EEG_L=EEGCondi2L; EEG_R=EEGCondi2R; end 
    if condi==3; EEG_L=EEGCondi3L; EEG_R=EEGCondi3R; end 
    if condi==4; EEG_L=EEGCondi4L; EEG_R=EEGCondi4R; end 
    
    Fs=2000;chan=15;
    
    % minimal number of sliding window (every 0.2*Inter100)
    increment=round(0.2*Inter100);
    nrepeat=(min([size(EEG_L,1) size(EEG_R,1)])-Inter100)/increment;
    nrepeat=floor(nrepeat);
    
    for i=1:nrepeat % i=15
    % segment into intervals
    y1(i).EEGint=zscore(EEG_L((increment*(i-1)+1):((increment*(i-1)+1+Inter100-1)),chan));
    y2(i).EEGint=zscore(EEG_R((increment*(i-1)+1):((increment*(i-1)+1+Inter100-1)),chan));
    
    % [freqs,fcoef,beta,xx,yy,FitValues] = oneoverf(y1(i).EEGint,Fs);
    y=y1(i).EEGint;
    [freqs,fcoef,beta1,xx,yy,FitValues] = oneoverf(y1(i).EEGint,Fs);
    plot(xx, yy,'bx');xlabel('Log10(f)');ylabel('Log10(power)');
    hold on;plot(xx',FitValues,'r--');
    
    % fft on the fly
    tempf = fft(y,[],1);
    for m = 1:50 % 50 frequencies
        freqs = (m-1)*10+1:m*10; % downsample freqs to only 500 samples
        power = sum(abs(tempf(freqs,:)),1);
        pow(j).EEG(m,:,k) = power;
    end
    
    % use allspetra
    
    
