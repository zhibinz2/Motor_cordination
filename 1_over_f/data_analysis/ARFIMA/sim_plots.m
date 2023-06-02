
%% PLOT-1 examine beta or alpha or the slope
freqs=[0.1:0.1:50];

figure;
alpha1=0.5;
P=1./(freqs.^alpha1);
plot(freqs,P,'g');% 1overf
xlabel('freqs');ylabel('1/f');
% xlim([0 10]);ylim([0 10])
hold on; 

alpha2=1;
P=1./(freqs.^alpha2);
plot(freqs,P,'b');% 1overf
xlabel('freqs');ylabel('1/f');
% xlim([0 10]);ylim([0 10])

alpha3=2;
P=1./(freqs.^alpha3);
plot(freqs,P,'r');% 1overf
xlabel('freqs');ylabel('1/f');
% xlim([0 4]);ylim([0 4])

alpha4=3;
P=1./(freqs.^alpha4);
plot(freqs,P,'m');% 1overf
xlabel('freqs');ylabel('1/f');
xlim([0 4]);ylim([0 4])

legend(['alpha1= ' num2str(alpha1)], ['alpha2= ' num2str(alpha2)],...
    ['alpha3= ' num2str(alpha3)], ['alpha4= ' num2str(alpha4)]);
%% 1.3Hz tapping data
cd /ssd/zhibin/1overf/20220609_2P/Segmented_data/Plots
% calculate std of the tapping intervals
cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing
%%  color scheme
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
ARFIMAcolors=[deepyellow;darkgreen;pink];
    
%% PLOT-2 Try the 2 functions fit the slope to find d
% Function 1
cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/ARFIMA
<<<<<<< HEAD
addpath(genpath('/home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/ARFIMA'));
addpath /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing
addpath /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/DFA
=======
addpath(genpath('/home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/ARFIMA'));
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/DFA

% XPS 17
addpath(genpath('C:\Users\zhouz\GitHub\Motor_cordination\1_over_f\data_analysis/ARFIMA'));
addpath C:\Users\zhouz\GitHub\Motor_cordination\1_over_f/data_analysis/2P_testing
addpath C:\Users\zhouz\GitHub\Motor_cordination\1_over_f/data_analysis/DFA

>>>>>>> ffc38334a8fbb643c76d8c2d62720c46ae129852
%%% ARFIMA(0,d,0)
figure('units','normalized','outerposition',[0 0 1 0.5]);
N=600;d=0.5;stdx=20;
[Z] = ARFIMA_SIM(N,[],[],d,stdx);
y=Z+750;
subplot(1,4,1); plot(y,'b');
subplot(1,4,2); autocorr(y,N-1);
subplot(1,4,3); Fs=1.3; [freqs,fcoef] = oneoverf(y,Fs);xlabel('Log(f)');ylabel('Log(power)');
subplot(1,4,4); [D,Alpha1]=DFA_main(y);

% Function 2
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/ARFIMA/dgp_arfima
% XPS 17
addpahh addpath C:\Users\zhouz\GitHub\Motor_cordination\1_over_f/ARFIMA/dgp_arfima
% ARFIMA(0,0.8,0) >>>>  Results=dgp_arfima(0,[],[],600,1,0.8);
figure('units','normalized','outerposition',[0 0 1 0.5]);
T=600;SD=20;d=0.5;F=0;
Z=dgp_arfima(0,[],[],T,SD,d,F);
y=Z+750;
subplot(1,4,1); plot(y,'b');
subplot(1,4,2); autocorr(y,T-1);
subplot(1,4,3); Fs=1.3; [freqs,fcoef] = oneoverf(y,Fs);xlabel('Log(f)');ylabel('Log(power)');
subplot(1,4,4); [D,Alpha1]=DFA_main(y);

%% PLOT-2-1 examine the 2 function code
clear
N=600;F=[];O=[];d=0.5;stdx=20;er=[];
[Z] = ARFIMA_SIM(N,F,O,d,stdx);

clear
c=0;ARo=[];MAo=[];T=100;SD=20;d=0.7;F=1;
sim=dgp_arfima(c,ARo,MAo,T,SD,d,F);

%% PLOT-3 Different d
figure('units','normalized','outerposition',[0 0 1 1]);
N=1000;stdx=20;
% ds=[-0.5 -0.2 0.5];
ds=[0.2 0.5 1];
for i=1:length(ds)
    [Z] = ARFIMA_SIM(N,[],[],ds(i),stdx);
    y=Z+750;
    subplot(length(ds),4,(i-1)*4+1); 
        plot(y,'b');title(['d = ' num2str(ds(i))]);
        ax = gca;
        ax.FontSize = 15; 
        xlabel('Taps');ylabel('Tapping interval (ms)')
        % ylim([680 840]);
    subplot(length(ds),4,(i-1)*4+2); 
        autocorr(y,N-1);title(['d = ' num2str(ds(i))]);
        ylabel('Autocorrelation (\rho)')
        ax = gca;
        ax.FontSize = 15; 
        % ylim([-0.5 1]);
    subplot(length(ds),4,(i-1)*4+3); 
        Fs=1.3; [freqs,fcoef,beta,xx,yy,FitValues] = oneoverf(y,Fs);
        plot(xx, yy,'bx');xlabel('Log10(f)');ylabel('Log10(power)');
        hold on;
        plot(xx',FitValues,'r--');
%         plotAdded(mdl);
        ax = gca;
        ax.FontSize = 15; 
        legend({'Data',['Fit (Beta=' num2str(beta) ')']},'Location','northeast');
        title(['d = ' num2str(ds(i)) '   (beta = d*2)']);
        % ylim([-2.5 1]);
    subplot(length(ds),4,(i-1)*4+4); 
        [D,Alpha1,n,F_n,FitValues]=DFA_main(y);
        plot(log10(n),log10(F_n),'bx');
        ax = gca;
        ax.FontSize = 15; 
        xlabel('Scale [log10(n)]'); ylabel('RMS [log10(F-n)]');
        hold on;
        plot(log10(n(1:end)),FitValues,'r--');
        legend({'Data',['Fit (H=' num2str(Alpha1) ')']},'Location','southeast');
        title(['d = ' num2str(ds(i)) '   (H = d+0.5)']);
        % ylim([0 3]);
end
set(gcf,'color','w'); % set background white for copying in ubuntu
% cd /usr/local/MATLAB/R2019a/toolbox/bioinfo/biodemos/ %suptitle.m
% cd /usr/local/MATLAB/R2022a/toolbox/bioinfo/biodemos/ % suptitle.m
% suptitle(['Different d with same std = ' num2str(stdx)]);
sgtitle(['Different d with same std = ' num2str(stdx)]);

figureName=['Different_d'];
saveas(gcf,figureName,'fig');

%% PLOT-3-1 add parcorr
figure('units','normalized','outerposition',[0 0 1 1]);
N=100;stdx=20;
ds=[0.2 0.5 1];
for i=1:length(ds)
    [Z] = ARFIMA_SIM(N,[],[],ds(i),stdx);
    y=Z+750;
    subplot(length(ds),5,(i-1)*5+1); 
        plot(y,'b');title(['d = ' num2str(ds(i))]);
        % ylim([680 840]);
    subplot(length(ds),5,(i-1)*5+2); 
        autocorr(y,N-1);title(['d = ' num2str(ds(i))]);
        ylabel('autocorr')
        % ylim([-0.5 1]);
    subplot(length(ds),5,(i-1)*5+3); 
        parcorr(y,N-1);title(['d = ' num2str(ds(i))]);
        ylabel('parcorr')
        % ylim([-0.5 1]);   
    subplot(length(ds),5,(i-1)*5+4); 
        Fs=1.3; [freqs,fcoef] = oneoverf(y,Fs);
        xlabel('Log(f)');ylabel('Log(power)');
        title(['d = ' num2str(ds(i)) '   (beta = d*2)']);
        % ylim([-2.5 1]);
    subplot(length(ds),5,(i-1)*5+5); 
        [D,Alpha1]=DFA_main(y);
        title(['d = ' num2str(ds(i)) '   (H = d+0.5)']);
        % ylim([0 3]);
end
% cd /usr/local/MATLAB/R2019a/toolbox/bioinfo/biodemos/ %suptitle.m
% cd /usr/local/MATLAB/R2022a/toolbox/bioinfo/biodemos/ % suptitle.m
suptitle(['Different d with same std = ' num2str(stdx)]);

%% PLOT-4 Different length
Ns=[100 500 1000 3000];d=1;stdx=20;
figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:length(Ns)
    [Z] = ARFIMA_SIM(Ns(i),[],[],d,stdx);
    y=Z+750;
    subplot(length(Ns),4,(i-1)*4+1); 
        plot(y,'b');
        title(['d = ' num2str(d) '    data length = ' num2str(Ns(i))]);
        % ylim([680 840]);
    subplot(length(Ns),4,(i-1)*4+2); 
        autocorr(y,Ns(i)-1);
        title(['d = ' num2str(d) '    data length = ' num2str(Ns(i))]);
        % ylim([-0.5 1]);
    subplot(length(Ns),4,(i-1)*4+3); 
        Fs=1.3; [freqs,fcoef] = oneoverf(y,Fs);
        xlabel('Log10(f)');ylabel('Log10(power)');
        title(['d = ' num2str(d) '    data length = ' num2str(Ns(i))]);
        % ylim([-2.5 1]);
    subplot(length(Ns),4,(i-1)*4+4); 
        [D,Alpha1]=DFA_main(y);
        title(['d = ' num2str(d) '    data length = ' num2str(Ns(i))]);
        % ylim([0 3]);
end
suptitle(['Different length with same d = ' num2str(d)]);

figureName=['Different_N'];
saveas(gcf,figureName,'fig');

%% PLOT-5 Segmenting 3000 taps into 100 taps/chuncks
Ns=[3000 1000 500 100];d=1;stdx=20;
figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:length(Ns)
    [Z] = ARFIMA_SIM(3000,[],[],d,stdx);
    y=Z+750;
    ys=y(1:Ns(i));
    subplot(length(Ns),4,(i-1)*4+1); 
        plot(ys,'b');
        title(['d = ' num2str(d) '    data length = ' num2str(Ns(i))]);
        % ylim([680 840]);
    subplot(length(Ns),4,(i-1)*4+2); 
        autocorr(ys,Ns(i)-1);
        title(['d = ' num2str(d) '    data length = ' num2str(Ns(i))]);
        % ylim([-0.5 1]);
    subplot(length(Ns),4,(i-1)*4+3); 
        Fs=1.3; [freqs,fcoef] = oneoverf(ys,Fs);
        xlabel('Log10(f)');ylabel('Log10(power)');
        title(['d = ' num2str(d) '    data length = ' num2str(Ns(i))]);
        % ylim([-2.5 1]);
    subplot(length(Ns),4,(i-1)*4+4); 
        [D,Alpha1]=DFA_main(ys);
        title(['d = ' num2str(d) '    data length = ' num2str(Ns(i))]);
        % ylim([0 3]);
end
suptitle(['The 1st Segment with different length from the 3000 taps']);

figureName=['Different_Segment'];
saveas(gcf,figureName,'fig');

%% PLOT-6-calculation for 100 or 600 taps see how stable the estimations are
%{ 
N=100; ds=[0.1:0.1:1]; repeat=100; stdx=20;
%}

% function [Betass,BetassErr,Hss,HssErr] = IterARFIMA(N,ds,repeat,stdx);
% 
% % calcualtion
% Betass=[];Hss=[];BetassErr=[];HssErr=[];
% for j=1:length(ds)
%     d=ds(j);
%     Betas=[];Hs=[];
%     for i=1:repeat 
%         % [Z] = ARFIMA_SIM(N,[],[],d,stdx);
%         Z=dgp_arfima(0,[],[],N,stdx,d,0);
%         y=Z+750;
%         % PSA
%         Fs=1.3; [~,~,beta] = oneoverf(y,Fs);
%         % DFA
%         [~,H]=DFA_main(y);
%         % append the values
%         Betas = [Betas beta];
%         Hs = [Hs H];  
%     end
% Betass=[Betass;Betas];
% BetassErr=[BetassErr; std(Betas)/sqrt(N)];
% Hss=[Hss;Hs];
% HssErr=[HssErr; std(Hs)/sqrt(N)];
% end
% end

%% PLOT-6-1 plot all 10 subplots
figure('units','normalized','outerposition',[0 0 1 1]);
for j=1:length(ds)
    subplot(2,10,j);    
    hist(Betass(j,:),repeat);
    xlabel('beta');ylabel('repeats');
    xlim([-2 5]);
    title(['d = ' num2str(ds(j))]);
    
    subplot(2,10,10+j);
    hist(Hss(j,:),repeat); 
    xlabel('H');ylabel('repeats');
    xlim([0 2]);
    title(['d = ' num2str(ds(j))]);
end
% suptitle(['Beta and H with different d repeating 100 times (ARFIMA-SIM ([],d,[]); length of ' num2str(N) ')']);
suptitle(['Beta and H with different d repeating 100 times (dgp-arfima ([],d,[]); length of ' num2str(N) ')']);

%% PLOT-6-2 : 2 plots with error bar
figure;
subplot(1,2,1);
errorbar(ds,mean(Betass,2),BetassErr,'.')
xlabel('d');ylabel('Beta');xlim([0 1.1]);
title(['Beta with standard error']);
B=2.*ds; hold on; plot(ds,B,'color',purple);
legend({'B','B=2*d'},'Location','southeast');
subplot(1,2,2);
errorbar(ds,mean(Hss,2),HssErr,'.')
xlabel('d');ylabel('H');xlim([0 1.1]);
title(['H with standard error']);
H=ds+0.5; hold on; plot(ds,H,'color',purple);
legend({'H','H=d+0.5'},'Location','southeast');
suptitle(['Beta and H with different d repeating 100 times (ARFIMA ([],d,[]); length of ' num2str(N) ')']);

%% PLOT-6-2-1 : 2 plots with error bar for 3 different length
addpath C:\Users\zhouz\GitHub\Motor_cordination\1_over_f\data_analysis\MSE-VARFI
Ns=[100 600 1000]; ds=[0.1:0.1:1]; repeat=10; stdx=20;
figure;
tic;
for i=1:3
    [Betass,BetassErr,Hss,HssErr] = IterARFIMA(Ns(i),ds,repeat,stdx);
    subplot(1,2,1);
        hold on;
        errorbar(ds,mean(Betass,2),BetassErr,'.','color',ARFIMAcolors(i,:));
        xlabel('d');ylabel('alpha');xlim([0 1.1]);ylim([0 2.5])
        title(['alpha with standard error']);
        hold off;
        if i==3;
            B=2.*ds; hold on; plot(ds,B,'color',purple);hold off;
            legend({'N=100','N=600','N=1000','Alpha=2*d'},'Location','southeast');
        end
    subplot(1,2,2);
        hold on;
        errorbar(ds,mean(Hss,2),HssErr,'.','color',ARFIMAcolors(i,:))
        xlabel('d');ylabel('H');xlim([0 1.1]);ylim([0.5 1.6])
        title(['H with standard error']);
        hold off;
        if i==3;
            H=ds+0.5; hold on; plot(ds,H,'color',purple);hold off;
            legend({'N=100','N=600','N=1000','H=d+0.5'},'Location','southeast');
        end
end
% suptitle(['alpha and H with different d repeating 100 times (ARFIMA([],d,[]); length of ' num2str(Ns) ')']);
% sgtitle(['alpha and H with different d repeating 100 times (ARFIMA([],d,[]); length of ' num2str(Ns) ')']);
toc

%Elapsed time is 640.222034 seconds.
figureName=['InterARFIMA_gdp20'];
figureName=['InterARFIMA_SIM100'];
saveas(gcf,figureName,'fig');
    
%% PLOT-6-3 : 2 violin plots
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/violin
figure;
subplot(1,2,1);
[h,L,MX,MED]=violin(Betass');
ylabel('\beta','FontSize',14)
xticks([1:10]);
xticklabels(string(ds));
xlabel('d');
title(['Beta']);
subplot(1,2,2);
[h,L,MX,MED]=violin(Hss');
ylabel('H','FontSize',14)
xticks([1:10]);
xticklabels(string(ds));
xlabel('d');
suptitle(['Beta and H with different d repeating 100 times (ARFIMA ([],d,[]); length of ' num2str(N) ')']);

%% PLOT-7 simulation plot for the army 
N=100; d=0.7; repeat=100; stdx=20;

Z=dgp_arfima(0,[],[],N,stdx,d,0);
y1=round(Z'+750); 
Y1=[];
for i=1:length(y1)
    Y1=[Y1 1 zeros(1,y1(i))];
end

Z=dgp_arfima(0,[],[],N,stdx,d,0);
y2=round(Z'+750);
Y2=[];
for i=1:length(y2)
    Y2=[Y2 1 zerxcos(1,y2(i))];
end
Y2=[zeros(1,750/2) Y2]; % add the offset
% figure;plot(y1,'r');hold on;plot(y2,'b');

figure('units','normalized','outerposition',[0 0 0.6 1]);
subplot(4,1,1);
plot(Y1,'r');hold on;plot(Y2,'b');
xlabel('time (ms)');
yticks(1);yticklabels('button pressed');
xlim([750*0 750*40]); ylim([0 1.5]);
legend('Player L','Player R');
title('Syncopated tapping sequences of 2 players');

subplot(4,2,3);
MinLength=min([length(y1) length(y2)]);
[acf,lags,bounds] = autocorr(y1(1:MinLength),MinLength-1);
stem(lags,acf,'r'); xlabel('Lag');ylabel('\rho(k)');
hold on;
h = line(lags,bounds(1)*ones(length(acf),1));
h1 = line(lags,bounds(2)*ones(length(acf),1));
set(h,'color',[1 0 1]);
set(h1,'color',[1 0 1]);
title('Autocorr: Player L');
subplot(4,2,5);
[acf,lags,bounds] = autocorr(y2(1:MinLength),MinLength-1);
stem(lags,acf,'b'); xlabel('Lag');ylabel('\rho(k)');
hold on;
h = line(lags,bounds(1)*ones(length(acf),1));
h1 = line(lags,bounds(2)*ones(length(acf),1));
set(h,'color',[1 0 1]);
set(h1,'color',[1 0 1]);
title('Autocorr: Player R');
subplot(4,2,7);
yy1=y1-mean(y1);yy2=y2-mean(y2);
[r,lags]=xcorr(yy1(1:MinLength),yy2(1:MinLength),MinLength-1,'normalized');
plot(lags,r,'color',[0 0.5 0]);xlabel('Lag');ylabel('\rho(k)');yline(0,'color',[1 0.8 0.2]);
title('Xcorr: Player L & R');

subplot(4,2,[4 6 8]);
[D,Alpha1,n,F_n,FitValues]=DFA_main(y1);
plot(log10(n),log10(F_n),'rx');
hold on; plot(log10(n),FitValues,'r');
[D,Alpha2,n,F_n,FitValues]=DFA_main(y2);
plot(log10(n),log10(F_n),'bx');
hold on; plot(log10(n),FitValues,'b');
legend({'Player L',['H=' num2str(Alpha1)],'Player R',['H=' num2str(Alpha2)],},'Location','southeast');
xlabel('Scale [log10(n)]') % win_lengths 
ylabel('RMS [log10(F-n)]') % RMS values
title('DFA')

figureName=['ZBPlot20220624'];
saveas(gcf,figureName,'jpg');
    
%% PLOT-7-1 simulation plot for the army 
%parameters
N=100; 
d=0.7; 
nrepeat=10; 
stdx=10; 
meanx = 750;

acfmean1 = zeros(1,N/2+1);
acfmean2 = zeros(1,N/2+1);
rmean12 = zeros(1,N+1);
boundsmean1 = zeros(1,2);
boundsmean2 = zeros(1,2);

for j = 1:nrepeat 

    Z=dgp_arfima(0,[],[],N,stdx,d,0);
    y1=round(Z'+meanx); 
    
    Z=dgp_arfima(0,[],[],N,stdx,d,0);
    y2=round(Z'+meanx);
    
    Y1=[];
    for i=1:length(y1)
        Y1 =[Y1 1 zeros(1,y1(i))];
    end
    Y2=[];
    for i=1:length(y2)
        Y2 =[Y2 1 zeros(1,y2(i))];
    end
    Y2 =[zeros(1,750/2) Y2]; % add the offset
    
    MinLength=min([length(y1) length(y2)]);
    
    y1 = y1-mean(y1);
    y2 = y2-mean(y2);
    
    [acf1,lags1,bounds1] = autocorr(y1(1:MinLength),N/2);
    acfmean1 = acfmean1+ acf1;
    boundsmean1 = boundsmean1 + bounds1;
    [acf2,lags2,bounds2] = autocorr(y2(1:MinLength),N/2);
    acfmean2 = acfmean2 + acf2;
    [r12,lags12]=xcorr(y1(1:MinLength),y2(1:MinLength),N/2,'normalized');
    rmean12 = rmean12 + r12;
    boundsmean2 = boundsmean2 + bounds2;
end;

acfmean1 = acfmean1/nrepeat;
acfmean2 = acfmean2/nrepeat;
rmean12 = rmean12/nrepeat;
boundsmean1=boundsmean1/nrepeat;
boundsmean2=boundsmean2/nrepeat;

% Plotting is afterwards
figure('units','normalized','outerposition',[0 0 0.6 1]);

subplot(4,1,1);
plot(Y1,'r');hold on;plot(Y2,'b');
xlabel('time (ms)');
yticks(1);yticklabels('button pressed');
xlim([750*0 750*40]); ylim([0 1.5]);
legend('Player L','Player R');
title('Syncopated tapping sequences of 2 players');

subplot(4,2,3);
stem(lags1,acfmean1,'r'); xlabel('Lag');ylabel('\rho(k)');
% plot(lags1,acfmean1,'r'); xlabel('Lag');ylabel('\rho(k)');
axis([0 50 -0.2 1]) 
hold on;
h = line(lags1,boundsmean1(1)*ones(length(acfmean1),1));
h1 = line(lags1,boundsmean1(2)*ones(length(acfmean1),1));
set(h,'color',[1 0 1]);
set(h1,'color',[1 0 1]);
title('Autocorr: Player L');

subplot(4,2,5);
stem(lags2,acfmean2,'b'); xlabel('Lag');ylabel('\rho(k)');
% plot(lags2,acfmean2,'b');xlabel('Lag');ylabel('\rho(k)');
axis([0 50 -0.2 1]) 
hold on;
h = line(lags2,boundsmean2(1)*ones(length(acfmean2),1));
h1 = line(lags2,boundsmean2(2)*ones(length(acfmean2),1));
set(h,'color',[1 0 1]);
set(h1,'color',[1 0 1]);
title('Autocorr: Player R');

subplot(4,2,7);
plot(lags12,rmean12,'color',[0 0.5 0]);xlabel('Lag');ylabel('\rho(k)');
yline(0,'color',[1 0.8 0.2]);
axis([-50 50 -0.2 1]);
title('Xcorr: Player L & R');

subplot(4,2,[4 6 8]);
[D,Alpha1,n,F_n,FitValues]=DFA_main(y1);
plot(log10(n),log10(F_n),'rx');
hold on; plot(log10(n),FitValues,'r');
[D,Alpha2,n,F_n,FitValues]=DFA_main(y2);
plot(log10(n),log10(F_n),'bx');
hold on; plot(log10(n),FitValues,'b');
legend({'Player L',['H=' num2str(Alpha1)],'Player R',['H=' num2str(Alpha2)],},'Location','southeast');
xlabel('Scale [log10(n)]') % win_lengths 
ylabel('RMS [log10(F-n)]') % RMS values
title('DFA')

figureName=['RMPlot20220625'];
saveas(gcf,figureName,'jpg');

%% PLOT-7-2 syncopation plot for the army 
cd D:\360MoveData\Users\alienware\Documents\GitHub\Motor_cordination\1_over_f\data_analysis\ARFIMA\Sim_plots
addpath D:\360MoveData\Users\alienware\Documents\GitHub\Motor_cordination\1_over_f\data_analysis\2P_testing
cd D:\360MoveData\Users\alienware\Documents\GitHub\Motor_cordination\1_over_f\data_analysis\2P_testing
addpath D:\360MoveData\Users\alienware\Documents\GitHub\Motor_cordination\1_over_f\data_analysis\DFA
cd G:\1overf\20220610_2P
cd G:\1overf\20220610_2P\Segmented_data\1_50Hz_ICAautomized
clear
load('BP20220610.mat')
load('FB20220610.mat')

condi=1; BP_L=BPCondi1L; FB_L=FBCondi1L; BP_R=BPCondi1R;FB_R=FBCondi1R;
condi=2; BP_L=BPCondi2L; FB_L=FBCondi2L; BP_R=BPCondi2R;FB_R=FBCondi2R;
condi=3; BP_L=BPCondi3L; FB_L=FBCondi3L; BP_R=BPCondi3R;FB_R=FBCondi3R;
condi=4; BP_L=BPCondi4L; FB_L=FBCondi4L; BP_R=BPCondi4R;FB_R=FBCondi4R;
% syncopation expt
conditionNames={'uncoupled' 'L-lead' 'R-lead' 'mutual-1.3Hz'};
% cut off from the last bottom press
BP_L=cutofflast(BP_L);
BP_R=cutofflast(BP_R);

y1=Calinterval(BP_L')./2;
y2=Calinterval(BP_R')./2;

% time (ms) for 100 taps
Inter100=round(100/1.3,3);

sr=2000;
nrepeat=5;
y1=[];y2=[];
MinLength=74; % choose an even number
acfmean1=zeros(1,round(MinLength/2)+1);
boundsmean1=zeros(1,2);
acfmean2=zeros(1,round(MinLength/2)+1);
boundsmean2=zeros(1,2);
rmean12=zeros(1,round(MinLength)+1);


for i=1:nrepeat
    y1(i).BPint=Calinterval(BP_L((Inter100*sr*(i-1)+1):(Inter100*sr*i))')./2;
    y2(i).BPint=Calinterval(BP_R((Inter100*sr*(i-1)+1):(Inter100*sr*i))')./2;
    
    % MinLength=min([length(y1(i).BPint) length(y2(i).BPint)]);
    
    y1(i).BPint = y1(i).BPint-mean(y1(i).BPint);
    y2(i).BPint = y2(i).BPint-mean(y2(i).BPint);
    
    [acf1,lags1,bounds1] = autocorr(y1(i).BPint(1:MinLength),round(MinLength/2));
    acfmean1 = acfmean1+ acf1;
    boundsmean1 = boundsmean1 + bounds1;
    [acf2,lags2,bounds2] = autocorr(y2(i).BPint(1:MinLength),round(MinLength/2));
    acfmean2 = acfmean2 + acf2;
    boundsmean2 = boundsmean2 + bounds2;
    [r12,lags12]=xcorr(y1(i).BPint(1:MinLength),y2(i).BPint(1:MinLength),round(MinLength/2),'normalized');
    rmean12 = rmean12 + r12;
    
end;

acfmean1 = acfmean1/nrepeat;
acfmean2 = acfmean2/nrepeat;
rmean12 = rmean12/nrepeat;
boundsmean1=boundsmean1/nrepeat;
boundsmean2=boundsmean2/nrepeat;


% calculate H in the 5 window
Alphas1=[];Alphas2=[];% Hs1Err=[];Hs2Err=[];  [HssErr; std(Hs)/sqrt(N)];
for i=1:nrepeat
    [D,Alpha1,n,F_n,FitValues]=DFA_main(y1(i).BPint(1:MinLength));
    Alphas1=[Alphas1 Alpha1];
    % Hs1Err=[Hs1Err std(Alphas1)/sqrt(N)];
    [D,Alpha2,n,F_n,FitValues]=DFA_main(y2(i).BPint(1:MinLength));
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
xlim([sr*0 sr*Inter100/3]); ylim([0 1.5]);
legend('Player L','Player R');
title(['Syncopated tapping sequences of 2 players' conditionNames(condi)]);

subplot(4,2,3);
stem(lags1,acfmean1,'r'); xlabel('Lag');ylabel('\rho(k)');
% plot(lags1,acfmean1,'r'); xlabel('Lag');ylabel('\rho(k)');
axis([0 MinLength/2 -0.3 1]) 
hold on;
h = line(lags1,boundsmean1(1)*ones(length(acfmean1),1));
h1 = line(lags1,boundsmean1(2)*ones(length(acfmean1),1));
set(h,'color',[1 0 1]);
set(h1,'color',[1 0 1]);
title('Autocorr: Player L');

subplot(4,2,5);
stem(lags2,acfmean2,'b'); xlabel('Lag');ylabel('\rho(k)');
% plot(lags2,acfmean2,'b');xlabel('Lag');ylabel('\rho(k)');
axis([0 MinLength/2 -0.3 1]) 
hold on;
h = line(lags2,boundsmean2(1)*ones(length(acfmean2),1));
h1 = line(lags2,boundsmean2(2)*ones(length(acfmean2),1));
set(h,'color',[1 0 1]);
set(h1,'color',[1 0 1]);
title('Autocorr: Player R');

subplot(4,2,7);
plot(lags12,rmean12,'color',[0 0.5 0]);xlabel('Lag');ylabel('\rho(k)');
yline(0,'color',[1 0.8 0.2]);xline(0,'color',[1 0.8 0.2]);
axis([-MinLength/2 MinLength/2 -0.3 1]);
title('Xcorr: Player L & R');

subplot(4,2,[4 6]);
i=2;
[D,Alpha1,n,F_n,FitValues]=DFA_main(y1(i).BPint(1:MinLength));
plot(log10(n),log10(F_n),'rx');
hold on; plot(log10(n),FitValues,'r');
[D,Alpha2,n,F_n,FitValues]=DFA_main(y2(i).BPint(1:MinLength));
plot(log10(n),log10(F_n),'bx');
hold on; plot(log10(n),FitValues,'b');
legend({'Player L',['H=' num2str(Alpha1)],'Player R',['H=' num2str(Alpha2)],},'Location','southeast');
xlabel('Scale [log10(n)]') % win_lengths 
ylabel('RMS [log10(F-n)]') % RMS values
title('DFA')

subplot(4,2,[8]);
plot(1:nrepeat,Alphas1,'r');hold on; plot(1:nrepeat,Alphas2,'b')
xlabel('window');ylabel('H');xlim([0 6]);
xticks(1:5);xticklabels(string(1:5));
title(['H values in 5 consecutive window']);
legend({'Player L','Player R'},'Location','southeast');

% figureName=['RMPlot20220625'];
% saveas(gcf,figureName,'jpg');

