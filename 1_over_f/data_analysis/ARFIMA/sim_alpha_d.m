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

%% PLOT-2 Try the 2 functions fit the slope to find d
% Function 1
cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/ARFIMA
addpath(genpath('/home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/ARFIMA'));
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/DFA
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
% ARFIMA(0,0.8,0) >>>>  Results=dgp_arfima(0,[],[],600,1,0.8);
figure('units','normalized','outerposition',[0 0 1 0.5]);
T=600;SD=20;d=0.5;F=0;
Z=dgp_arfima(0,[],[],T,SD,d,F);
y=Z+750;
subplot(1,4,1); plot(y,'b');
subplot(1,4,2); autocorr(y,T-1);
subplot(1,4,3); Fs=1.3; [freqs,fcoef] = oneoverf(y,Fs);xlabel('Log(f)');ylabel('Log(power)');
subplot(1,4,4); [D,Alpha1]=DFA_main(y);

%% Different d
figure('units','normalized','outerposition',[0 0 1 1]);
N=1000;stdx=20;
% ds=[-0.5 -0.2 0.5];
ds=[0.2 0.5 1];
for i=1:length(ds)
    [Z] = ARFIMA_SIM(N,[],[],ds(i),stdx);
    y=Z+750;
    subplot(length(ds),4,(i-1)*4+1); 
        plot(y,'b');title(['d = ' num2str(ds(i))]);
        % ylim([680 840]);
    subplot(length(ds),4,(i-1)*4+2); 
        autocorr(y,N-1);title(['d = ' num2str(ds(i))]);
        % ylim([-0.5 1]);
    subplot(length(ds),4,(i-1)*4+3); 
        Fs=1.3; [freqs,fcoef,beta,xx,yy,FitValues] = oneoverf(y,Fs);
        plot(xx, yy,'bx');xlabel('Log10(f)');ylabel('Log10(power)');
        hold on;
        plot(xx',FitValues,'r--');
%         plotAdded(mdl);
        legend({'Data',['Fit (Beta=' num2str(beta) ')']},'Location','northeast');
        title(['d = ' num2str(ds(i)) '   (beta = d*2)']);
        % ylim([-2.5 1]);
    subplot(length(ds),4,(i-1)*4+4); 
        [D,Alpha1,n,F_n,FitValues]=DFA_main(y);
        plot(log10(n),log10(F_n),'bx');
        xlabel('Scale [log10(n)]'); ylabel('RMS [log10(F-n)]');
        hold on;
        plot(log10(n(1:end)),FitValues,'r--');
        legend({'Data',['Fit (H=' num2str(Alpha1) ')']},'Location','southeast');
        title(['d = ' num2str(ds(i)) '   (H = d+0.5)']);
        % ylim([0 3]);
end
% cd /usr/local/MATLAB/R2019a/toolbox/bioinfo/biodemos/ %suptitle.m
% cd /usr/local/MATLAB/R2022a/toolbox/bioinfo/biodemos/ % suptitle.m
suptitle(['Different d with same std = ' num2str(stdx)]);

figureName=['Different_d'];
saveas(gcf,figureName,'fig');

%% add parcorr
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

%% Different length
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

%% Segmenting 3000 taps into 100 taps/chuncks
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

%% PLOT for 100 or 600 taps see how stable the estimations are
N=100; ds=[0.1:0.1:1]; repeat=100; stdx=20;

% calcualtion
Betass=[];Hss=[];BetassErr=[];HssErr=[];
for j=1:length(ds)
    d=ds(j);
    Betas=[];Hs=[];
    for i=1:repeat 
        % [Z] = ARFIMA_SIM(N,[],[],d,stdx);
        Z=dgp_arfima(0,[],[],N,stdx,d,0);
        y=Z+750;
        % PSA
        Fs=1.3; [~,~,beta] = oneoverf(y,Fs);
        % DFA
        [~,H]=DFA_main(y);
        % append the values
        Betas = [Betas beta];
        Hs = [Hs H];  
    end
Betass=[Betass;Betas];
BetassErr=[BetassErr; std(Betas)/sqrt(N)];
Hss=[Hss;Hs];
HssErr=[HssErr; std(Hs)/sqrt(N)];
end

% plot all 10 subplots
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
% 2 plots with error bar
figure;
subplot(1,2,1);
errorbar(ds,mean(Betass,2),BetassErr,'.')
xlabel('d');ylabel('Beta');xlim([0 1.1]);
title(['Beta with standard error']);
subplot(1,2,2);
errorbar(ds,mean(Hss,2),HssErr,'.')
xlabel('d');ylabel('H');xlim([0 1.1]);
title(['H with standard error']);
suptitle(['Beta and H with different d repeating 100 times (ARFIMA ([],d,[]); length of ' num2str(N) ')']);

% 2 violin plots
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

%% simulation plot for the army 
N=100; d=1; repeat=100; stdx=20;

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
    Y2=[Y2 1 zeros(1,y2(i))];
end
Y2=[zeros(1,750/2) Y2];
% figure;plot(y1,'r');hold on;plot(y2,'b');

figure;
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
% % PSA
% Fs=1.3; [~,~,beta] = oneoverf(y,Fs);
% DFA
[D,Alpha1,n,F_n,FitValues]=DFA_main(y1);
plot(log10(n),log10(F_n),'rx');
hold on; plot(log10(n),FitValues,'r--');
[D,Alpha2,n,F_n,FitValues]=DFA_main(y2);
plot(log10(n),log10(F_n),'bx');
hold on; plot(log10(n),FitValues,'b--');
legend({'Player L',['H=' num2str(Alpha1)],'Player R',['H=' num2str(Alpha2)],},'Location','southeast');
xlabel('Scale [log10(n)]') % win_lengths 
ylabel('RMS [log10(F-n)]') % RMS values
title('DFA')

figureName=['ZBPlot20220624'];
saveas(gcf,figureName,'jpg');
    


