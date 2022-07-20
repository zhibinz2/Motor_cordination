%% PLOT 1: Auto&Xcorr and DFA and Granger Causality for each trial

cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/ARFIMA/Sim_plots
cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/DFA
cd /ssd/zhibin/1overf/20220713_2P/Segmented_data
load('data_variables20220713.mat')

% d estimate and removal
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/MSE-VARFI
% Granger Causality
addpath(genpath('/home/zhibin/Documents/GitHub/MVGC1'));
cd /home/zhibin/Documents/GitHub/MVGC1
run startup
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/VAR_Granger

% Calculation
y1=[];y2=[];y12=[];XcorrR=[];
nTrials=12;
for i=1:nTrials
    y1(i).BPint=Calinterval((BP(1).BP{i})')./sr; % in second
    y2(i).BPint=Calinterval((BP(2).BP{i})')./sr; % in second
    
    MinLength=min([length(y1(i).BPint) length(y2(i).BPint)]);
%     if mod(MinLength,2) == 1; % if odd number, make it even 
%         MinLength = MinLength-1;
%     end
    
    % remove the mean
    y1(i).BPint = y1(i).BPint-mean(y1(i).BPint);
    y2(i).BPint = y2(i).BPint-mean(y2(i).BPint);
    
    % Autocorr of y1
    acf=[];lags=[];bounds=[];
    [acf,lags,bounds] = autocorr(y1(i).BPint(1:MinLength),MinLength-1);
    y1(i).acf = acf;
    y1(i).bounds = bounds;
    y1(i).lags = lags;
    % Autocorr of y2
    acf=[];lags=[];bounds=[];
    [acf,lags,bounds] = autocorr(y2(i).BPint(1:MinLength),MinLength-1);
    y2(i).acf = acf;
    y2(i).bounds = bounds;
    y2(i).lags = lags;
    % Xcorr between y1 and y2
    r12=[];lags12=[];
    [r12,lags12]=xcorr(y1(i).BPint(1:MinLength),y2(i).BPint(1:MinLength),30,'normalized');
    y12(i).r12 = r12;
    y12(i).lags12 = lags12;
    % calculate H (or alpha) in DFA of y1 and y2
    n=[];F_n=[];FitValues=[];
    [~,Alpha,n,F_n,FitValues]=DFA_main(y1(i).BPint(1:MinLength));
    y1(i).Alpha=Alpha;
    y1(i).n=n;
    y1(i).F_n=F_n;
    y1(i).FitValues=FitValues;
    n=[];F_n=[];FitValues=[];
    [~,Alpha,n,F_n,FitValues]=DFA_main(y2(i).BPint(1:MinLength));
    y2(i).Alpha=Alpha;
    y2(i).n=n;
    y2(i).F_n=F_n;
    y2(i).FitValues=FitValues;   

    % Granger causality
    F=[];pval=[];,sig=[];
    %[F,pval,sig] = mystatespace(y1(i).BPint,y2(i).BPint); % Statespace method
    [F,pval,sig] = myautocov(y1(i).BPint,y2(i).BPint); % Autocov method
    y12(i).F=F;
    y12(i).pval=pval;
    y12(i).sig=sig;
end

% plotting
cd /ssd/zhibin/1overf/20220713_2P/Segmented_data/Plots/Corr_DFA

for i=1:nTrials
    
    figure('units','normalized','outerposition',[0 0 1 1]);

    subplot(4,1,1);
    plot([1:length(BP(1).BP{i})]./sr,BP(1).BP{i},'r');
    hold on;
    plot([1:length(BP(2).BP{i})]./sr,BP(2).BP{i},'b');
    xlabel('time (s)');
    yticks(1); yticklabels('button pressed');ylim([0 1.5]);
    hold off;
    legend('Player L','Player R');
    title(['Tapping sequences of 2 players in Trial - ' num2str(i) ': Condition - ' conditionNames{allPerm(i)} ]);

    subplot(4,2,3);
    stem(y1(i).lags,y1(i).acf,'r'); xlabel('Lag');ylabel('\rho(k)');
    ylim([-0.3 1]);% axis([0 MinLength/2 -0.3 1]) 
    hold on;
    h = line(y1(i).lags,(y1(i).bounds)'*ones(size(y1(i).lags)));
    set(h,'color',[1 0 1]);
    title('Autocorr: Player L');

    subplot(4,2,5);
    stem(y2(i).lags,y2(i).acf,'b'); xlabel('Lag');ylabel('\rho(k)');
    ylim([-0.3 1]);% axis([0 MinLength/2 -0.3 1]) 
    hold on;
    h = line(y2(i).lags,(y2(i).bounds)'*ones(size(y2(i).lags)));
    set(h,'color',[1 0 1]);
    title('Autocorr: Player R');

    subplot(4,2,7);
    plot(y12(i).lags12,y12(i).r12,'color',[0 0.5 0]);xlabel('Lag');ylabel('\rho(k)');
    yline(0,'color',[1 0.8 0.2]);xline(0,'color',[1 0.8 0.2]);
    ylim([-0.3 1]); % axis([-MinLength/2 MinLength/2 -0.3 1]);
    title('Xcorr: Player L & R');

    subplot(4,2,[4 6]);
    plot(log10(y1(i).n),log10(y1(i).F_n),'rx');
    hold on; plot(log10(y1(i).n),y1(i).FitValues,'r');
    plot(log10(y2(i).n),log10(y2(i).F_n),'bx');
    hold on; plot(log10(y2(i).n),y2(i).FitValues,'b');
    % axis([0.6 2 -1.8 0]);
    legend({'Player L',['H=' num2str(y1(i).Alpha)],'Player R',['H=' num2str(y2(i).Alpha)],},'Location','southeast');
    xlabel('Scale [log10(n)]') % win_lengths 
    ylabel('RMS [log10(F-n)]') % RMS values
    title('DFA')
    
    subplot(4,6,[22]);
    plot_pw(y12(i).F);
    title('Pairwise-conditional GC');
    % xticks([1 2]); xticklabels([{'L'} {'R'}]);% tickangle(90);
    set(gca,'XTick',1:2);
    set(gca,'XTickLabel',[{'L'} {'R'}]);
    set(gca,'YTick',1:2);
    set(gca,'YTickLabel',[{'L'} {'R'}]);
    
    subplot(4,6,[23]);
    plot_pw(y12(i).pval);
    title(['p-values (' tstat '-test)']);
    set(gca,'XTick',1:2);
    set(gca,'XTickLabel',[{'L'} {'R'}]);
    set(gca,'YTick',1:2);
    set(gca,'YTickLabel',[{'L'} {'R'}]);
    
    subplot(4,6,[24]);
    plot_pw(y12(i).sig);
    title(['Significant at \alpha = ' num2str(alpha)]);
    set(gca,'XTick',1:2);
    set(gca,'XTickLabel',[{'L'} {'R'}]);
    set(gca,'YTick',1:2);
    set(gca,'YTickLabel',[{'L'} {'R'}]);

%     figureName=['Corr_DFA-trial-' num2str(i)];
%     saveas(gcf,figureName,'jpg');

end
close all;

%% PLOT 2 