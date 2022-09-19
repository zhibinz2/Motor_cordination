%% PLOT 1: Auto&Xcorr and DFA and Granger Causality for each trial

cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/ARFIMA/Sim_plots
cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/DFA

cd /ssd/zhibin/1overf/20220713_2P/Segmented_data
cd /ssd/zhibin/1overf/20220721_2P/Segmented_data

load('data_variables20220713.mat')
load('data_variables20220721.mat')

% d estimate and removal
% http://www.lucafaes.net/LMSE-MSE_VARFI.html
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
    % [F,pval,sig] = mystatespace(y1(i).BPint,y2(i).BPint); % Statespace method
    [F,pval,sig] = myautocov(y1(i).BPint,y2(i).BPint); % Autocov method
    y12(i).F=F;
    y12(i).pval=pval;
    y12(i).sig=sig;
end

%% PLOT 1:  plotting
cd /ssd/zhibin/1overf/20220713_2P/Segmented_data/Plots/Corr_DFA
cd /ssd/zhibin/1overf/20220804_2P/Segmented_data/Plots/Corr_DFA

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
    stem(y12(i).lags12,y12(i).r12,'color',[0 0.5 0]);xlabel('Lag');ylabel('\rho(k)');
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

    figureName=['Corr_DFA-trial-' num2str(i)];
    saveas(gcf,figureName,'jpg');

end
close all;

%% PLOT 2: Auto&Xcorr and DFA and Autocov Granger Causality and PSA for each trial
cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/ARFIMA/Sim_plots
cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/DFA

cd /ssd/zhibin/1overf/20220713_2P/Segmented_data
cd /ssd/zhibin/1overf/20220721_2P/Segmented_data
cd /ssd/zhibin/1overf/20220804_2P/Segmented_data
load('data_variables20220713.mat')
load('data_variables20220721.mat')

% d estimate and removal
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/MSE-VARFI
% Granger Causality
addpath(genpath('/home/zhibin/Documents/GitHub/MVGC1'));
% rmpath('/home/zhibin/Documents/GitHub/MVGC1'); % remove duplecate startup file from path

% Structure viewer not showing properly after running the following 2 lines
% No need to run the following 2 lines
% cd /home/zhibin/Documents/GitHub/MVGC1
% run startup

% remove this path to recover structure viewing
% rmpath(genpath('/home/zhibin/Documents/GitHub/MVGC1')); 

addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/VAR_Granger

% parameters for powers
% open /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing/spectra2
maxfreq=50; EEGwin=10; % 10 Second with df of 0.1Hz

% Calculation on BP
y1=[];y2=[];y12=[];XcorrR=[];
nTrials=12;
alpha = 0.05; 
d_min=-0.5; d_max=1;
for i=1:nTrials
    % extract BP01
    y1(1,i).BP01=(BP(1).BP{i})'; % in boolean
    y2(1,i).BP01=(BP(2).BP{i})'; % in boolean 
   
    % extract BP intervals
    y1(1,i).BPint=Calinterval(y1(1,i).BP01)./sr; % in second
    y2(1,i).BPint=Calinterval(y2(1,i).BP01)./sr; % in second
    % MinLength=min([length(y1(i).BPint) length(y2(i).BPint)]);

    % Resample method 
    Maxlength=max([length(y1(1,i).BPint) length(y2(1,i).BPint)]);
    y1(1,i).resam_BP_inter=(resample(y1(1,i).BPint,Maxlength,length(y1(1,i).BPint)));
    y2(1,i).resam_BP_inter=(resample(y2(1,i).BPint,Maxlength,length(y2(1,i).BPint)));
        
    % remove the mean
    y1(1,i).resam_BP_inter = y1(i).resam_BP_inter-mean(y1(i).resam_BP_inter);
    y2(1,i).resam_BP_inter = y2(i).resam_BP_inter-mean(y2(i).resam_BP_inter);

    % Autocorr of y1
    acf=[];lags=[];bounds=[];
    [acf,lags,bounds] = autocorr(y1(i).resam_BP_inter,20);
    y1(i).acf = acf;
    y1(i).bounds = bounds;
    y1(i).lags = lags;
    
    % Autocorr of y2
    acf=[];lags=[];bounds=[];
    [acf,lags,bounds] = autocorr(y2(i).resam_BP_inter,20);
    y2(i).acf = acf;
    y2(i).bounds = bounds;
    y2(i).lags = lags;
    
    % Xcorr between y1 and y2
    r12=[];lags12=[];
    % d estimate
    est_d1=d_estimation(y1(1,i).resam_BP_inter,d_min,d_max);
    est_d2=d_estimation(y2(1,i).resam_BP_inter,d_min,d_max);
    %d filter
    [Ytmp1,~]=remove_d(y1(1,i).resam_BP_inter',est_d1); %filtered data
    [Ytmp2,~]=remove_d(y2(1,i).resam_BP_inter',est_d2); %filtered data
    % Xcorr
    [r12,lags12]=xcorr(Ytmp1',Ytmp2',10,'normalized');
    y12(i).r12 = r12;
    y12(i).lags12 = lags12;
    
    % calculate H (or alpha) in DFA of y1 and y2
    n=[];F_n=[];FitValues=[];
    [~,Alpha,n,F_n,FitValues]=DFA_main(y1(i).resam_BP_inter);
    y1(i).Alpha=Alpha;
    y1(i).n=n;
    y1(i).F_n=F_n;
    y1(i).FitValues=FitValues;
    n=[];F_n=[];FitValues=[];
    [~,Alpha,n,F_n,FitValues]=DFA_main(y2(i).resam_BP_inter);
    y2(i).Alpha=Alpha;
    y2(i).n=n;
    y2(i).F_n=F_n;
    y2(i).FitValues=FitValues;   

    % Granger causality
    F=[];pval=[];,sig=[];
    %[F,pval,sig] = mystatespace(y1(i).BPint,y2(i).BPint); % Statespace method
    [F,pval,sig] = myautocov(Ytmp1,Ytmp2); % Autocov method
    y12(i).F=F;
    y12(i).pval=pval;
    y12(i).sig=sig;
end

% Calculate EEG Power of the whole trial
for i=1:nTrials
%     y1(1,i).EEG=zscore(EEG(1).EEG{i},[],1);
%     y2(1,i).EEG=zscore(EEG(2).EEG{i},[],1);
    % extract EEG in each condition and normalized with pacing
    y1(1,i).EEG=zscore(EEG(1).EEG{i}./repmat(mean(abs(EEG_pacing(1).EEG_pacing{i}),1),size(EEG(1).EEG{i},1),1),[],1);
    y2(1,i).EEG=zscore(EEG(2).EEG{i}./repmat(mean(abs(EEG_pacing(2).EEG_pacing{i}),1),size(EEG(2).EEG{i},1),1),[],1);
    
    % create a function that reorganize EEG into time x chan x chunks
    % feed it to allspectra (for 2 s chunks) or spectra3 (for 10s chunks)
    [~,freqs, ~, y1(i).pow, ~] = spectra3(y1(i).EEG,sr,maxfreq,EEGwin);
    [~,freqs, ~, y2(i).pow, ~] = spectra3(y2(i).EEG,sr,maxfreq,EEGwin);
    % later examine pow eppow corr cprod in imagesc and create myallspectra
end
    
%% PLOT 2: plotting
cd /ssd/zhibin/1overf/20220713_2P/Segmented_data/Plots/Corr_DFA_GC_PSA
cd /ssd/zhibin/1overf/20220721_2P/Segmented_data/Plots/Corr_DFA_GC_PSA
cd /ssd/zhibin/1overf/20220808_2P/Segmented_data/Plots/Corr_DFA_GC_PSA
cd /ssd/zhibin/1overf/20220811_2P/Segmented_data/Plots/Corr_DFA_GC_PSA
cd /ssd/zhibin/1overf/20220816_2P/Segmented_data/Plots/Corr_DFA_GC_PSA

% check condiSeq for synchronization
condiSeq=allPerm;

for i=1:nTrials
    
    % figure('units','normalized','outerposition',[0 0 1 1]);
    canvas(0.5,0.5);

    subplot(4,1,1);
    plot([1:length(BP(1).BP{i})]./sr,BP(1).BP{i},'r');
    hold on;
    plot([1:length(BP(2).BP{i})]./sr,BP(2).BP{i},'b');
    xlabel('time (s)');
    yticks(1); yticklabels('button pressed');ylim([0 1.5]);
    hold off;
    legend('Player L','Player R');
    title(['Tapping sequences of 2 players in Trial - ' num2str(i) ': Condition - ' conditionNames{condiSeq(i)}]);

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
    stem(y12(i).lags12,y12(i).r12,'color',[0 0.5 0]);xlabel('Lag');ylabel('\rho(k)');
    yline(0,'color',[1 0.8 0.2]);xline(0,'color',[1 0.8 0.2]);
    ylim([-0.3 1]); % axis([-MinLength/2 MinLength/2 -0.3 1]);
    title('Xcorr: Player L & R');

    subplot(4,6,[10 11]);
    plot(log10(y1(i).n),log10(y1(i).F_n),'rx');
    hold on; plot(log10(y1(i).n),y1(i).FitValues,'r');
    plot(log10(y2(i).n),log10(y2(i).F_n),'bx');
    hold on; plot(log10(y2(i).n),y2(i).FitValues,'b');
    % axis([0.6 2 -1.8 0]);
    legend({'Player L',['H=' num2str(y1(i).Alpha)],'Player R',['H=' num2str(y2(i).Alpha)],},'Location','southeast','NumColumns',2);
    xlabel('Scale [log10(n)]') % win_lengths 
    ylabel('RMS [log10(F-n)]') % RMS values
    title('DFA')
    
    subplot(4,6,[12]);
    b=bar([y12(i).F(2) y12(i).F(3)]);set(gca, 'XTickLabel', {'L -> R' 'R -> L'});
    b.FaceColor='flat';b.CData(1,:) = [1 0 0];b.CData(2,:) = [0 0 1];
    % axis([0.25 2.75 0 max([y12(i).F(2) y12(i).F(3)])+0.05]); 
    xlim([0.25 2.75]); ylim([0 2]);
    if length(find(y12(i).sig==1))==0;
        tempF=y12(i).F(find(y12(i).sig==0));
        tempP=round(y12(i).pval(find(y12(i).sig==0)),3);
        text((2-1-0.5), tempF(1)+0.1, ['P = ' num2str(tempP(1))],'Color',[1 0 0]);
        text((3-1-0.5), tempF(2)+0.1, ['P = ' num2str(tempP(2))],'Color',[0 0 1]);
    end
    if length(find(y12(i).sig==1))==1;
        text((find(y12(i).sig==1)-1-0.5), y12(i).F(find(y12(i).sig==1))+0.1, ['* P = ' num2str(round(y12(i).pval(find(y12(i).sig==1)),3))],'Color',[1 0 0]);
        text((find(y12(i).sig==0)-1-0.5), y12(i).F(find(y12(i).sig==0))+0.1, ['P = ' num2str(round(y12(i).pval(find(y12(i).sig==0)),3))],'Color',[0 0 1]);
    end
    if length(find(y12(i).sig==1))==2;
        tempF=y12(i).F(find(y12(i).sig==1));
        tempP=round(y12(i).pval(find(y12(i).sig==1)),3);
        text((2-1-0.5), tempF(1)+0.1, ['* P = ' num2str(tempP(1))],'Color',[1 0 0]);
        text((3-1-0.5), tempF(2)+0.1, ['* P = ' num2str(tempP(2))],'Color',[0 0 1]);
    end
    
    title('Pairwise-conditional GC');
    xlabel(['Significant at \alpha = ' num2str(alpha)]);

    subplot(4,12,[31:33 43:45]);
    imagesc(y1(i).pow);colorbar;caxis([0 1e6]);
    yticks([0:10:500]);yticklabels(string(0:1:50));
    title('Power');ylabel('Frequencies (Hz)');xlabel('Channels (L player)');
    subplot(4,12,[34:36 46:48]);
    imagesc(y2(i).pow);colorbar;caxis([0 1e6]);
    yticks([0:10:500]);yticklabels(string(0:1:50));
    title('Power');xlabel('Channels (R player)'); %ylabel('Frequencies');
    
    figureName=['Corr_DFA_GC_PSA-trial-' num2str(i)];
    saveas(gcf,figureName,'jpg');
    close all;
end

%% PLOT 2-1: plotting
cd /ssd/zhibin/1overf/20220804_2P/Segmented_data/Plots/Plot2_1/
mkdir Plot2_1
cd Plot2_1
% Plots - color scheme
red   = [1 0 0];
pink  = [1 0.65 0.75];
blue  = [0 0 1];
mediumblue = [0 0.4 0.7];
green = [0 1 0];
darkgreen = [0 0.5 0];
grey  = [0.5 0.5 0.5];
yellow  = [1 1 0];
deepyellow  = [1 0.8 0.2];
megenta = [1 0 1];% fill([0 1 1 0],[0 0 1 1],megenta)
cyan = [0 1 1]; % fill([0 1 1 0],[0 0 1 1],cc)
purple = [0.6 0.1 0.9];
condicolors=[darkgreen;red;blue;megenta;purple;purple];
HNLcolors = [darkgreen; deepyellow; pink];
% % test color
% showcolor=pink;
% imagesc(cat(3,showcolor(1),showcolor(2),showcolor(3)));

% check condiSeq for synchronization
% condiSeq=allPerm;

nTrials=12;
alpha = 0.05;


% Autocorr L R and DFA, 3X12 subplots
% figure('units','normalized','outerposition',[0 0 1 1]);
canvas(0.5,0.3);
for i=1:nTrials
    subplot(3,12,i); % Autocorr L
    stem(y1(i).lags,y1(i).acf,'r'); xlabel('Lag');ylabel('\rho(k)');
    ylim([-0.3 1]);% axis([0 MinLength/2 -0.3 1]) 
    hold on;
    h = line(y1(i).lags,(y1(i).bounds)'*ones(size(y1(i).lags)));
    set(h,'color',[1 0 1]);
    title([conditionNames{condiSeq(i)}],'Color',condicolors(condiSeq(i),:));
    subtitle('Autocorr: Player L','Color',condicolors(condiSeq(i),:));

    subplot(3,12,12+i); % Autocorr R
    stem(y2(i).lags,y2(i).acf,'b'); xlabel('Lag');ylabel('\rho(k)');
    ylim([-0.3 1]);% axis([0 MinLength/2 -0.3 1]) 
    hold on;
    h = line(y2(i).lags,(y2(i).bounds)'*ones(size(y2(i).lags)));
    set(h,'color',[1 0 1]);
    subtitle('Autocorr: Player R','Color',condicolors(condiSeq(i),:));

    subplot(3,12,24+i); % DFA
    plot(log10(y1(i).n),log10(y1(i).F_n),'rx');
    hold on; plot(log10(y1(i).n),y1(i).FitValues,'r');
    plot(log10(y2(i).n),log10(y2(i).F_n),'bx');
    hold on; plot(log10(y2(i).n),y2(i).FitValues,'b');
    % axis([0.6 2 -1.8 0]);
    % legend({'L',['H=' num2str(y1(i).Alpha)],'R',['H=' num2str(y2(i).Alpha)],},'Location','southeast','NumColumns',2);
    xlabel('Scale [log10(n)]') % win_lengths 
    ylabel('RMS [log10(F-n)]') % RMS values
    title('DFA','Color',condicolors(condiSeq(i),:));
    subtitle(['H_L=' num2str(round(y1(i).Alpha,2)) ' H_R=' num2str(round(y2(i).Alpha,2))],'Color',condicolors(condiSeq(i),:))
    ylim([-2 0.2]);
end
sgtitle(['Session - ' num2str(seed) ' ' sessionTypes{session}]);
figureName=['Autocorr_DFA'];
saveas(gcf,figureName,'fig');



% Xcorr & GC 2X12 subplots
canvas(0.5,0.2);
for i=1:nTrials
    subplot(2,12,i);
    stem(y12(i).lags12,y12(i).r12,'color',[0 0.5 0]);xlabel('Lag');ylabel('\rho(k)');
    yline(0,'color',[1 0.8 0.2]);xline(0,'color',[1 0.8 0.2]);
    ylim([-0.3 1]); % axis([-MinLength/2 MinLength/2 -0.3 1]);
    title([conditionNames{condiSeq(i)}],'Color',condicolors(condiSeq(i),:));
    subtitle(['Xcorr: L & R' ],'Color',condicolors(condiSeq(i),:));

    subplot(2,12,12+i);
    b=bar([y12(i).F(2) y12(i).F(3)]);set(gca, 'XTickLabel', {'L -> R' 'R -> L'});
    b.FaceColor='flat';b.CData(1,:) = [1 0 0];b.CData(2,:) = [0 0 1];
    % axis([0.25 2.75 0 max([y12(i).F(2) y12(i).F(3)])+0.05]); 
    xlim([0.25 2.75]); ylim([0 1.5]);
    if length(find(y12(i).sig==1))==0;
        tempF=y12(i).F(find(y12(i).sig==0));
        tempP=round(y12(i).pval(find(y12(i).sig==0)),3);
        text((2-1-0.5), tempF(1)+0.1, [num2str(round(tempP(1),2))],'Color',[1 0 0]);
        text((3-1-0.5), tempF(2)+0.1, [num2str(round(tempP(2),2))],'Color',[0 0 1]);
    end
    if length(find(y12(i).sig==1))==1;
        text((find(y12(i).sig==1)-1-0.5), y12(i).F(find(y12(i).sig==1))+0.1, ['*' num2str(round(y12(i).pval(find(y12(i).sig==1)),2))],'Color',[1 0 0]);
        text((find(y12(i).sig==0)-1-0.5), y12(i).F(find(y12(i).sig==0))+0.1, [num2str(round(y12(i).pval(find(y12(i).sig==0)),2))],'Color',[0 0 1]);
    end
    if length(find(y12(i).sig==1))==2;
        tempF=y12(i).F(find(y12(i).sig==1));
        tempP=round(y12(i).pval(find(y12(i).sig==1)),3);
        text((2-1-0.5), tempF(1)+0.1, ['*' num2str(round(tempP(1),2))],'Color',[1 0 0]);
        text((3-1-0.5), tempF(2)+0.1, ['*' num2str(round(tempP(2),2))],'Color',[0 0 1]);
    end
    title([conditionNames{condiSeq(i)}],'Color',condicolors(condiSeq(i),:));
    subtitle('GC','Color',condicolors(condiSeq(i),:));
    xlabel(['Significant at \alpha = ' num2str(alpha)]);
end
sgtitle(['Session - ' num2str(seed) ' ' sessionTypes{session}]);
figureName=['Xcorr_GC'];
saveas(gcf,figureName,'fig');


    
% EEG power spectra
canvas(0.55,0.15);
for i=1:nTrials
    subplot(2,12,i);
%     imagesc(y1(i).pow);colorbar;caxis([0 1e6]);colormap('jet');
%     yticks([0:100:500]);yticklabels(string(0:10:50));ylabel('Frequencies (Hz)');
    logimagesc((y1(i).pow),freqs);% caxis([0 0.8e7]);
    title([conditionNames{condiSeq(i)}],'Color',condicolors(condiSeq(i),:));
    subtitle('Power','Color',condicolors(condiSeq(i),:));
    xlabel('Chans (L)');

    subplot(2,12,12+i);
%     imagesc(y2(i).pow);colorbar;caxis([0 1e6]);colormap('jet');
%     yticks([0:10:500]);yticklabels(string(0:1:50));ylabel('Frequencies (Hz)');
    logimagesc((y2(i).pow),freqs);% caxis([0 0.8e7]);
    title([conditionNames{condiSeq(i)}],'Color',condicolors(condiSeq(i),:));
    subtitle('Power','Color',condicolors(condiSeq(i),:));
    xlabel('Chans (R)'); 
end
sgtitle(['Session - ' num2str(seed) ' ' sessionTypes{session}]);
figureName=['EEG_spectra'];
saveas(gcf,figureName,'fig');


%% PLOT 3: slinding window: DFA and Autocov Granger Causality for each trial
 
cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/ARFIMA/Sim_plots
cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/DFA

cd /ssd/zhibin/1overf/20220713_2P/Segmented_data
load('data_variables20220713.mat')

cd /ssd/zhibin/1overf/20220721_2P/Segmented_data
load('data_variables20220721.mat')


% d estimate and removal
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/MSE-VARFI
% Granger Causality
addpath(genpath('/home/zhibin/Documents/GitHub/MVGC1'));
% cd /home/zhibin/Documents/GitHub/MVGC1
% run startup
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/VAR_Granger

% remove this path to recover structure viewing
% rmpath(genpath('/home/zhibin/Documents/GitHub/MVGC1')); 

% parameters for powers
% open /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing/spectra2
maxfreq=50; win=10; % 10 Second with df of 0.1Hz

% Sliding win on BP01 for DFA and GC
y1;y2;
% Calculation of sliding window
% window of time (ms) for 100 taps
winsize = 80*sr; % 80 Second with df of 0.0125 Hz, about 104 taps
% winsize = 10*sr; % 10 Second with df of 0.1Hz
% winsize = 2*sr; % 2 Second with df of 0.5 Hz
overlapsize = round(winsize*0.3); % number of samples for each overlapping window of3.3 seconds

% Matric Structure (no heading in variable viewer)
ys1=[];ys2=[];ys12=[];
nTrials=12;
alpha = 0.05; 
for i=1:nTrials
    % BP01
%     y1(i).BP01; % in boolean
%     y2(i).BP01; % in boolean 
    BP(1).BP{i}; % in boolean
    BP(2).BP{i}; % in boolean 
    
    % loop through each sliding window
    % determine how many sliding windows needed for EEG
%     MinLength=min([length(BP(i).BP01) length(y2(i).BP01)]); % cut to same length
    MinLength=min([length(BP(1).BP{i}) length(BP(2).BP{i})]);
    nwin=floor((MinLength-winsize)/overlapsize);

     for k=1:nwin;
        samples = ((k-1)*overlapsize+1):(k*overlapsize+winsize);
        
        % extract BP intervals
%         temp1=Calinterval(y1(i).BP01(samples))./sr; % in second
%         temp2=Calinterval(y2(i).BP01(samples))./sr; % in second
        temp1=Calinterval((BP(1).BP{i})')./sr; % in second
        temp2=Calinterval((BP(2).BP{i})')./sr; % in second
        
        % remove the mean
        ys1(i,k).BPint = temp1-mean(temp1);
        ys2(i,k).BPint = temp2-mean(temp2);
        
        % cut to same length again
        MinLengthys=min([length(temp1) length(temp2)]);

        % calculate H (or alpha) in DFA of y1 and y2
        n=[];F_n=[];FitValues=[];
        [~,Alpha,n,F_n,FitValues]=DFA_main(ys1(i,k).BPint(1:MinLengthys));
        ys1(i,k).Alpha=Alpha;
        ys1(i,k).n=n;
        ys1(i,k).F_n=F_n;
        ys1(i,k).FitValues=FitValues;
        n=[];F_n=[];FitValues=[];
        [~,Alpha,n,F_n,FitValues]=DFA_main(ys2(i,k).BPint(1:MinLengthys));
        ys2(i,k).Alpha=Alpha;
        ys2(i,k).n=n;
        ys2(i,k).F_n=F_n;
        ys2(i,k).FitValues=FitValues;   

        % Granger causality
        F=[];pval=[];sig=[];
        %[F,pval,sig] = mystatespace(y1(i).BPint,y2(i).BPint); % Statespace method
        [F,pval,sig] = myautocov(ys1(i,k).BPint,ys2(i,k).BPint); % Autocov method
        ys12(i,k).F=F;
        ys12(i,k).pval=pval;
        ys12(i,k).sig=sig;
     end
end

% remove this path to recover structure viewing
rmpath(genpath('/home/zhibin/Documents/GitHub/MVGC1')); 
open ys1;
open ys12;
open ys2;
addpath(genpath('/home/zhibin/Documents/GitHub/MVGC1')); 

tic
% Heading structure (preferred)
ys1=[];ys2=[];ys12=[];
nTrials=12;
alpha = 0.05; 
for i=1:nTrials
    % BP01
    temp01=BP(1).BP{i}; % in boolean
    temp02=BP(2).BP{i}; % in boolean 
    
    % loop through each sliding window
    % determine how many sliding windows needed for EEG
    MinLength=min([length(temp01) length(temp02)]); % cut to same length
    nwin=floor((MinLength-winsize)/overlapsize);
    
    % EEG
    EEGtemp1=zscore(EEG(1).EEG{i},[],1);
    EEGtemp2=zscore(EEG(2).EEG{i},[],1);

     for k=1:nwin;
        samples = ((k-1)*overlapsize+1):(k*overlapsize+winsize);
        
        % extract BP intervals
        temp1=Calinterval(temp01(samples)')./sr; % in second
        temp2=Calinterval(temp02(samples)')./sr; % in second
        % remove the mean
        ys1(1,i).BPint{k} = temp1-mean(temp1);
        ys2(1,i).BPint{k} = temp2-mean(temp2);
        

        % extract EEG wins
        EEGtemp11=EEGtemp1(samples,:);
        EEGtemp22=EEGtemp2(samples,:);
        % remove the mean
        ys1(1,i).EEG{k} = EEGtemp11-repmat(mean(EEGtemp11,1),size(EEGtemp11,1),1);
        ys2(1,i).EEG{k} = EEGtemp22-repmat(mean(EEGtemp22,1),size(EEGtemp22,1),1);
        
        % cut to same length again
        MinLengthys=min([length(temp1) length(temp2)]);

        % calculate H (or alpha) in DFA of y1 and y2
        n=[];F_n=[];FitValues=[];
        [~,Alpha,n,F_n,FitValues]=DFA_main(ys1(i).BPint{k}(1:MinLengthys));
        ys1(1,i).Alpha{k}=Alpha;
        n=[];F_n=[];FitValues=[];
        [~,Alpha,n,F_n,FitValues]=DFA_main(ys2(i).BPint{k}(1:MinLengthys));
        ys2(1,i).Alpha{k}=Alpha;

        % Granger causality
        F=[];pval=[];sig=[];
        [F,pval,sig] = myautocov(ys1(i).BPint{k}(1:MinLengthys),ys2(i).BPint{k}(1:MinLengthys)); % Autocov method
        ys12(1,i).F{k}=F;
        ys12(1,i).pval{k}=pval;
        ys12(1,i).sig{k}=sig;

        % GC for EEG
        F=[];pval=[];sig=[];
        chan=5; % choose a channel
        [F,pval,sig] = myautocov(ys1(1,i).EEG{k}(:,chan),ys2(1,i).EEG{k}(:,chan)); % Autocov method
        ys12(1,i).EEG_F{k}=F;
        ys12(1,i).EEG_pval{k}=pval;
        ys12(1,i).EEG_sig{k}=sig;
     end
end
toc

% remove this path to recover structure viewing
rmpath(genpath('/home/zhibin/Documents/GitHub/MVGC1')); 
open ys1;
open ys12;
open ys2;
addpath(genpath('/home/zhibin/Documents/GitHub/MVGC1')); 



% sum EEG power
%  for m = 1:50 % 50 frequencies
%     freqs = (m-1)*10+1:m*10; % downsample freq
%     power_L = sum(abs(tempf_L(freqs,:)),1);
%     power_R = sum(abs(tempf_R(freqs,:)),1);
%     pow_L(j).EEG(m,:,k) = power_L; % freq x chan x time
%     pow_R(j).EEG(m,:,k) = power_R;
%  end

    
%% PLOT 3: plotting
cd /ssd/zhibin/1overf/20220713_2P/Segmented_data/Plots/Corr_DFA_GC_PSA
cd /ssd/zhibin/1overf/20220721_2P/Segmented_data/Plots/

for i=1:nTrials
    
    figure('units','normalized','outerposition',[0 0 1 1]);

    for k=1:nwin;

    end

    figureName=['Corr_DFA_GC_SPA-trial-' num2str(i)];
    saveas(gcf,figureName,'jpg');
    close all;
end

%% PLOT 4: 500ms before BP
% Calculation on BP
y1=[];y2=[];y12=[];
nTrials=12;
BPwin=20;
EEGwin=0.5;
for i=1:nTrials
    % extract BP01
    y1(1,i).BP01=(BP(1).BP{i})'; % in boolean
    y2(1,i).BP01=(BP(2).BP{i})'; % in boolean 

    % extract the indices for time point of each BP 
    y1(1,i).BPind=find(y1(1,i).BP01); % BPind=y1(1,i).BPind; plot(BPind,ones(1,length(BPind)),'bo');
    y2(1,i).BPind=find(y2(1,i).BP01); % BPind=y2(1,i).BPind; plot(BPind,ones(1,length(BPind)),'bo');

    % count number of taps
    y1(1,i).numtaps=size(y1(i).BPind,2);
    y2(1,i).numtaps=size(y2(i).BPind,2);
    % y12(1,i).MinTaps=min([y1(1,i).numtaps y2(1,i).numtaps]);
    % y12(1,i).MaxTaps=max([y1(1,i).numtaps y2(1,i).numtaps]);

    % open P2_Syncopate_plots.m
    % cd /ssd/zhibin/1overf/20220610_2P/Segmented_data/Plots
    % samples=1:min([length(y1(1,i).BP01) length(y2(1,i).BP01)]);
    % open sync.m % corr method?
    % y12(1,i).CorrBPint=corr(Calinterval(y1(1,i).BP01(samples)),Calinterval(y2(1,i).BP01(samples)));

    % convert to intervals (length=numtaps-1)
    y1(1,i).BPint=Calinterval(y1(1,i).BP01)'; 
    y2(1,i).BPint=Calinterval(y2(1,i).BP01)'; 

    % Resample first then smooth (length=Maxlength of numtaps-1)
    Maxlength=max([size(y1(1,i).BPint,1) size(y2(1,i).BPint,1)]);
    y1(1,i).BPint_resam=resample(y1(1,i).BPint,Maxlength,size(y1(1,i).BPint,1));
    y2(1,i).BPint_resam=resample(y2(1,i).BPint,Maxlength,size(y2(1,i).BPint,1));
    % smooth (length=Maxlength-BPwin)
    y1(1,i).BPwin=BPwin; y2(1,i).BPwin=BPwin;
    y1(1,i).BPint_resam_smoo=smoothing(((y1(1,i).BPint_resam)./2),BPwin); 
    y2(1,i).BPint_resam_smoo=smoothing(((y2(1,i).BPint_resam)./2),BPwin);

%     figure;
%     plot(y1(1,i).resam_BP_inter1,'r'); hold on; plot(y2(1,i).resam_BP_inter2,'b'); hold off;
%     ylabel('interval (ms)','color','m'); xlabel('taps (resampled)');
%     title({['Condi ' conditionNames{condiSeq(i)}],['BPwin ' num2str(BPwin)]},'Color',condicolors(condiSeq(i),:));
%     legend({'BP-L','BP-R'});

    % calculate correlation coeffiences win the same smoothing window
    BPint_corrSeries=[];
    for k=1:(Maxlength-BPwin-BPwin+1)
        BPint_corrSeries(k)=corr(y1(1,i).BPint_resam_smoo(k:k+BPwin-1),y2(1,i).BPint_resam_smoo(k:k+BPwin-1));
    end
    y12(1,i).BPint_corrSeries=BPint_corrSeries; % (length=Maxlength-BPwin-BPwin+1)

%     figure; 
%     plot(corrSeries,'k'); 
%     ylabel('corr coef'); xlabel('taps (resampled)');
%     corrSeries = corrSeries(~isnan(corrSeries));
%     median(corrSeries)
%     title(['corr coef smooth win ' num2str(BPwin)],'Color',condicolors(condiSeq(i),:));
%     ylim([-1 1]);
%     yline(0,'color',deepyellow);


    % extract EEG in each condition and normalized with pacing
    y1(1,i).EEG=zscore(EEG(1).EEG{i}./repmat(mean(EEG_pacing(1).EEG_pacing{i},1),size(EEG(1).EEG{i},1),1),[],1);
    y2(1,i).EEG=zscore(EEG(2).EEG{i}./repmat(mean(EEG_pacing(2).EEG_pacing{i},1),size(EEG(2).EEG{i},1),1),[],1);

    % extract EEG from 500 ms before each tap
    y1(1,i).EEGwin=EEGwin; y2(1,i).EEGwin=EEGwin;
    for j=2:y1(1,i).numtaps % start from the 2nd tap (length=numtap-1) (first element is empty)
        BPind=y1(1,i).BPind(j);
        y1(1,i).EEG500ms(1,j).EEG500ms=y1(1,i).EEG(BPind-999:BPind,:);
        delta=[];theta=[];alpha=[];beta=[];gamma=[];
        [delta, theta, alpha, beta, gamma]=sum5band(y1(1,i).EEG500ms(1,j).EEG500ms,sr,EEGwin);
        y1(1,i).EEG500ms(1,j).delta=delta;
        y1(1,i).EEG500ms(1,j).theta=theta;
        y1(1,i).EEG500ms(1,j).alpha=alpha;
        y1(1,i).EEG500ms(1,j).beta=beta;
        y1(1,i).EEG500ms(1,j).gamma=gamma;
    end

    for j=2:y2(1,i).numtaps
        BPind=y2(1,i).BPind(j);
        y2(1,i).EEG500ms(1,j).EEG500ms=y2(1,i).EEG(BPind-999:BPind,:);
        delta=[];theta=[];alpha=[];beta=[];gamma=[];
        [delta, theta, alpha, beta, gamma]=sum5band(y2(1,i).EEG500ms(1,j).EEG500ms,sr,EEGwin);
        y2(1,i).EEG500ms(1,j).delta=delta;
        y2(1,i).EEG500ms(1,j).theta=theta;
        y2(1,i).EEG500ms(1,j).alpha=alpha;
        y2(1,i).EEG500ms(1,j).beta=beta;
        y2(1,i).EEG500ms(1,j).gamma=gamma;
    end
end

% Combine power from each tap into one matrix
eeg1_mat=[];eeg2_mat=[];
for i=1:nTrials
    for j=2:y1(1,i).numtaps
        eeg1_mat(1,i).EEG500ms_mat(j,:,1)=y1(1,i).EEG500ms(1,j).delta;
        eeg1_mat(1,i).EEG500ms_mat(j,:,2)=y1(1,i).EEG500ms(1,j).theta;
        eeg1_mat(1,i).EEG500ms_mat(j,:,3)=y1(1,i).EEG500ms(1,j).alpha;
        eeg1_mat(1,i).EEG500ms_mat(j,:,4)=y1(1,i).EEG500ms(1,j).beta;
        eeg1_mat(1,i).EEG500ms_mat(j,:,5)=y1(1,i).EEG500ms(1,j).gamma;
    end
    for j=2:y2(1,i).numtaps
        eeg2_mat(1,i).EEG500ms_mat(j,:,1)=y2(1,i).EEG500ms(1,j).delta;
        eeg2_mat(1,i).EEG500ms_mat(j,:,2)=y2(1,i).EEG500ms(1,j).theta;
        eeg2_mat(1,i).EEG500ms_mat(j,:,3)=y2(1,i).EEG500ms(1,j).alpha;
        eeg2_mat(1,i).EEG500ms_mat(j,:,4)=y2(1,i).EEG500ms(1,j).beta;
        eeg2_mat(1,i).EEG500ms_mat(j,:,5)=y2(1,i).EEG500ms(1,j).gamma;
    end
end
% Downsample (not good) power to the same length as BPint_corrSerie (length=Maxlength-BPwin-BPwin+1) 
% Resample the MaxTaps according to results from previous step
for i=1:nTrials
    % Resam_length=length(BPint_corrSeries); % this will downsample, not good
    MaxTaps=max([y1(1,i).numtaps-1 y2(1,i).numtaps-1]); % start from 2nd tap
    eeg1_mat(1,i).EEG500ms_resam=resample(eeg1_mat(1,i).EEG500ms_mat(2:end,:,:),MaxTaps,y1(1,i).numtaps-1);
    eeg2_mat(1,i).EEG500ms_resam=resample(eeg2_mat(1,i).EEG500ms_mat(2:end,:,:),MaxTaps,y2(1,i).numtaps-1);
    eeg1_mat(1,i).EEG500ms_resam_smoo=smoothing(eeg1_mat(1,i).EEG500ms_resam,BPwin);
    eeg2_mat(1,i).EEG500ms_resam_smoo=smoothing(eeg2_mat(1,i).EEG500ms_resam,BPwin);
end

%% PLOT 4: 500ms before BP (examine one eeg band at a time)
mkdir('Plot4');cd Plot4

% Plot behaviral correlation
canvas(0.5,0.3);
for i=1:nTrials
    subplot(3,12,i); % BPint
    plot(y1(1,i).BPint_resam_smoo,'r.'); hold on; plot(y2(1,i).BPint_resam_smoo,'b.'); hold off;
    ylabel('interval (ms)'); xlabel('taps (resampled)');
    ylim([400 900]);xlim([0 size(y1(1,i).BPint_resam_smoo,1)]);
    title({[conditionNames{condiSeq(i)}],['BPwin ' num2str(BPwin)]},'Color',condicolors(condiSeq(i),:));
    % legend({'L','R'});
    
    subplot(3,12,12+i); % corrSeries
    plot(y12(1,i).BPint_corrSeries,'k'); 
    ylabel('corr coef'); xlabel(['taps (sum ' num2str(round(sum(y12(1,i).BPint_corrSeries,2))) ')'])
    % corrSeries = corrSeries(~isnan(corrSeries)); median(corrSeries)
    title({['corr coef'], ['smoo win ' num2str(BPwin)]},'Color',condicolors(condiSeq(i),:));
    ylim([-1 1]); yline(0,'color',deepyellow);
    xlim([0 size(y12(1,i).BPint_corrSeries,2)]);

    subplot(3,12,24+i); % BPint difference
    plot(abs(y1(1,i).BPint_resam_smoo-y2(1,i).BPint_resam_smoo),'k'); 
    ylabel('difference (ms)'); xlabel(['taps (sum ' num2str(round(sum(abs(y1(1,i).BPint_resam_smoo-y2(1,i).BPint_resam_smoo)))) ')'])
    ylim([-10 200]); yline(0,'color',deepyellow);xlim([0 size(y1(1,i).BPint_resam_smoo,1)]);
    title({[conditionNames{condiSeq(i)}],['BPwin ' num2str(BPwin)]},'Color',condicolors(condiSeq(i),:));
end
sgtitle(['session ' num2str(seed) ' BPint corr']);
    
figureName=['Plot4: BPint coor'];
saveas(gcf,figureName,'fig');



% Plot each power band
eegbandNames={'delta','theta','alpha','beta','gamma'};
selected_chan=[1:32]; % selected_chan=15;
ylims=[500 700 1400 1800 2000];
% seleted_band=4;
canvas(0.5,0.2);
for seleted_band=1:5
    for i=1:nTrials
        subplot(2,12,i); % power of one eeg band in 32 channels evolved over time in Player L
        plotx(eeg1_mat(1,i).EEG500ms_resam_smoo(:,selected_chan,seleted_band));
        ylim([0 ylims(seleted_band)]);
        xlim([0 size(eeg1_mat(1,i).EEG500ms_resam_smoo,1)]);
        title('L power','Color',condicolors(condiSeq(i),:));
        
        subplot(2,12,12+i); % power of one eeg band in 32 channels evolved over time in Player R
        plotx(eeg2_mat(1,i).EEG500ms_resam_smoo(:,selected_chan,seleted_band));
        ylim([0 ylims(seleted_band)]);
        xlim([0 size(eeg2_mat(1,i).EEG500ms_resam_smoo,1)]);
        title('R power','Color',condicolors(condiSeq(i),:));
    end
    sgtitle({['session ' num2str(seed) ' power band ' eegbandNames{seleted_band} ],[' chan ' labels{selected_chan}]});
    
    figureName=['Plot4: power band ' eegbandNames{seleted_band} ' chan ' labels{selected_chan}];
    saveas(gcf,figureName,'fig');
end

%% PLOT 5: GC between BPint and EEGpow
GC1=[];GC2=[];
eegbandNames={'delta','theta','alpha','beta','gamma'};
% selected_chan=15; % selected_chan=15;
% seleted_band=4;
tic 
for i=1:nTrials
    % L behaviral -> R EEG
    Ytmp1=[];Ytmp2=[];
    Ytmp1=y1(1,i).BPint; % length=numtaps-1
    Ytmp2=eeg2_mat(1,i).EEG500ms_mat(2:end,:,:); % start from 2nd tap
    MaxTaps=max([y1(1,i).numtaps-1 y2(1,i).numtaps-1]);
    % resample then smooth
    Ytmp1_resam=resample(Ytmp1,MaxTaps,length(Ytmp1));
    Ytmp2_resam=resample(Ytmp2,MaxTaps,length(Ytmp2));
    GC1(1,i).Ytmp1_resam_smoo=smoothing(Ytmp1_resam,BPwin);
    GC1(1,i).Ytmp2_resam_smoo=smoothing(Ytmp2_resam,BPwin);
    for seleted_band=1:5
        for selected_chan=1:32
            % d estimate and remove d from them
            Ytmp1=[];Ytmp2=[];
            [Ytmp1,Ytmp2] = rm_est_d(GC1(1,i).Ytmp1_resam_smoo,GC1(1,i).Ytmp2_resam_smoo(:,selected_chan,seleted_band));
            % Granger causality between BPint and EEGpow
            F=[];pval=[];,sig=[];
            [F,pval,sig] = myautocov(Ytmp1,Ytmp2); % Autocov method
            GC1(1,i).F(1,seleted_band).F(1,selected_chan).F=F;
            GC1(1,i).pval(1,seleted_band).pval(1,selected_chan).pval=pval;
            GC1(1,i).sig(1,seleted_band).sig(1,selected_chan).sig=sig;
        end
    end

    
    % R behaviral -> L EEG
    Ytmp1=[];Ytmp2=[];
    Ytmp1=y2(1,i).BPint; % length=numtaps-1
    Ytmp2=eeg1_mat(1,i).EEG500ms_mat(2:end,:,:); % start from 2nd tap
    MaxTaps=max([y1(1,i).numtaps-1 y2(1,i).numtaps-1]);
    % resample then smooth
    Ytmp1_resam=resample(Ytmp1,MaxTaps,length(Ytmp1));
    Ytmp2_resam=resample(Ytmp2,MaxTaps,length(Ytmp2));
    GC2(1,i).Ytmp1_resam_smoo=smoothing(Ytmp1_resam,BPwin);
    GC2(1,i).Ytmp2_resam_smoo=smoothing(Ytmp2_resam,BPwin);
    for seleted_band=1:5
        for selected_chan=1:32
            % d estimate and remove d from them
            Ytmp1=[];Ytmp2=[];
            [Ytmp1,Ytmp2] = rm_est_d(GC2(1,i).Ytmp1_resam_smoo,GC2(1,i).Ytmp2_resam_smoo(:,selected_chan,seleted_band));
            % Granger causality between BPint and EEGpow
            F=[];pval=[];,sig=[];
            [F,pval,sig] = myautocov(Ytmp1,Ytmp2); % Autocov method
            GC2(1,i).F(1,seleted_band).F(1,selected_chan).F=F;
            GC2(1,i).pval(1,seleted_band).pval(1,selected_chan).pval=pval;
            GC2(1,i).sig(1,seleted_band).sig(1,selected_chan).sig=sig;
        end
    end

end
toc

% organize all F(2) into matrics of 12trialx5bandx32chan
GC1mat=zeros(12,5,32);
for i=1:nTrials
    for seleted_band=1:5
        for selected_chan=1:32
            GC1mat(i,seleted_band,selected_chan)=GC1(i).F(seleted_band).F(selected_chan).F(2);
            GC2mat(i,seleted_band,selected_chan)=GC2(i).F(seleted_band).F(selected_chan).F(2);
        end
    end
end

%% PLOT 5: GC between BPint and EEGpow
mkdir('Plot5');cd Plot5

cd (['/ssd/zhibin/1overf/' num2str(seed) '_2P/Segmented_data/Plots/Plot5']);

canvas(0.5,0.5);
for i=1:nTrials
    subplot(2,12,i); % power of one eeg band in 32 channels evolved over time in Player L
    b=bar([GC1(i).F(2) GC1(i).F(3)]);set(gca, 'XTickLabel', {'L BP -> R EEG' 'R EEG -> L BP'});
    b.FaceColor='flat';b.CData(1,:) = pink;b.CData(2,:) = yellow;
    xlim([0.25 2.75]); 
    ylim([0 0.07]);
    title(conditionNames{condiSeq(i)},'Color',condicolors(condiSeq(i),:));
    if length(find(GC1(i).sig==1))==0;
        tempF=GC1(i).F(find(GC1(i).sig==0));
        tempP=round(GC1(i).pval(find(GC1(i).sig==0)),3);
        text((2-1-0.5), tempF(1), [num2str(round(tempP(1),2))],'Color',[1 0 0]);
        text((3-1-0.5), tempF(2), [num2str(round(tempP(2),2))],'Color',[0 0 1]);
    end
    if length(find(GC1(i).sig==1))==1;
        text((find(GC1(i).sig==1)-1-0.5), GC1(i).F(find(GC1(i).sig==1)), ['*' num2str(round(GC1(i).pval(find(GC1(i).sig==1)),2))],'Color',[1 0 0]);
        text((find(GC1(i).sig==0)-1-0.5), GC1(i).F(find(GC1(i).sig==0)), [num2str(round(GC1(i).pval(find(GC1(i).sig==0)),2))],'Color',[0 0 1]);
    end
    if length(find(GC1(i).sig==1))==2;
        tempF=GC1(i).F(find(GC1(i).sig==1));
        tempP=round(GC1(i).pval(find(GC1(i).sig==1)),3);
        text((2-1-0.5), tempF(1), ['*' num2str(round(tempP(1),2))],'Color',[1 0 0]);
        text((3-1-0.5), tempF(2), ['*' num2str(round(tempP(2),2))],'Color',[0 0 1]);
    end

    subplot(2,12,12+i); % power of one eeg band in 32 channels evolved over time in Player R
    b=bar([GC2(i).F(2) GC2(i).F(3)]);set(gca, 'XTickLabel', {'R BP -> L EEG' 'L EEG -> R BP'});
    b.FaceColor='flat';b.CData(1,:) = pink;b.CData(2,:) = yellow;
    xlim([0.25 2.75]); 
    ylim([0 0.07]);
    title(conditionNames{condiSeq(i)},'Color',condicolors(condiSeq(i),:));
    if length(find(GC2(i).sig==1))==0;
        tempF=GC2(i).F(find(GC2(i).sig==0));
        tempP=round(GC2(i).pval(find(GC2(i).sig==0)),3);
        text((2-1-0.5), tempF(1), [num2str(round(tempP(1),2))],'Color',[1 0 0]);
        text((3-1-0.5), tempF(2), [num2str(round(tempP(2),2))],'Color',[0 0 1]);
    end
    if length(find(GC2(i).sig==1))==1;
        text((find(GC2(i).sig==1)-1-0.5), GC2(i).F(find(GC2(i).sig==1)), ['*' num2str(round(GC2(i).pval(find(GC2(i).sig==1)),2))],'Color',[1 0 0]);
        text((find(GC2(i).sig==0)-1-0.5), GC2(i).F(find(GC2(i).sig==0)), [num2str(round(GC2(i).pval(find(GC2(i).sig==0)),2))],'Color',[0 0 1]);
    end
    if length(find(GC2(i).sig==1))==2;
        tempF=GC2(i).F(find(GC2(i).sig==1));
        tempP=round(GC2(i).pval(find(GC2(i).sig==1)),3);
        text((2-1-0.5), tempF(1), ['*' num2str(round(tempP(1),2))],'Color',[1 0 0]);
        text((3-1-0.5), tempF(2), ['*' num2str(round(tempP(2),2))],'Color',[0 0 1]);
    end
end
sgtitle({['session ' num2str(seed) ' power band ' eegbandNames{seleted_band} ],[' chan ' labels{selected_chan}]});

figureName=['Plot5: GC ' eegbandNames{seleted_band} ' chan ' labels{selected_chan}];
saveas(gcf,figureName,'fig');

% bar plot of 12 trials
% select a band and a chan for bar plot
seleted_band=3;selected_chan=15;
GC1vec=GC1mat(:,seleted_band,selected_chan);
GC2vec=GC1mat(:,seleted_band,selected_chan);
GC12=[GC1vec GC2vec];
figure;
b=bar(GC12);
b(1).FaceColor=red;b(2).FaceColor=blue;
ylabel('GC index');
xticks([1:12]);
xticklabels(conditionNames(condiSeq));
title('red: L BPint -> R EEG;      blue: R BPint -> L EEG');

% Topoplot 5 bands x 12trials (2 figures, one for each direction)
cd /ssd/zhibin/1overf/20220713_2P/Segmented_data/Plots/Plot5
for d=1:2
    if d==1
        GCmat=GC1mat;
    else
        GCmat=GC2mat;
    end
    canvas(0.8,0.5);
    for seleted_band=1:5
        for i=1:nTrials
            subplot(5,12,12*(seleted_band-1)+i)
            topoplot(squeeze(GCmat(i,seleted_band,:)),chaninfo,'nosedir','+X');
            title(['\' eegbandNames{seleted_band} ' ' conditionNames{condiSeq(i)}],'Color',condicolors(condiSeq(i),:));
            colorbar; caxis([-0.1 0.1]);
        end
    end
    
    if d==1
        sgtitle({['session ' num2str(seed) ' ' sessionTypes{session}  ' L BPint -> R EEG']},'color','red');
        figureName=['Plot5: GC  L BPint -> R EEG'];
    else
        sgtitle({['session ' num2str(seed) ' ' sessionTypes{session}  ' R BPint -> L EEG']},'color','blue');
        figureName=['Plot5: GC R BPint -> L EEG'];
    end
    saveas(gcf,figureName,'fig')
end


%% PLOT 6 DFA H values for all sessions based on matched intervals (sorted order)
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816];
sessions={'synch','synco','synch','synco','synch','synco','synch','synco'};
H=zeros(2,12,8);
cd /ssd/zhibin/1overf/
for r=1:8;
    clear intervals conditions sortorders
    runid = num2str(seeds(r,:));
    path = [runid '_2P/Cleaned_data/'];
    load([path  'clean_' runid '.mat']);
    % sort order and plot results
    [x,sortorder]=sort(conditions);
    for j=1:12
        clear intL_good_dmean intR_good_dmean
        % After aligment
        intervals{sortorder(j)}(:,1);
        intervals{sortorder(j)}(:,2);
        % remove the mean
        intL_good_dmean=intervals{sortorder(j)}(:,1)-mean(intervals{sortorder(j)}(:,1));
        intR_good_dmean=intervals{sortorder(j)}(:,2)-mean(intervals{sortorder(j)}(:,2));
        [~,H(1,j,r)]=DFA_main(intL_good_dmean);
        [~,H(2,j,r)]=DFA_main(intR_good_dmean);
    end
end

% print H to a table
Trials={'uncouple','uncouple','uncouple','L-lead','L-lead','L-lead',...
    'R-lead','R-lead','R-lead','mutual','mutual','mutual'}';
H_L_synch1=H(1,:,1)';H_R_synch1=H(2,:,1)';
H_L_synco1=H(1,:,2)';H_R_synco1=H(2,:,2)';
H_L_synch2=H(1,:,3)';H_R_synch2=H(2,:,3)';
H_L_synco2=H(1,:,4)';H_R_synco2=H(2,:,4)';
H_L_synch3=H(1,:,5)';H_R_synch3=H(2,:,5)';
H_L_synco3=H(1,:,6)';H_R_synco3=H(2,:,6)';
H_L_synch4=H(1,:,7)';H_R_synch4=H(2,:,7)';
H_L_synco4=H(1,:,8)';H_R_synco4=H(2,:,8)';
Table_H=table(Trials, H_L_synch1, H_R_synch1, H_L_synco1, H_R_synco1,...
    H_L_synch2, H_R_synch2, H_L_synco2, H_R_synco2,...
    H_L_synch3, H_R_synch3, H_L_synco3, H_R_synco3,...
    H_L_synch4, H_R_synch4, H_L_synco4, H_R_synco4);

% plot
canvas(0.5, 0.5);
for r=1:8
    subplot(8,1,r);
    y=[];
    y=[H(1,:,r)' H(2,:,r)'];
    b=bar(y);
    b(1).FaceColor = [1 0 0]; b(2).FaceColor = [0 0 1];
    ylabel('H');ylim([0.2 1.5]);
    xticks(1:12);xticklabels(Trials);
    title([num2str(seeds(r,:)) '-' sessions{r}]);
    hold on;
    yline(1,'color',deepyellow);
    yline(0.5,'color',deepyellow);
    legend('L','R','1','0.5','location','eastoutside');
end
sgtitle('H from DFA in all sessions')

% average all trials of the same condition
H_L_uncouple_synch = reshape(H(1,[1:3],[1:2:7]),1,[])  % L,uncouple,synch
H_L_Llead_synch    = reshape(H(1,[4:6],[1:2:7]),1,[])  % L,L lead,synch
H_L_Rlead_synch    = reshape(H(1,[7:9],[1:2:7]),1,[])  % L,R lead,synch
H_L_mutual_synch   = reshape(H(1,[10:12],[1:2:7]),1,[])% L, mutual,synch
H_R_uncouple_synch = reshape(H(2,[1:3],[1:2:7]),1,[])  % R,uncouple,synch
H_R_Llead_synch    = reshape(H(2,[4:6],[1:2:7]),1,[])  % R,L lead,synch
H_R_Rlead_synch    = reshape(H(2,[7:9],[1:2:7]),1,[])  % R,R lead,synch
H_R_mutual_synch   = reshape(H(2,[10:12],[1:2:7]),1,[])% R, mutual,synch
H_LR_synch_mean = [mean(H_L_uncouple_synch) mean(H_R_uncouple_synch);...
    mean(H_L_Llead_synch) mean(H_R_Llead_synch);...
    mean(H_L_Rlead_synch) mean(H_R_Rlead_synch);...
    mean(H_L_mutual_synch) mean(H_R_mutual_synch)];


H_L_uncouple_synco = reshape(H(1,[1:3],[2:2:8]),1,[])  % R,uncouple,synco
H_L_Llead_synco    = reshape(H(1,[4:6],[2:2:8]),1,[])  % R,L lead,synco
H_L_Rlead_synco    = reshape(H(1,[7:9],[2:2:8]),1,[])  % R,R lead,synco
H_L_mutual_synco   = reshape(H(1,[10:12],[2:2:8]),1,[])% R, mutual,synco
H_R_uncouple_synco = reshape(H(2,[1:3],[2:2:8]),1,[])  % R,uncouple,synco
H_R_Llead_synco    = reshape(H(2,[4:6],[2:2:8]),1,[])  % R,L lead,synco
H_R_Rlead_synco    = reshape(H(2,[7:9],[2:2:8]),1,[]) % R,R lead,synco
H_R_mutual_synco   = reshape(H(2,[10:12],[2:2:8]),1,[])% R, mutual,synco
H_LR_synco_mean = [mean(H_L_uncouple_synco) mean(H_R_uncouple_synco);...
    mean(H_L_Llead_synco) mean(H_R_Llead_synco);...
    mean(H_L_Rlead_synco) mean(H_R_Rlead_synco);...
    mean(H_L_mutual_synco) mean(H_R_mutual_synco)];

figure;
subplot(2,1,1);
b=bar(H_LR_synch_mean);b(1).FaceColor = [1 0 0]; b(2).FaceColor = [0 0 1];
xticks(1:4);xticklabels({'uncouple','L-lead','R-lead','mutual'});
ylabel('H');
hold on;
yline(1,'color',deepyellow);
yline(0.5,'color',deepyellow);
title('mean H in all synch sessions')
subplot(2,1,2);
b=bar(H_LR_synco_mean);b(1).FaceColor = [1 0 0]; b(2).FaceColor = [0 0 1];
xticks(1:4);xticklabels({'uncouple','L-lead','R-lead','mutual'});
ylabel('H');
hold on;
yline(1,'color',deepyellow);
yline(0.5,'color',deepyellow);
title('mean H in all synco sessions')

% combine synch and synco
H_LR_syn_mean=[mean([H_L_uncouple_synch H_L_uncouple_synco H_R_uncouple_synch H_R_uncouple_synco]);...
    mean([H_L_Llead_synch H_L_Llead_synco H_R_Rlead_synch H_R_Rlead_synco]);... 
    mean([H_R_Llead_synch H_R_Llead_synco H_L_Rlead_synch H_L_Rlead_synco]);...
    mean([H_L_mutual_synch H_L_mutual_synco H_R_mutual_synch H_R_mutual_synco])];
figure;
b=bar(H_LR_syn_mean);
xticks(1:4);xticklabels({'uncouple','Leading','Following','mutual'});
ylabel('H');
hold on;
ylim([0 0.9])
yline(0.5,'color',deepyellow);
title('mean H in all sessions');

% violin plots
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/violin
canvas(0.4,0.3);
subplot(2,1,1);
Y=[H_L_uncouple_synch', H_R_uncouple_synch', H_L_Llead_synch',H_R_Llead_synch',H_L_Rlead_synch',H_R_Rlead_synch',H_L_mutual_synch',H_R_mutual_synch'];
v=violin(Y,...
    'facecolor',[1 0 0;0 0 1;1 0 0;0 0 1;1 0 0;0 0 1;1 0 0;0 0 1],...
    'edgecolor','b',...
'bw',0.3,...
'mc','k',...
'medc','r--');
set(gca,'XtickLabel',{'uncouple','uncouple','L-lead','L-lead','R-lead','R-lead','mutual','mutual'});
ylabel('H');
title('H synch (L - red; R - blue)');
grid on;

subplot(2,1,2);
Y=[H_L_uncouple_synco', H_R_uncouple_synco', H_L_Llead_synco',H_R_Llead_synco',H_L_Rlead_synco',H_R_Rlead_synco',H_L_mutual_synco',H_R_mutual_synco'];
v=violin(Y,...
    'facecolor',[1 0 0;0 0 1;1 0 0;0 0 1;1 0 0;0 0 1;1 0 0;0 0 1],...
    'edgecolor','b',...
'bw',0.3,...
'mc','k',...
'medc','r--');
set(gca,'XtickLabel',{'uncouple','uncouple','L-lead','L-lead','R-lead','R-lead','mutual','mutual'});
ylabel('H');
title('H synco (L - red; R - blue)');
grid on;

% combine synch and synco
figure;
Y={[H_L_uncouple_synch H_L_uncouple_synco H_R_uncouple_synch H_R_uncouple_synco]',...
    [H_L_Llead_synch H_L_Llead_synco H_R_Rlead_synch H_R_Rlead_synco]',... 
    [H_R_Llead_synch H_R_Llead_synco H_L_Rlead_synch H_L_Rlead_synco]',...
    [H_L_mutual_synch H_L_mutual_synco H_R_mutual_synch H_R_mutual_synco]};
v=violin(Y,...
    'facecolor',[1 0 0;0 0 1;1 0 0;0 0 1;1 0 0;0 0 1;1 0 0;0 0 1],...
    'edgecolor','b',...
'bw',0.3,...
'mc','k',...
'medc','r--');
set(gca,'XtickLabel',{'uncouple','Leading','Following','mutual'});
ylabel('H');
title('mean H in all sessions');
grid on;


%% PLOT 6-1 DFA H values for all sessions based on original intervals (sorted order)
clear
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816];
sessions={'synch','synco','synch','synco','synch','synco','synch','synco'};

H=zeros(2,12,8);
cd /ssd/zhibin/1overf/
for r=1:8;
    clear intervals conditions sortorders
    runid = num2str(seeds(r,:));
    path = [runid '_2P/Segmented_data/'];
    load([path  'data_variables' runid '.mat']);
    % sort order and plot results
    [x,sortorder]=sort(condiSeq);
    for j = 1:12
        clear BPL BPR sampL sampR intL intR intL_dmean intR_dmean
	    BPL = [BP(1).BP{sortorder(j)}];
	    BPR = [BP(2).BP{sortorder(j)}];
        sampL = find(BPL==1);
        sampR = find(BPR==1);
        intL=diff(sampL);
        intR=diff(sampR);
        % remove the mean
        intL_dmean=intL-mean(intL);
        intR_dmean=intR-mean(intR);
        [~,H(1,j,r)]=DFA_main(intL_dmean);
        [~,H(2,j,r)]=DFA_main(intR_dmean);
    end
end

% run plots in previous section

%% PLOT 7 Autocorr summary in all sessions (sorted order)
clear
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816];
sessions={'synch','synco','synch','synco','synch','synco','synch','synco'};

Autocorrs=zeros(2,12,8);
cd /ssd/zhibin/1overf/
for r=1:8;
    clear intervals conditions sortorders
    runid = num2str(seeds(r,:));
    path = [runid '_2P/Segmented_data/'];
    load([path  'data_variables' runid '.mat']);
    % sort order and plot results
    [x,sortorder]=sort(condiSeq);
    for j = 1:12
        clear BPL BPR sampL sampR intL intR intL_dmean intR_dmean
	    BPL = [BP(1).BP{sortorder(j)}];
	    BPR = [BP(2).BP{sortorder(j)}];
        sampL = find(BPL==1);
        sampR = find(BPR==1);
        intL=diff(sampL);
        intR=diff(sampR);
        % remove the mean
        intL_dmean=intL-mean(intL);
        intR_dmean=intR-mean(intR);
        acf=[]; [acf] = autocorr(intL_dmean,20); Autocorrs(1,j,r)=sum(acf(2:end));
        acf=[]; [acf] = autocorr(intR_dmean,20); Autocorrs(2,j,r)=sum(acf(2:end));
    end
end

% print Autocorrs to a table
Trials={'uncouple','uncouple','uncouple','L-lead','L-lead','L-lead',...
    'R-lead','R-lead','R-lead','mutual','mutual','mutual'}';
Autocorrs_L_synch1=Autocorrs(1,:,1)';Autocorrs_R_synch1=Autocorrs(2,:,1)';
Autocorrs_L_synco1=Autocorrs(1,:,2)';Autocorrs_R_synco1=Autocorrs(2,:,2)';
Autocorrs_L_synch2=Autocorrs(1,:,3)';Autocorrs_R_synch2=Autocorrs(2,:,3)';
Autocorrs_L_synco2=Autocorrs(1,:,4)';Autocorrs_R_synco2=Autocorrs(2,:,4)';
Autocorrs_L_synch3=Autocorrs(1,:,5)';Autocorrs_R_synch3=Autocorrs(2,:,5)';
Autocorrs_L_synco3=Autocorrs(1,:,6)';Autocorrs_R_synco3=Autocorrs(2,:,6)';
Autocorrs_L_synch4=Autocorrs(1,:,7)';Autocorrs_R_synch4=Autocorrs(2,:,7)';
Autocorrs_L_synco4=Autocorrs(1,:,8)';Autocorrs_R_synco4=Autocorrs(2,:,8)';
Table_Autocorrs=table(Trials, Autocorrs_L_synch1, Autocorrs_R_synch1, Autocorrs_L_synco1, Autocorrs_R_synco1,...
    Autocorrs_L_synch2, Autocorrs_R_synch2, Autocorrs_L_synco2, Autocorrs_R_synco2,...
    Autocorrs_L_synch3, Autocorrs_R_synch3, Autocorrs_L_synco3, Autocorrs_R_synco3,...
    Autocorrs_L_synch4, Autocorrs_R_synch4, Autocorrs_L_synco4, Autocorrs_R_synco4);

% run similar plots in previous section (PLOT 6)
% plot
canvas(0.5, 0.5);
for r=1:8
    subplot(8,1,r);
    y=[];
    y=[Autocorrs(1,:,r)' Autocorrs(2,:,r)'];
    b=bar(y);
    b(1).FaceColor = [1 0 0]; b(2).FaceColor = [0 0 1];
    ylabel('Autocorr'); ylim([0 15]);
    xticks(1:12);xticklabels(Trials);
    title([num2str(seeds(r,:)) '-' sessions{r}]);
%     hold on;
%     yline(1,'color',deepyellow);
%     yline(0.5,'color',deepyellow);
%     legend('L','R','1','0.5','location','eastoutside');
end
sgtitle('Sum of Autocorr (20 lags) in all sessions')

% average all trials of the same condition
Autocorrs_L_uncouple_synch = reshape(Autocorrs(1,[1:3],[1:2:7]),1,[])  % L,uncouple,synch
Autocorrs_L_Llead_synch    = reshape(Autocorrs(1,[4:6],[1:2:7]),1,[])  % L,L lead,synch
Autocorrs_L_Rlead_synch    = reshape(Autocorrs(1,[7:9],[1:2:7]),1,[])  % L,R lead,synch
Autocorrs_L_mutual_synch   = reshape(Autocorrs(1,[10:12],[1:2:7]),1,[])% L, mutual,synch
Autocorrs_R_uncouple_synch = reshape(Autocorrs(2,[1:3],[1:2:7]),1,[])  % R,uncouple,synch
Autocorrs_R_Llead_synch    = reshape(Autocorrs(2,[4:6],[1:2:7]),1,[])  % R,L lead,synch
Autocorrs_R_Rlead_synch    = reshape(Autocorrs(2,[7:9],[1:2:7]),1,[])  % R,R lead,synch
Autocorrs_R_mutual_synch   = reshape(Autocorrs(2,[10:12],[1:2:7]),1,[])% R, mutual,synch
Autocorrs_LR_synch_mean = [mean(Autocorrs_L_uncouple_synch) mean(Autocorrs_R_uncouple_synch);...
    mean(Autocorrs_L_Llead_synch) mean(Autocorrs_R_Llead_synch);...
    mean(Autocorrs_L_Rlead_synch) mean(Autocorrs_R_Rlead_synch);...
    mean(Autocorrs_L_mutual_synch) mean(Autocorrs_R_mutual_synch)];


Autocorrs_L_uncouple_synco = reshape(Autocorrs(1,[1:3],[2:2:8]),1,[])  % R,uncouple,synco
Autocorrs_L_Llead_synco    = reshape(Autocorrs(1,[4:6],[2:2:8]),1,[])  % R,L lead,synco
Autocorrs_L_Rlead_synco    = reshape(Autocorrs(1,[7:9],[2:2:8]),1,[])  % R,R lead,synco
Autocorrs_L_mutual_synco   = reshape(Autocorrs(1,[10:12],[2:2:8]),1,[])% R, mutual,synco
Autocorrs_R_uncouple_synco = reshape(Autocorrs(2,[1:3],[2:2:8]),1,[])  % R,uncouple,synco
Autocorrs_R_Llead_synco    = reshape(Autocorrs(2,[4:6],[2:2:8]),1,[])  % R,L lead,synco
Autocorrs_R_Rlead_synco    = reshape(Autocorrs(2,[7:9],[2:2:8]),1,[]) % R,R lead,synco
Autocorrs_R_mutual_synco   = reshape(Autocorrs(2,[10:12],[2:2:8]),1,[])% R, mutual,synco
Autocorrs_LR_synco_mean = [mean(Autocorrs_L_uncouple_synco) mean(Autocorrs_R_uncouple_synco);...
    mean(Autocorrs_L_Llead_synco) mean(Autocorrs_R_Llead_synco);...
    mean(Autocorrs_L_Rlead_synco) mean(Autocorrs_R_Rlead_synco);...
    mean(Autocorrs_L_mutual_synco) mean(Autocorrs_R_mutual_synco)];

figure;
subplot(2,1,1);
b=bar(Autocorrs_LR_synch_mean);b(1).FaceColor = [1 0 0]; b(2).FaceColor = [0 0 1];
xticks(1:4);xticklabels({'uncouple','L-lead','R-lead','mutual'});
ylabel('\rho(k)');
title('mean autocorr sum for 20 lags in all synch sessions')
subplot(2,1,2);
b=bar(Autocorrs_LR_synco_mean);b(1).FaceColor = [1 0 0]; b(2).FaceColor = [0 0 1];
xticks(1:4);xticklabels({'uncouple','L-lead','R-lead','mutual'});
ylabel('\rho(k)');
title('mean autocorr sum for 20 lags in all synco sessions')

% combine synch and synco
Autocorrs_LR_syn_mean=[mean([Autocorrs_L_uncouple_synch Autocorrs_L_uncouple_synco Autocorrs_R_uncouple_synch Autocorrs_R_uncouple_synco]);...
    mean([Autocorrs_L_Llead_synch Autocorrs_L_Llead_synco Autocorrs_R_Rlead_synch Autocorrs_R_Rlead_synco]);... 
    mean([Autocorrs_R_Llead_synch Autocorrs_R_Llead_synco Autocorrs_L_Rlead_synch Autocorrs_L_Rlead_synco]);...
    mean([Autocorrs_L_mutual_synch Autocorrs_L_mutual_synco Autocorrs_R_mutual_synch Autocorrs_R_mutual_synco])];
figure;
b=bar(Autocorrs_LR_syn_mean);
xticks(1:4);xticklabels({'uncouple','Leading','Following','mutual'});
ylabel('\rho(k)');
hold on;
title('mean autocorr sum for 20 lags in all sessions');

% combine synch and synco
figure;
Y={[Autocorrs_L_uncouple_synch Autocorrs_L_uncouple_synco Autocorrs_R_uncouple_synch Autocorrs_R_uncouple_synco]',...
    [Autocorrs_L_Llead_synch Autocorrs_L_Llead_synco Autocorrs_R_Rlead_synch Autocorrs_R_Rlead_synco]',... 
    [Autocorrs_R_Llead_synch Autocorrs_R_Llead_synco Autocorrs_L_Rlead_synch Autocorrs_L_Rlead_synco]',...
    [Autocorrs_L_mutual_synch Autocorrs_L_mutual_synco Autocorrs_R_mutual_synch Autocorrs_R_mutual_synco]};
v=violin(Y,...
    'edgecolor','b',...
'bw',0.3,...
'mc','k',...
'medc','r--');
set(gca,'XtickLabel',{'uncouple','Leading','Following','mutual'});
ylabel('\rho(k)');
title('mean autocorr sum for 20 lags in all sessions');
grid on;

%% PLOT 8 DFA after d removed (sorted order)
clear
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816];
sessions={'synch','synco','synch','synco','synch','synco','synch','synco'};

H=zeros(2,12,8);
cd /ssd/zhibin/1overf/
for r=1:8;
    clear intervals conditions sortorders
    runid = num2str(seeds(r,:));
    path = [runid '_2P/Segmented_data/'];
    load([path  'data_variables' runid '.mat']);
    % sort order and plot results
    [x,sortorder]=sort(condiSeq);
    for j = 1:12
        clear BPL BPR sampL sampR intL intR intL_dmean intR_dmean
	    BPL = [BP(1).BP{sortorder(j)}];
	    BPR = [BP(2).BP{sortorder(j)}];
        sampL = find(BPL==1);
        sampR = find(BPR==1);
        intL=diff(sampL);
        intR=diff(sampR);
        % remove the mean
        intL_dmean=intL-mean(intL);
        intR_dmean=intR-mean(intR);
        % estimate H
        [~,H(1,j,r)]=DFA_main(intL_dmean);
        [~,H(2,j,r)]=DFA_main(intR_dmean);
    end
end
d1=H-0.5;% estimate d using DFA method (method1: based on DFA, d=h-0.5)
d_min=-0.5; d_max=1; % (method2: open d_estimation)

% d removal
H_est1=zeros(2,12,8);
H_est2=zeros(2,12,8);
int_dmean_drm1=cell(2,12,8);
cd /ssd/zhibin/1overf/
for r=1:8;
    clear intervals conditions sortorders
    runid = num2str(seeds(r,:));
    path = [runid '_2P/Segmented_data/'];
    load([path  'data_variables' runid '.mat']);
    % sort order and plot results
    [x,sortorder]=sort(condiSeq);
    for j = 1:12
        clear BPL BPR sampL sampR intL intR intL_dmean intR_dmean
	    BPL = [BP(1).BP{sortorder(j)}];
	    BPR = [BP(2).BP{sortorder(j)}];
        sampL = find(BPL==1);
        sampR = find(BPR==1);
        intL=diff(sampL);
        intR=diff(sampR);
        % remove the mean
        intL_dmean=intL-mean(intL);
        intR_dmean=intR-mean(intR);
        % d removal  (method1: based on DFA, d=h-0.5)
        [int_dmean_drm1{1,j,r}]=remove_d(intL_dmean,d1(1,j,r));
        [int_dmean_drm1{2,j,r}]=remove_d(intR_dmean,d1(2,j,r));
        % % d removal  (method2: open d_estimation)
        [int_dmean_drm2{1,j,r}]=remove_d(intL_dmean,d_estimation(intL_dmean,d_min,d_max));
        [int_dmean_drm2{2,j,r}]=remove_d(intR_dmean,d_estimation(intR_dmean,d_min,d_max));
        % H estimate again using DFA method
        [~,H_est1(1,j,r)]=DFA_main(int_dmean_drm1{1,j,r});
        [~,H_est1(2,j,r)]=DFA_main(int_dmean_drm1{2,j,r});
        [~,H_est2(1,j,r)]=DFA_main(int_dmean_drm2{1,j,r});
        [~,H_est2(2,j,r)]=DFA_main(int_dmean_drm2{2,j,r});
    end
end
% d removal (based on DFA, d=h-0.5)


% print H_est2 to a table
Trials={'uncouple','uncouple','uncouple','L-lead','L-lead','L-lead',...
    'R-lead','R-lead','R-lead','mutual','mutual','mutual'}';
H_est2_L_synch1=H_est2(1,:,1)';H_est2_R_synch1=H_est2(2,:,1)';
H_est2_L_synco1=H_est2(1,:,2)';H_est2_R_synco1=H_est2(2,:,2)';
H_est2_L_synch2=H_est2(1,:,3)';H_est2_R_synch2=H_est2(2,:,3)';
H_est2_L_synco2=H_est2(1,:,4)';H_est2_R_synco2=H_est2(2,:,4)';
H_est2_L_synch3=H_est2(1,:,5)';H_est2_R_synch3=H_est2(2,:,5)';
H_est2_L_synco3=H_est2(1,:,6)';H_est2_R_synco3=H_est2(2,:,6)';
H_est2_L_synch4=H_est2(1,:,7)';H_est2_R_synch4=H_est2(2,:,7)';
H_est2_L_synco4=H_est2(1,:,8)';H_est2_R_synco4=H_est2(2,:,8)';
Table_H_est2=table(Trials, H_est2_L_synch1, H_est2_R_synch1, H_est2_L_synco1, H_est2_R_synco1,...
    H_est2_L_synch2, H_est2_R_synch2, H_est2_L_synco2, H_est2_R_synco2,...
    H_est2_L_synch3, H_est2_R_synch3, H_est2_L_synco3, H_est2_R_synco3,...
    H_est2_L_synch4, H_est2_R_synch4, H_est2_L_synco4, H_est2_R_synco4);

% plot H_est2 after
canvas(0.5, 0.5);
for r=1:8
    subplot(8,1,r);
    y=[];
    y=[H_est2(1,:,r)' H_est2(2,:,r)'];
    b=bar(y);
    b(1).FaceColor = [1 0 0]; b(2).FaceColor = [0 0 1];
    ylabel('H');ylim([0 1.5]);
    xticks(1:12);xticklabels(Trials);
    title([num2str(seeds(r,:)) '-' sessions{r}]);
    hold on;
    yline(1,'color',deepyellow);
    yline(0.5,'color',deepyellow);
    legend('L','R','1','0.5','location','eastoutside');
end
sgtitle('H from DFA after d removal in all sessions (using Leites d estimation function)');


%% PLOT 9 Xcorr after d removal and save new intervals (sorted order) 
clear
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816];
sessions={'synch','synco','synch','synco','synch','synco','synch','synco'};
Trials={'uncouple','uncouple','uncouple','L-lead','L-lead','L-lead',...
    'R-lead','R-lead','R-lead','mutual','mutual','mutual'}';

% d estimation from H (DFA method) based on the good matched intervals
H=zeros(2,12,8);XcorrPeakLag=zeros(12,8);XcorrPeak=zeros(12,8);
cd /ssd/zhibin/1overf/
for r=1:8;
    clear intervals conditions sortorders
    runid = num2str(seeds(r,:));
    path = [runid '_2P/Cleaned_data/'];
    load([path  'clean_' runid '.mat']);
    % sort order and plot results
    [x,sortorder]=sort(conditions);
    for j = 1:12
        clear intL_good_dmean intR_good_dmean
        % After aligment
        % intervals{sortorder(j)}(:,1);
        % intervals{sortorder(j)}(:,2);
        % remove the mean
        intL_good_dmean=intervals{sortorder(j)}(:,1)-mean(intervals{sortorder(j)}(:,1));
        intR_good_dmean=intervals{sortorder(j)}(:,2)-mean(intervals{sortorder(j)}(:,2));
        % estimate the d
        [~,H(1,j,r)]=DFA_main(intL_good_dmean);
        [~,H(2,j,r)]=DFA_main(intR_good_dmean);
        % Xcorr based on int_dmean_drm
        r12=[];lags12=[];
        [r12,lags12]=xcorr(intL_good_dmean,intR_good_dmean,10,'normalized');
        XcorrPeakLag(j,r)=lags12(find(r12==max(r12)));
        XcorrPeak(j,r)=max(r12);
    end
end
eval(['save ' '/ssd/zhibin/1overf/all_session20220713_0816/H.mat H'])
d1=H-0.5;% estimate d using DFA method (method1: based on DFA, d=h-0.5)
        
        
% d removal
H_est2=zeros(2,12,8);
int_dmean_drm=cell(2,12,8);
cd /ssd/zhibin/1overf/
for r=1:8;
    clear intervals conditions sortorders
    runid = num2str(seeds(r,:));
    path = [runid '_2P/Cleaned_data/'];
    load([path  'clean_' runid '.mat']);
    % sort order and plot results
    [x,sortorder]=sort(conditions);
    for j = 1:12
        clear intL_good_dmean intR_good_dmean
        % After aligment
        % intervals{sortorder(j)}(:,1);
        % intervals{sortorder(j)}(:,2);
        % remove the mean
        intL_good_dmean=intervals{sortorder(j)}(:,1)-mean(intervals{sortorder(j)}(:,1));
        intR_good_dmean=intervals{sortorder(j)}(:,2)-mean(intervals{sortorder(j)}(:,2));
        % d removal  (method1: based on DFA, d=h-0.5)
        [int_dmean_drm{1,j,r}]=remove_d(intL_good_dmean,d1(1,j,r));
        [int_dmean_drm{2,j,r}]=remove_d(intR_good_dmean,d1(2,j,r));
        % H estimate again using DFA method
        [~,H_est2(1,j,r)]=DFA_main(int_dmean_drm{1,j,r});
        [~,H_est2(2,j,r)]=DFA_main(int_dmean_drm{2,j,r});
    end
end
% save the "int_dmean_drm" under [/ssd/zhibin/1overf/all_session20220713_0816]
eval(['save ' '/ssd/zhibin/1overf/all_session20220713_0816/int_dmean_drm.mat int_dmean_drm'])
% plot H all session and compare like in PLOT 8

% Xcorr based on int_dmean_drm
XcorrPeakLag=zeros(12,8);XcorrPeak=zeros(12,8);
for r=1:8;
    for j = 1:12
        r12=[];lags12=[];
        [r12,lags12]=xcorr(int_dmean_drm{1,j,r},int_dmean_drm{2,j,r},10,'normalized');
        XcorrPeakLag(j,r)=lags12(find(r12==max(r12)));
        XcorrPeak(j,r)=max(r12);
    end
end

% Plot XcorrPeakLag in all sessions
canvas(0.3, 0.4);
for r=1:8
    subplot(8,1,r);
    plot(1:12,XcorrPeakLag(:,r),'ko')
    ylabel('lag');ylim([-10 10]);
    yline(0,'color',deepyellow);set(gca, 'YDir','reverse')
    xticks(1:12);xticklabels(Trials);
    title([num2str(seeds(r,:)) '-' sessions{r}]);
end
sgtitle('peak lag from Xcorr in all sessions (-lag = L leading; +lag = R leading)')

% print XcorrPeakLag to a table
XcorrPeakLag_synch1=XcorrPeakLag(:,1);
XcorrPeakLag_synco1=XcorrPeakLag(:,2);
XcorrPeakLag_synch2=XcorrPeakLag(:,3);
XcorrPeakLag_synco2=XcorrPeakLag(:,4);
XcorrPeakLag_synch3=XcorrPeakLag(:,5);
XcorrPeakLag_synco3=XcorrPeakLag(:,6);
XcorrPeakLag_synch4=XcorrPeakLag(:,7);
XcorrPeakLag_synco4=XcorrPeakLag(:,8);
Table_XcorrPeakLag=table(Trials, XcorrPeakLag_synch1,XcorrPeakLag_synco1,...
    XcorrPeakLag_synch2, XcorrPeakLag_synco2,...
    XcorrPeakLag_synch3, XcorrPeakLag_synco3,...
    XcorrPeakLag_synch4, XcorrPeakLag_synco4);

% average all trials of the same condition (lag of peak)
XcorrPeakLag_uncouple_synch = reshape(XcorrPeakLag([1:3],[1:2:7]),1,[])  % uncouple,synch
XcorrPeakLag_Llead_synch    = reshape(XcorrPeakLag([4:6],[1:2:7]),1,[])  % L lead,synch
XcorrPeakLag_Rlead_synch    = reshape(XcorrPeakLag([7:9],[1:2:7]),1,[])  % R lead,synch
XcorrPeakLag_mutual_synch   = reshape(XcorrPeakLag([10:12],[1:2:7]),1,[])% mutual,synch
XcorrPeakLag_LR_synch_mean = [mean(XcorrPeakLag_uncouple_synch);...
    mean(XcorrPeakLag_Llead_synch);...
    mean(XcorrPeakLag_Rlead_synch);...
    mean(XcorrPeakLag_mutual_synch)];
XcorrPeakLag_LR_synch_std = [std(XcorrPeakLag_uncouple_synch);...
    std(XcorrPeakLag_Llead_synch);...
    std(XcorrPeakLag_Rlead_synch);...
    std(XcorrPeakLag_mutual_synch)];

XcorrPeakLag_uncouple_synco = reshape(XcorrPeakLag([1:3],[2:2:8]),1,[])  % R,uncouple,synco
XcorrPeakLag_Llead_synco    = reshape(XcorrPeakLag([4:6],[2:2:8]),1,[])  % R,L lead,synco
XcorrPeakLag_Rlead_synco    = reshape(XcorrPeakLag([7:9],[2:2:8]),1,[])  % R,R lead,synco
XcorrPeakLag_mutual_synco   = reshape(XcorrPeakLag([10:12],[2:2:8]),1,[])% R, mutual,synco
XcorrPeakLag_LR_synco_mean = [mean(XcorrPeakLag_uncouple_synco);...
    mean(XcorrPeakLag_Llead_synco);...
    mean(XcorrPeakLag_Rlead_synco);...
    mean(XcorrPeakLag_mutual_synco)];
XcorrPeakLag_LR_synco_std = [std(XcorrPeakLag_uncouple_synco);...
    std(XcorrPeakLag_Llead_synco);...
    std(XcorrPeakLag_Rlead_synco);...
    std(XcorrPeakLag_mutual_synco)];
% plot average lag of peak
canvas(0.2, 0.3);
subplot(2,1,1);
errorbar(1:4,XcorrPeakLag_LR_synch_mean,XcorrPeakLag_LR_synch_std,'ko');
xticks(1:4);xticklabels({'uncouple','L-lead','R-lead','mutual'});
ylabel('average lag of the peak'); % ylabel('\rho(k)');
ylim([-10 10]);yline(0,'color',deepyellow);set(gca, 'YDir','reverse');
xlim([0 5]);
title('average lag of xcorr peak in all synch sessions (-lag: L leading; +lag: R leading)')
subplot(2,1,2);
errorbar(1:4,XcorrPeakLag_LR_synco_mean,XcorrPeakLag_LR_synco_std,'ko');
xticks(1:4);xticklabels({'uncouple','L-lead','R-lead','mutual'});
ylabel('average lag of the peak'); % ylabel('\rho(k)');
ylim([-10 10]);yline(0,'color',deepyellow);set(gca, 'YDir','reverse')
xlim([0 5]);
title('average lag of xcorr peak in all synco sessions (-lag: L leading; +lag: R leading)')

% average all trials of the same condition (lag of peak)
XcorrPeak_uncouple_synch = reshape(XcorrPeak([1:3],[1:2:7]),1,[])  % uncouple,synch
XcorrPeak_Llead_synch    = reshape(XcorrPeak([4:6],[1:2:7]),1,[])  % L lead,synch
XcorrPeak_Rlead_synch    = reshape(XcorrPeak([7:9],[1:2:7]),1,[])  % R lead,synch
XcorrPeak_mutual_synch   = reshape(XcorrPeak([10:12],[1:2:7]),1,[])% mutual,synch
XcorrPeak_LR_synch_mean = [mean(XcorrPeak_uncouple_synch);...
    mean(XcorrPeak_Llead_synch);...
    mean(XcorrPeak_Rlead_synch);...
    mean(XcorrPeak_mutual_synch)];
XcorrPeak_LR_synch_std = [std(XcorrPeak_uncouple_synch);...
    std(XcorrPeak_Llead_synch);...
    std(XcorrPeak_Rlead_synch);...
    std(XcorrPeak_mutual_synch)];

XcorrPeak_uncouple_synco = reshape(XcorrPeak([1:3],[2:2:8]),1,[])  % R,uncouple,synco
XcorrPeak_Llead_synco    = reshape(XcorrPeak([4:6],[2:2:8]),1,[])  % R,L lead,synco
XcorrPeak_Rlead_synco    = reshape(XcorrPeak([7:9],[2:2:8]),1,[])  % R,R lead,synco
XcorrPeak_mutual_synco   = reshape(XcorrPeak([10:12],[2:2:8]),1,[])% R, mutual,synco
XcorrPeak_LR_synco_mean = [mean(XcorrPeak_uncouple_synco);...
    mean(XcorrPeak_Llead_synco);...
    mean(XcorrPeak_Rlead_synco);...
    mean(XcorrPeak_mutual_synco)];
XcorrPeak_LR_synco_std = [std(XcorrPeak_uncouple_synco);...
    std(XcorrPeak_Llead_synco);...
    std(XcorrPeak_Rlead_synco);...
    std(XcorrPeak_mutual_synco)];

% plot average peak
canvas(0.2, 0.3);
subplot(2,1,1);
errorbar(1:4,XcorrPeak_LR_synch_mean,XcorrPeak_LR_synch_std,'ko');
xticks(1:4);xticklabels({'uncouple','L-lead','R-lead','mutual'});
ylabel('\rho(k)');
ylim([-0.1 1]);yline(0,'color',red);xlim([0 5]);
title('average peak of xcorr in all synch sessions');
subplot(2,1,2);
errorbar(1:4,XcorrPeak_LR_synco_mean,XcorrPeak_LR_synco_std,'ko');
xticks(1:4);xticklabels({'uncouple','L-lead','R-lead','mutual'});
ylabel('\rho(k)');
ylim([-0.1 1]);yline(0,'color',red);xlim([0 5]);
title('average peak of xcorr in all synco sessions');


%% PLOT 10 redo H and intervals_H_removed separately (based on matched invervals)
clear
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816];
sessions={'synch','synco','synch','synco','synch','synco','synch','synco'};
cd /ssd/zhibin/1overf/
H=[];intervals_H_removed=[];
for s=1:8
    clear intervals H intervals_H_removed
    runid=num2str(seeds(s,:));
    load(['/ssd/zhibin/1overf/' runid '_2P/Cleaned_data/clean_' runid '.mat' ])   
    for b=1:12
%         % After aligment
%         intervals{b}(:,1);
%         intervals{b}(:,2);
%         % remove the mean (no need for DFA)
%         intL_good_dmean=intervals{b}(:,1)-mean(intervals{b}(:,1));
%         intR_good_dmean=intervals{b}(:,2)-mean(intervals{b}(:,2));
        [~,H(1,b)]=DFA_main(intervals{b}(:,1));
        [~,H(2,b)]=DFA_main(intervals{b}(:,2));
        % d removal
        [intervals_H_removed{b}(:,1)]=remove_d(intervals{b}(:,1),H(1,b)-0.5);
        [intervals_H_removed{b}(:,2)]=remove_d(intervals{b}(:,2),H(2,b)-0.5);
        % save H and intervals_H_removed
        eval(['save ' '/ssd/zhibin/1overf/' runid '_2P/Cleaned_data/d_removal' runid '.mat H intervals_H_removed'])
    end
end
%% PLOT 10-1 get H from all sessions (based on matched invervals)
clear
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816];
sessions={'synch','synco','synch','synco','synch','synco','synch','synco'};
cd /ssd/zhibin/1overf/
H_all=[];
for s=1:8
    runid=num2str(seeds(s,:));
    clear H
    load(['/ssd/zhibin/1overf/' runid '_2P/Cleaned_data/d_removal' runid '.mat' ])   
    for b=1:12
        H_all(:,s,b)=H(:,b);
    end
end

% H-int matching for L and R as in Marmelat 2012
H_all; % 2x8x12
H_Lall=squeeze(H_all(1,:,:));
H_Lall=reshape(H_Lall',[],1);
H_Rall=squeeze(H_all(2,:,:));
H_Rall=reshape(H_Rall',[],1);
Inds4(:,c);
uncoupleInd(1:24);
leadingInd(1:24);
followingInd(1:24);
mutualInd(1:24);

canvas(0.6,0.2);
% canvas(0.2,0.1);
% subplot(1,2,1);
subplot(1,3,1);
plot(H_Lall(uncoupleInd(1:24)),H_Rall(uncoupleInd(1:24)),'ko');
xlabel('DFA exponenet, Pacipant L');ylabel('DFA exponenet, Pacipant R');
title('uncouple');
A=[];Alpha1=[];FitValues=[];
A=polyfit(H_Lall(uncoupleInd(1:24)),H_Rall(uncoupleInd(1:24)),1);
Alpha1=A(1); % the slope, the first order polynomial coefficient from polyfit (Hurst Componenet >1 ?)
FitValues=polyval(A,H_Lall(uncoupleInd(1:24)));
hold on; plot(H_Lall(uncoupleInd(1:24)),FitValues,'k-');
xlim([0.4 1.4]);ylim([0.4 1.4]);plot([0 1.4], [0 1.4],'m--');hold off;
subplot(1,3,2);
plot([H_Lall(leadingInd(1:24)); H_Lall(followingInd(1:24))], [H_Rall(leadingInd(1:24)); H_Rall(followingInd(1:24))],'ko');
xlabel('DFA exponenet, Pacipant L');ylabel('DFA exponenet, Pacipant R');
title('uni-directional coupling (leading or following)');
A=[];Alpha1=[];FitValues=[];
A=polyfit([H_Lall(leadingInd(1:24)); H_Lall(followingInd(1:24))], [H_Rall(leadingInd(1:24)); H_Rall(followingInd(1:24))],1);
Alpha1=A(1); % the slope, the first order polynomial coefficient from polyfit (Hurst Componenet >1 ?)
FitValues=polyval(A,[H_Lall(leadingInd(1:24)); H_Lall(followingInd(1:24))]);
hold on; plot([H_Lall(leadingInd(1:24)); H_Lall(followingInd(1:24))],FitValues,'k-');
xlim([0.4 1.4]);ylim([0.4 1.4]);plot([0 1.4], [0 1.4],'m--');hold off;
% plot(H_Lall(followingInd(1:24)),H_Rall(followingInd(1:24)),'ko');
% xlabel('DFA exponenet, Pacipant L');ylabel('DFA exponenet, Pacipant R');
% title('following');
% A=[];Alpha1=[];FitValues=[];
% A=polyfit(H_Lall(followingInd(1:24)),H_Rall(followingInd(1:24)),1);
% Alpha1=A(1); % the slope, the first order polynomial coefficient from polyfit (Hurst Componenet >1 ?)
% FitValues=polyval(A,H_Lall(followingInd(1:24)));
% hold on; plot(H_Lall(followingInd(1:24)),FitValues,'k-');
% subplot(1,2,2);
subplot(1,3,3);
plot(H_Lall(mutualInd(1:24)),H_Rall(mutualInd(1:24)),'ko');
xlabel('DFA exponenet, Pacipant L');ylabel('DFA exponenet, Pacipant R');
title('mutual');
A=[];Alpha1=[];FitValues=[];
A=polyfit(H_Lall(mutualInd(1:24)),H_Rall(mutualInd(1:24)),1);
Alpha1=A(1); 
FitValues=polyval(A,H_Lall(mutualInd(1:24)));
hold on; plot(H_Lall(mutualInd(1:24)),FitValues,'k-');
xlim([0.4 1.4]);ylim([0.4 1.4]);plot([0 1.4], [0 1.4],'m--');hold off;
%% PLOT 11 Xcorr slinding window of 20 intervals (based on the matched invervals)
clear
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816];
sessions={'synch','synco','synch','synco','synch','synco','synch','synco'};
% sliding win
BPwin=20;
% pick one session to analyze
s=2;
runid=num2str(seeds(s,:));
load(['/ssd/zhibin/1overf/' runid '_2P/Cleaned_data/clean_' runid '.mat' ])
load(['/ssd/zhibin/1overf/' runid '_2P/Cleaned_data/d_removal' runid '.mat' ])

BPint_corrSeries={};intNum=[];BPint_xcorrSeries={}
for b=1:12
    % remove the mean (no need for DFA)
    intL_good_dmean=intervals{b}(:,1)-mean(intervals{b}(:,1));
    intR_good_dmean=intervals{b}(:,2)-mean(intervals{b}(:,2));
    % sliding window of 20 to compute xcorr at lag 0
    intNum(b)=size(intervals{b},1);
    for k=1:(intNum(b)-BPwin+1)
        % corr
        BPint_corrSeries{b}(k)=corr(intL_good_dmean(k:k+BPwin-1),intR_good_dmean(k:k+BPwin-1));
        % xcorr
        r12=[];lags12=[];
        [r12,lags12]=xcorr(intL_good_dmean(k:k+BPwin-1),intR_good_dmean(k:k+BPwin-1),10,'normalized');
        BPint_xcorrSeries{b}(k)=r12(11);
    end
end


% Plot behaviral correlation for one subject
canvas(1,0.5);
for b=1:12
    subplot(4,12,b); % BPint
    plot(intervals{b}(:,1)/2,'r'); hold on; plot(intervals{b}(:,2)/2,'b'); hold off;
    ylabel('interval (ms)'); xlabel('taps (matched)');
    ylim([400 1100]);xlim([0 intNum(b)]);
    title({[conditionNames{conditions(b)}],['BPwin ' num2str(BPwin)]},'Color',condicolors(conditions(b),:));
    % legend({'L','R'});
    
    subplot(4,12,12+b); % corrSeries
    plot(BPint_corrSeries{b},'k'); 
    ylabel('corr coef'); xlabel(['taps (sum ' num2str(round(sum(BPint_corrSeries{b},2))) ')'])
    % corrSeries = corrSeries(~isnan(corrSeries)); median(corrSeries)
    title({['corr coef'], ['smoo win ' num2str(BPwin)]},'Color',condicolors(conditions(b),:));
    ylim([-1 1]); yline(0,'color',deepyellow);
    xlim([0 intNum(b)]);

    subplot(4,12,24+b); % XcorrSeries
    plot(BPint_xcorrSeries{b},'k'); 
    ylabel('xcorr coef'); xlabel(['taps (sum ' num2str(round(sum(BPint_xcorrSeries{b},2))) ')'])
    % corrSeries = corrSeries(~isnan(corrSeries)); median(corrSeries)
    title({['xcorr lag 0'], ['smoo win ' num2str(BPwin)]},'Color',condicolors(conditions(b),:));
    ylim([-1 1]); yline(0,'color',deepyellow);
    xlim([0 intNum(b)]);

    subplot(4,12,36+b); % BPint difference
    plot(abs(intervals{b}(:,1)/2-intervals{b}(:,2)/2),'k'); 
    ylabel('difference (ms)'); 
    xlabel(['taps (sum ' num2str(round(sum(abs(intervals{b}(:,1)/2-intervals{b}(:,2)/2)))) ')'])
    ylim([-10 400]); yline(0,'color',deepyellow);
    xlim([0 intNum(b)]);
    title({[conditionNames{conditions(b)}],['BPwin ' num2str(BPwin)]},'Color',condicolors(conditions(b),:));
end
sgtitle(['session ' num2str(seed) ' ' sessionTypes{session} ' BPint corr xcorr diff']);


%% PLOT 12 500ms EEG 5 bands
clear
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816];
sessions={'synch','synco','synch','synco','synch','synco','synch','synco'};
EEGwin=0.5; % second

delta_L=cell(8,12);theta_L=cell(8,12);alpha_L=cell(8,12);beta_L=cell(8,12);gamma_L=cell(8,12);
delta_R=cell(8,12);theta_R=cell(8,12);alpha_R=cell(8,12);beta_R=cell(8,12);gamma_R=cell(8,12);
zEEG500_L=cell(8,12);zEEG500_R=cell(8,12);

cd /ssd/zhibin/1overf/
for s=1:8 % each session
    
    runid=num2str(seeds(s,:));
    clear dataL dataR
    load(['/ssd/zhibin/1overf/' runid '_2P/Cleaned_data/clean_' runid '.mat' ])
    zEEG_L={};zEEG_R={};
    for b=1:12 % each block
        
        zEEG_L{b}=zscore(dataL{b}(:,1:32),[],1);
        zEEG_R{b}=zscore(dataR{b}(:,1:32),[],1);

        for i=1:size(samples{b},1) % each matched interval
            
            zEEG500_L{s,b}{i}=zEEG_L{b}(samples{b}(i,1)-999:samples{b}(i,1),:); % start from first matched tap (if not enough samples before, then consider starting from 2nd tap)
            zEEG500_R{s,b}{i}=zEEG_R{b}(samples{b}(i,2)-999:samples{b}(i,2),:);
            [delta_L{s,b}(i,:), theta_L{s,b}(i,:), alpha_L{s,b}(i,:), beta_L{s,b}(i,:), gamma_L{s,b}(i,:)]...
                =sum5band(zEEG500_L{s,b}{i},sr,EEGwin);
            [delta_R{s,b}(i,:), theta_R{s,b}(i,:), alpha_R{s,b}(i,:), beta_R{s,b}(i,:), gamma_R{s,b}(i,:)]...
                =sum5band(zEEG500_R{s,b}{i},sr,EEGwin);
        end
    end
end

%% PLOT 13 corr between sum of EEG power and Hurst all sessions
zEEG500_L; zEEG500_R; % for all sessions from PLOT 12
delta_L;theta_L;alpha_L;beta_L;gamma_L;
delta_R;theta_R;alpha_R;beta_R;gamma_R;
% sum the band power in each channel in each block (8x12=96)
delta_L_sum=[];theta_L_sum=[];alpha_L_sum=[];beta_L_sum=[];gamma_L_sum=[];
delta_R_sum=[];theta_R_sum=[];alpha_R_sum=[];beta_R_sum=[];gamma_R_sum=[];
for s=1:8
    for b=1:12
        delta_L_sum(s,b,:)=sum([delta_L{s,b}]);
        theta_L_sum(s,b,:)=sum([theta_L{s,b}]);
        alpha_L_sum(s,b,:)=sum([alpha_L{s,b}]);
        beta_L_sum(s,b,:)=sum([beta_L{s,b}]);
        gamma_L_sum(s,b,:)=sum([gamma_L{s,b}]);
        delta_R_sum(s,b,:)=sum([delta_R{s,b}]);
        theta_R_sum(s,b,:)=sum([theta_R{s,b}]);
        alpha_R_sum(s,b,:)=sum([alpha_R{s,b}]);
        beta_R_sum(s,b,:)=sum([beta_R{s,b}]);
        gamma_R_sum(s,b,:)=sum([gamma_R{s,b}]);
    end
end
% squeeze into 32 vectors for corr with H
delta_L_chan=[];theta_L_chan=[];alpha_L_chan=[];beta_L_chan=[];gamma_L_chan=[];
delta_R_chan=[];theta_R_chan=[];alpha_R_chan=[];beta_R_chan=[];gamma_R_chan=[];
for c=1:32
    delta_L_chan(:,c)=reshape(delta_L_sum(:,:,c)',[],1); % 96x32 (each element from one block in time sequence) 
    theta_L_chan(:,c)=reshape(theta_L_sum(:,:,c)',[],1);
    alpha_L_chan(:,c)=reshape(alpha_L_sum(:,:,c)',[],1);
    beta_L_chan(:,c)=reshape(beta_L_sum(:,:,c)',[],1);
    gamma_L_chan(:,c)=reshape(gamma_L_sum(:,:,c)',[],1);
    delta_R_chan(:,c)=reshape(delta_R_sum(:,:,c)',[],1);
    theta_R_chan(:,c)=reshape(theta_R_sum(:,:,c)',[],1);
    alpha_R_chan(:,c)=reshape(alpha_R_sum(:,:,c)',[],1);
    beta_R_chan(:,c)=reshape(beta_R_sum(:,:,c)',[],1);
    gamma_R_chan(:,c)=reshape(gamma_R_sum(:,:,c)',[],1);
end
% Combine L and R
delta_LR_chan=[delta_L_chan;delta_R_chan];
theta_LR_chan=[theta_L_chan;theta_R_chan];
alpha_LR_chan=[alpha_L_chan;alpha_R_chan];
beta_LR_chan=[beta_L_chan;beta_R_chan];
gamma_LR_chan=[gamma_L_chan;gamma_R_chan];

% Organize H_all
H_all; % (8x12) for all sessions from PLOT 10-1 
H_all_L=squeeze(H_all(1,:,:));
H_all_R=squeeze(H_all(2,:,:));
% squeeze into 1 vector from the 96 blocks for each subject, for corr with pow in each chan
H_all_L=reshape(H_all_L',[],1);% 96x1 (each element from one block in time sequence) 
H_all_R=reshape(H_all_R',[],1);
% Combine L and R
H_all_LR=[H_all_L;H_all_R];

% Do the correlation between pow and H
for c=1:32
    delta_L_H_L_corr(c)=corr(delta_L_chan(:,c),H_all_L);
    theta_L_H_L_corr(c)=corr(theta_L_chan(:,c),H_all_L);
    alpha_L_H_L_corr(c)=corr(alpha_L_chan(:,c),H_all_L);
    beta_L_H_L_corr(c)=corr(beta_L_chan(:,c),H_all_L);
    gamma_L_H_L_corr(c)=corr(gamma_L_chan(:,c),H_all_L);
    delta_R_H_R_corr(c)=corr(delta_R_chan(:,c),H_all_R);
    theta_R_H_R_corr(c)=corr(theta_R_chan(:,c),H_all_R);
    alpha_R_H_R_corr(c)=corr(alpha_R_chan(:,c),H_all_R);
    beta_R_H_R_corr(c)=corr(beta_R_chan(:,c),H_all_R);
    gamma_R_H_R_corr(c)=corr(gamma_R_chan(:,c),H_all_R);
end
% Combine L and R for correlation;
for c=1:32
    delta_LR_H_LR_corr(c)=corr(delta_LR_chan(:,c),H_all_LR);
    theta_LR_H_LR_corr(c)=corr(theta_LR_chan(:,c),H_all_LR);
    alpha_LR_H_LR_corr(c)=corr(alpha_LR_chan(:,c),H_all_LR);
    beta_LR_H_LR_corr(c)=corr(beta_LR_chan(:,c),H_all_LR);
    gamma_LR_H_LR_corr(c)=corr(gamma_LR_chan(:,c),H_all_LR);
end
% 4 states
delta_LR_H_LR_4corr=zeros(4,32);
theta_LR_H_LR_4corr=zeros(4,32);
alpha_LR_H_LR_4corr=zeros(4,32);
beta_LR_H_LR_4corr=zeros(4,32);
gamma_LR_H_LR_4corr=zeros(4,32);
for s=1:4
    for c=1:32
        delta_LR_H_LR_4corr(s,c)=corr(delta_LR_chan(Inds4(:,s),c),H_all_LR(Inds4(:,s)));
        theta_LR_H_LR_4corr(s,c)=corr(theta_LR_chan(Inds4(:,s),c),H_all_LR(Inds4(:,s)));
        alpha_LR_H_LR_4corr(s,c)=corr(alpha_LR_chan(Inds4(:,s),c),H_all_LR(Inds4(:,s)));
        beta_LR_H_LR_4corr(s,c)=corr(beta_LR_chan(Inds4(:,s),c),H_all_LR(Inds4(:,s)));
        gamma_LR_H_LR_4corr(s,c)=corr(gamma_LR_chan(Inds4(:,s),c),H_all_LR(Inds4(:,s)));
    end
end

% Topoplots (5x1)
canvas(0.2,0.5);
cmin=-0.2;cmax=0.5;
subplot(5,2,1);topoplot(delta_L_H_L_corr,channels,'nosedir','+X');title('delta-L & H-L');colorbar;colormap('jet');clim([cmin cmax]);
subplot(5,2,3);topoplot(theta_L_H_L_corr,channels,'nosedir','+X');title('theta-L & H-L');colorbar;colormap('jet');clim([cmin cmax]);
subplot(5,2,5);topoplot(alpha_L_H_L_corr,channels,'nosedir','+X');title('alpha-L & H-L');colorbar;colormap('jet');clim([cmin cmax]);
subplot(5,2,7);topoplot(beta_L_H_L_corr,channels,'nosedir','+X');title('beta-L & H-L');colorbar;colormap('jet');clim([cmin cmax]);
subplot(5,2,9);topoplot(gamma_L_H_L_corr,channels,'nosedir','+X');title('gamma-L & H-L');colorbar;colormap('jet');clim([cmin cmax]);
subplot(5,2,2);topoplot(delta_R_H_R_corr,channels,'nosedir','+X');title('delta-R & H-R');colorbar;colormap('jet');clim([cmin cmax]);
subplot(5,2,4);topoplot(theta_R_H_R_corr,channels,'nosedir','+X');title('theta-R & H-R');colorbar;colormap('jet');clim([cmin cmax]);
subplot(5,2,6);topoplot(alpha_R_H_R_corr,channels,'nosedir','+X');title('alpha-R & H-R');colorbar;colormap('jet');clim([cmin cmax]);
subplot(5,2,8);topoplot(beta_R_H_R_corr,channels,'nosedir','+X');title('beta-R & H-R');colorbar;colormap('jet');clim([cmin cmax]);
subplot(5,2,10);topoplot(gamma_R_H_R_corr,channels,'nosedir','+X');title('gamma-R & H-R');colorbar;colormap('jet');clim([cmin cmax]);
sgtile('corr of EEG and interval_H')

% Combine L and R in topoplots (1x5)
canvas(0.5,0.1);
cmin=-0.4;cmax=0.4;
subplot(1,5,1);topoplot(delta_LR_H_LR_corr,channels,'nosedir','+X');title('delta & H');colorbar;colormap('jet');clim([cmin cmax]);
subplot(1,5,2);topoplot(theta_LR_H_LR_corr,channels,'nosedir','+X');title('theta & H');colorbar;colormap('jet');clim([cmin cmax]);
subplot(1,5,3);topoplot(alpha_LR_H_LR_corr,channels,'nosedir','+X');title('alpha & H');colorbar;colormap('jet');clim([cmin cmax]);
subplot(1,5,4);topoplot(beta_LR_H_LR_corr,channels,'nosedir','+X');title('beta & H');colorbar;colormap('jet');clim([cmin cmax]);
subplot(1,5,5);topoplot(gamma_LR_H_LR_corr,channels,'nosedir','+X');title('gamma & H');colorbar;colormap('jet');clim([cmin cmax]);
sgtitle('corr of EEG and interval-H')

% Combine L and R in 4 states(4x5)
canvas(0.5,0.4);
cmin=-0.4;cmax=0.4;
for s=1:4
    subplot(4,5,5*(s-1)+1);
    topoplot(delta_LR_H_LR_4corr(s,:),channels,'nosedir','+X');
    title([states4names{s} ': delta & H'],'Color',condicolors(s,:));
    colorbar;colormap('jet');clim([cmin cmax]);
    subplot(4,5,5*(s-1)+2);
    topoplot(theta_LR_H_LR_4corr(s,:),channels,'nosedir','+X');
    title([states4names{s} ': theta & H'],'Color',condicolors(s,:));
    colorbar;colormap('jet');clim([cmin cmax]);
    subplot(4,5,5*(s-1)+3);
    topoplot(alpha_LR_H_LR_4corr(s,:),channels,'nosedir','+X');
    title([states4names{s} ': alpha & H'],'Color',condicolors(s,:));
    colorbar;colormap('jet');clim([cmin cmax]);
    subplot(4,5,5*(s-1)+4);
    topoplot(beta_LR_H_LR_4corr(s,:),channels,'nosedir','+X');
    title([states4names{s} ': beta & H'],'Color',condicolors(s,:));
    colorbar;colormap('jet');clim([cmin cmax]);
    subplot(4,5,5*(s-1)+5);
    topoplot(gamma_LR_H_LR_4corr(s,:),channels,'nosedir','+X');
    title([states4names{s} ': gamma & H'],'Color',condicolors(s,:));
    colorbar;colormap('jet');clim([cmin cmax]);
end
sgtitle('4states: Corr of EEG & H-int');
%% PLOT 13-1 corr between H-power and H-interval all sessions & 4 states
% H derived from EEG power
delta_L;theta_L;alpha_L;beta_L;gamma_L;
delta_R;theta_R;alpha_R;beta_R;gamma_R;
H_delta_L=[];H_theat_L=[];H_alpha_L=[];H_beta_L=[];H_gamma_L=[];
H_delta_R=[];H_theat_R=[];H_alpha_R=[];H_beta_R=[];H_gamma_R=[];
for s=1:8
    for b=1:12
        for c=1:32
            [~,H_delta_L(s,b,c)]=DFA_main(delta_L{s,b}(:,c));
            [~,H_theta_L(s,b,c)]=DFA_main(theta_L{s,b}(:,c));
            [~,H_alpha_L(s,b,c)]=DFA_main(alpha_L{s,b}(:,c));
            [~,H_beta_L(s,b,c)]=DFA_main(beta_L{s,b}(:,c));
            [~,H_gamma_L(s,b,c)]=DFA_main(gamma_L{s,b}(:,c));
            [~,H_delta_R(s,b,c)]=DFA_main(delta_R{s,b}(:,c));
            [~,H_theta_R(s,b,c)]=DFA_main(theta_R{s,b}(:,c));
            [~,H_alpha_R(s,b,c)]=DFA_main(alpha_R{s,b}(:,c));
            [~,H_beta_R(s,b,c)]=DFA_main(beta_R{s,b}(:,c));
            [~,H_gamma_R(s,b,c)]=DFA_main(gamma_R{s,b}(:,c)); 
        end
    end
end
% squeeze into 32 vectors for corr with H
H_delta_L_chan=[];H_theat_L_chan=[];H_alpha_L_chan=[];H_beta_L_chan=[];H_gamma_L_chan=[];
H_delta_R_chan=[];H_theat_R_chan=[];H_alpha_R_chan=[];H_beta_R_chan=[];H_gamma_R_chan=[];
for c=1:32
    H_delta_L_chan(:,c) = reshape(H_delta_L(:,:,c)',[],1);
    H_theta_L_chan(:,c)  = reshape(H_theta_L(:,:,c)',[],1);
    H_alpha_L_chan(:,c)  = reshape(H_alpha_L(:,:,c)',[],1);
    H_beta_L_chan(:,c)   = reshape(H_beta_L(:,:,c)',[],1);
    H_gamma_L_chan(:,c)  = reshape(H_gamma_L(:,:,c)',[],1);
    H_delta_R_chan(:,c)  = reshape(H_delta_R(:,:,c)',[],1);
    H_theta_R_chan(:,c)  = reshape(H_theta_R(:,:,c)',[],1);
    H_alpha_R_chan(:,c)  = reshape(H_alpha_R(:,:,c)',[],1);
    H_beta_R_chan(:,c)   = reshape(H_beta_R(:,:,c)',[],1);
    H_gamma_R_chan(:,c)   = reshape(H_gamma_R(:,:,c)',[],1);
end
% combine L and R
H_delta_LR_chan=[H_delta_L_chan;H_delta_R_chan];
H_theta_LR_chan=[H_theta_L_chan;H_theta_R_chan];
H_alpha_LR_chan=[H_alpha_L_chan;H_alpha_R_chan];
H_beta_LR_chan=[H_beta_L_chan;H_beta_R_chan];
H_gamma_LR_chan=[H_gamma_L_chan;H_gamma_R_chan];

% check how many H-EEG > 0.5 (showing complexity)
% mean(H_delta_LR_chan)
% mean(H_theta_LR_chan)
% mean(H_alpha_LR_chan)
% mean(H_beta_LR_chan)
% mean(H_gamma_LR_chan)
canvas(0.5,0.1);
cmin=0.3;cmax=0.7;
subplot(1,5,1);topoplot(mean(H_delta_LR_chan),channels,'nosedir','+X');title('H-delta');colorbar;colormap('jet');clim([cmin cmax]);
subplot(1,5,2);topoplot(mean(H_theta_LR_chan),channels,'nosedir','+X');title('H-theta');colorbar;colormap('jet');clim([cmin cmax]);
subplot(1,5,3);topoplot(mean(H_alpha_LR_chan),channels,'nosedir','+X');title('H-alpha');colorbar;colormap('jet');clim([cmin cmax]);
subplot(1,5,4);topoplot(mean(H_beta_LR_chan),channels,'nosedir','+X');title('H-beta');colorbar;colormap('jet');clim([cmin cmax]);
subplot(1,5,5);topoplot(mean(H_gamma_LR_chan),channels,'nosedir','+X');title('H-gamma');colorbar;colormap('jet');clim([cmin cmax]);
sgtitle('mean H-EEG values')

% organize conditions in all sessions into a vector of 192
condition_all=[];
for s=1:8
    clear conditions
    runid=num2str(seeds(s,:));
    load(['/ssd/zhibin/1overf/' runid '_2P/Cleaned_data/clean_' runid '.mat' ]);
    condition_all(s,:)=conditions;
end
% reshape into a vector in time sequence
condition_all=reshape(condition_all',[],1);
% check how many H-EEG > 0.5 according to conditions (showing complexity)
states4names={'uncoupled','Leading','Following','mutual'};
% find the indices for each condition in L R conbined sequence (4 states)
uncoupleInd=repmat(find(condition_all==1),2,1);
leadingInd=[find(condition_all==2);96+find(condition_all==3)];
followingInd=[find(condition_all==3);96+find(condition_all==2)];
mutualInd=repmat(find(condition_all==4),2,1);
Inds4=[uncoupleInd leadingInd followingInd mutualInd];
canvas(0.5,0.5);
cmin=0.5;cmax=1;
for c=1:4
subplot(4,5,(c-1)*5+1);
topoplot(mean(H_delta_LR_chan(Inds4(:,c),:)),channels,'nosedir','+X');
title(['H-delta ' states4names{c}]);colorbar;colormap('jet');clim([cmin cmax]);
subplot(4,5,(c-1)*5+2);
topoplot(mean(H_theta_LR_chan(Inds4(:,c),:)),channels,'nosedir','+X');
title(['H-theta ' states4names{c}]);colorbar;colormap('jet');clim([cmin cmax]);
subplot(4,5,(c-1)*5+3);
topoplot(mean(H_alpha_LR_chan(Inds4(:,c),:)),channels,'nosedir','+X');
title(['H-alpha ' states4names{c}]);colorbar;colormap('jet');clim([cmin cmax]);
subplot(4,5,(c-1)*5+4);
topoplot(mean(H_beta_LR_chan(Inds4(:,c),:)),channels,'nosedir','+X');
title(['H-beta ' states4names{c}]);colorbar;colormap('jet');clim([cmin cmax]);
subplot(4,5,(c-1)*5+5);
topoplot(mean(H_gamma_LR_chan(Inds4(:,c),:)),channels,'nosedir','+X');
title(['H-gamma ' states4names{c}]);colorbar;colormap('jet');clim([cmin cmax]);
end
sgtitle('mean H-EEG values for each condition');

% Do the correlation between H-pow and H-interval
% Combine L and R for correlation;
for c=1:32
    H_delta_LR_H_LR_corr(c)=corr(H_delta_LR_chan(:,c),H_all_LR);
    H_theta_LR_H_LR_corr(c)=corr(H_theta_LR_chan(:,c),H_all_LR);
    H_alpha_LR_H_LR_corr(c)=corr(H_alpha_LR_chan(:,c),H_all_LR);
    H_beta_LR_H_LR_corr(c)=corr(H_beta_LR_chan(:,c),H_all_LR);
    H_gamma_LR_H_LR_corr(c)=corr(H_gamma_LR_chan(:,c),H_all_LR);
end

% Combine L and R in topoplots (1x5)
canvas(0.5,0.1);
cmin=-0.4;cmax=0.4;
subplot(1,5,1);topoplot(H_delta_LR_H_LR_corr,channels,'nosedir','+X');
title('H-delta & H-int');colorbar;colormap('jet');clim([cmin cmax]);
subplot(1,5,2);topoplot(H_theta_LR_H_LR_corr,channels,'nosedir','+X');
title('H-theta & H-int');colorbar;colormap('jet');clim([cmin cmax]);
subplot(1,5,3);topoplot(H_alpha_LR_H_LR_corr,channels,'nosedir','+X');
title('H-alpha & H-int');colorbar;colormap('jet');clim([cmin cmax]);
subplot(1,5,4);topoplot(H_beta_LR_H_LR_corr,channels,'nosedir','+X');
title('H-beta & H-int');colorbar;colormap('jet');clim([cmin cmax]);
subplot(1,5,5);topoplot(H_gamma_LR_H_LR_corr,channels,'nosedir','+X');
title('H-gamma & H-int');colorbar;colormap('jet');clim([cmin cmax]);
sgtitle('corr of H-EEG and H-interval')



%% PLOT 14 PLS regression
% % previously in sync.m
% open sync.m
% open syn_tryout.m
% % powforpls (28 matches x 3200 predictors) 50freqsx64chan=3200
% % cor(j).BP' (28 corr ceof x 1) 
% [R2,reg,ypred] = npls_pred(powforpls,cor(j).BP',1);
% weights = reshape(reg{1},50,64); % reg{1} (3200 x 1)
% R2 = R2;
% ypred = ypred;

pow5forpls=cat(3,H_delta_LR_chan,H_theta_LR_chan,H_alpha_LR_chan,H_beta_LR_chan,H_gamma_LR_chan);
pow5forpls=permute(pow5forpls,[3 2 1]);% freqxchanxtrials
pow5forpls=reshape(pow5forpls,5*32,192)'; % 192x160
% (ALL states: 5freq x 32chan = 160 predictors x 192 trials)
[R2,reg,ypred] = npls_pred(pow5forpls,H_all_LR,1);
weights = reshape(reg{1},5,32);
canvas(0.2,0.1)
imagesc(weights);colorbar; % caxis([-2 2]*10E-7);
yticks([1:5]);yticklabels({'delta','theta','alpha','beta','gamma'});
xticks([1:32]);xticklabels([labels]);xtickangle(90);
title('All 4 statues: PLS model');

% (Each of the 4 states: 5freq x 32 chan = 160 predictors x 48 trials)
canvas(0.4,0.2);
for c=1:4 % four states
    R2=[];reg=[];ypred=[];
    [R2,reg,ypred] = npls_pred(pow5forpls(Inds4(:,c),:),H_all_LR(Inds4(:,c)),1);
    plsmodel(c).weights=reshape(reg{1},5,32);
    plsmodel(c).R2 = R2;
    plsmodel(c).ypred = ypred;
    subplot(2,2,c);
    imagesc(plsmodel(c).weights);colorbar; 
    yticks([1:5]);yticklabels({'delta','theta','alpha','beta','gamma'});
    xticks([1:32]);xticklabels([labels]);xtickangle(90);
    clim([-0.1 0.1]); % caxis([-0.1 0.1]);
    title([states4names{c} ': PLS model'],'Color',condicolors(c,:));
end
sgtitle('PLS model: H-EEG(-500ms) -> H-int')


%% try out DFA on EEG
ans(:,15)
[D,Alpha1,n,F_n,FitValues]=DFA_main(ans(:,15));




%% saves all variables from the current workspace
tic
save([num2str(seed) 'workspace.mat']);
toc
% takes about 11 min


%% Get variable names from mat file (faster)
load("20220721.mat")

tic
load([num2str(seed) 'workspace.mat']);
toc
% takes about 3 min
varlist1 = who();

% varlist1 = who('-file','20220811.mat');

tic
% matObj1 = matfile('20220721workspace.mat');
% varlist1 = who(matObj1);
toc

varlist2 = who();

varlist_diff = setdiff(varlist2,varlist1);

% clear varlist1

% append new variables to the old mat files
% tic
% varlist_diff_string=string(varlist_diff);
% for i=1:length(varlist_diff_string)
%     save('20220811workspace.mat',varlist_diff_string(i),'-append');
% end
% toc

tic 
varlist_diff_string=string(varlist_diff);
save([num2str(seed) 'workspace.mat'],varlist_diff_string{:},'-append');
toc
% takes about 1.5 min

% clear all new added variables
clear(varlist_diff{:})

%% moving files 

seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816];
cd /ssd/zhibin/1overf/
for s=1:8
    mkdir(['/ssd/zhibin/1overf/' num2str(seeds(s,:)) '_2P/Cleaned_data/Plots'])
    cd(['/ssd/zhibin/1overf/' num2str(seeds(s,:)) '_2P/Cleaned_data'])
    movefile cleanup*.png Plots
    cd /ssd/zhibin/1overf/
end

%%
