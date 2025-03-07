%% SECT 1: Auto&Xcorr DFA GC for each trial

cd /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/ARFIMA/Sim_plots
cd /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing
addpath /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing
addpath /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/DFA

cd /ssd/zhibin/1overf/20220713_2P/Segmented_data
cd /ssd/zhibin/1overf/20220721_2P/Segmented_data
cd /ssd/zhibin/1overf/20221003_2P/Segmented_data
cd /ssd/zhibin/1overf/2022100401_2P/Segmented_data
cd /ssd/zhibin/1overf/2022100402_2P/Segmented_data
cd /ssd/zhibin/1overf/20221005_2P/Segmented_data

load('data_variables20220713.mat')
load('data_variables20220721.mat')

% d estimate and removal
% http://www.lucafaes.net/LMSE-MSE_VARFI.html
addpath /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/MSE-VARFI
% Granger Causality
addpath(genpath('/home/zhibinz2/Documents/GitHub/MVGC1'));
cd /home/zhibinz2/Documents/GitHub/MVGC1
run startup
addpath /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/VAR_Granger

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
%     [F,pval,sig] = mystatespace(y1(i).BPint,y2(i).BPint); % Statespace method
    [F,pval,sig] = myautocov(y1(i).BPint,y2(i).BPint); % Autocov method
    y12(i).F=F;
    y12(i).pval=pval;
    y12(i).sig=sig;
end

%% PLOT 1: Auto&Xcorr DFA GC for each trial
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

%% SECT 2: Auto&Xcorr DFA GC PSA for each trial
cd /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/ARFIMA/Sim_plots
cd /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing
addpath /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing
addpath /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/DFA

cd /ssd/zhibin/1overf/20220713_2P/Segmented_data
cd /ssd/zhibin/1overf/20220721_2P/Segmented_data
cd /ssd/zhibin/1overf/20220804_2P/Segmented_data
load('data_variables20220713.mat')
load('data_variables20220721.mat')

% d estimate and removal
addpath /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/MSE-VARFI
% Granger Causality
addpath(genpath('/home/zhibinz2/Documents/GitHub/MVGC1'));
% rmpath('/home/zhibinz2/Documents/GitHub/MVGC1'); % remove duplecate startup file from path

% Structure viewer not showing properly after running the following 2 lines
% No need to run the following 2 lines
% cd /home/zhibinz2/Documents/GitHub/MVGC1
% run startup

% remove this path to recover structure viewing
% rmpath(genpath('/home/zhibinz2/Documents/GitHub/MVGC1')); 

addpath /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/VAR_Granger

% parameters for powers
% open /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing/spectra2
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
    
%% PLOT 2: Auto&Xcorr DFA GC PSA for each trial
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

%% PLOT 2-1: (Color Scheme) 
% cd /ssd/zhibin/1overf/20220804_2P/Segmented_data/Plots/Plot2_1/
% mkdir Plot2_1
% cd Plot2_1
% Plots - color scheme
red   = [1 0 0];
pink  = [1 0.65 0.75];
black = [0 0 0];
white = [1 1 1];
blue  = [0 0 1];
mediumblue = [0 0.4 0.7];
green = [0 1 0];
darkgreen = [0 0.5 0];
grey  = [0.5 0.5 0.5];
yellow  = [1 1 0];
deepyellow  = [1 0.8 0.2];
gold = [212/255 175/255 55/255];
brown = [150/255 75/255 0];
megenta = [1 0 1];% fill([0 1 1 0],[0 0 1 1],megenta)
cyan = [0 1 1]; % fill([0 1 1 0],[0 0 1 1],cc)
purple = [0.6 0.1 0.9];
% https://www.mathworks.com/help/matlab/creating_plots/specify-plot-colors.html
matlab_blue=[0 0.4470 0.7410];
matlab_orange=[0.8500 0.3250 0.0980];
matlab_gold=[0.9290 0.6940 0.1250];
matlab_purple=[0.4940 0.1840 0.5560];
matlab_green=[0.4660 0.6740 0.1880];
matlab_cyan=[0.3010 0.7450 0.9330];
matlab_red=[0.6350 0.0780 0.1840];
% combine colors
condicolors=[darkgreen;red;blue;megenta;purple;purple];
dire3colors=[darkgreen;brown;megenta];
% syn2colors=[matlab_blue;matlab_orange];
syn2colors=[darkgreen;pink];
HNLcolors = [darkgreen; deepyellow; pink];

% test color
showcolor=black;
imagesc(cat(3,showcolor(1),showcolor(2),showcolor(3)));

% colormap
cd /home/zhibinz2/Documents/GitHub/matlab-archive/hnlcode/common/gen_code/color
addpath /home/zhibinz2/Documents/GitHub/matlab-archive/hnlcode/common/gen_code/color
cd 
hnc = hotncold(100);
hottermap=hotter(100);
wnc=warmncool(100);

% test colormap
mycolormap=hottermap;
mymat=repmat(-0.99:0.02:0.99,100,1); % transitional color
% mymat=(rand(100)-0.5)*2; % random from -1 to 1
imagesc(mymat);colorbar;
colormap(mycolormap);

%% PLOT 2-1: Auto&Xcorr DFA GC PSA EEG Spectra for each trial
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


%% SECT 3 (unfinished): slinding window: DFA Autocov GC for each trial
cd /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/ARFIMA/Sim_plots
cd /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing
addpath /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing
addpath /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/DFA

cd /ssd/zhibin/1overf/20220713_2P/Segmented_data
load('data_variables20220713.mat')

cd /ssd/zhibin/1overf/20220721_2P/Segmented_data
load('data_variables20220721.mat')


% d estimate and removal
addpath /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/MSE-VARFI
% Granger Causality
addpath(genpath('/home/zhibinz2/Documents/GitHub/MVGC1'));
% cd /home/zhibinz2/Documents/GitHub/MVGC1
% run startup
addpath /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/VAR_Granger

% remove this path to recover structure viewing
% rmpath(genpath('/home/zhibinz2/Documents/GitHub/MVGC1')); 

% parameters for powers
% open /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing/spectra2
maxfreq=50; win=10; % 10 Second with df of 0.1Hz

% Sliding win on BP01 for DFA and GC
y1;y2;
% Calculation of sliding window
% window of time (ms) for 100 taps
winsize = 80*sr; % 80 Second with df of 0.0125 Hz, about 104 taps
% winsize = 10*sr; % 10 Second with df of 0.1Hz
% winsize = 2*sr; % 2 Second with df of 0.5 Hz
overlapsize = round(winsize*0.3); % number of samples for each overlapping window of3.3 seconds

% Matrix Structure (no heading in variable viewer)
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
rmpath(genpath('/home/zhibinz2/Documents/GitHub/MVGC1')); 
open ys1;
open ys12;
open ys2;
addpath(genpath('/home/zhibinz2/Documents/GitHub/MVGC1')); 

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
rmpath(genpath('/home/zhibinz2/Documents/GitHub/MVGC1')); 
open ys1;
open ys12;
open ys2;
addpath(genpath('/home/zhibinz2/Documents/GitHub/MVGC1')); 


% sum EEG power
%  for m = 1:50 % 50 frequencies
%     freqs = (m-1)*10+1:m*10; % downsample freq
%     power_L = sum(abs(tempf_L(freqs,:)),1);
%     power_R = sum(abs(tempf_R(freqs,:)),1);
%     pow_L(j).EEG(m,:,k) = power_L; % freq x chan x time
%     pow_R(j).EEG(m,:,k) = power_R;
%  end

    
%% PLOT 3 (unfinished): slinding window: DFA Autocov GC for each trial
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

%% SECT 4: 500ms before BP
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
eegbandNames={'Delta','Theta','Alpha','Beta','Gamma'};
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

%% SECT 5: GC between BPint and EEGpow
GC1=[];GC2=[];
eegbandNames={'Delta','Theta','Alpha','Beta','Gamma'};
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
            topoplot(squeeze(GCmat(i,seleted_band,:)),chaninfo,'nosedir','+X','style','map');
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


%% SECT 6 DFA H values for all sessions (unmatched/matched int) (sorted order)
clear
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);
sessions={'synch','synco','synch','synco','synch','synco','synch','synco','synch','synco',...
    'synch','synco'};

H=zeros(2,12,numSes);
cd /ssd/zhibin/1overf/
tic
for r=1:numSes;
    clear intervals conditions sortorders
    runid = num2str(seeds(r,:));
    % load unmatched int.
    % path = [runid '_2P/Segmented_data/'];
    % load([path  'data_variables' runid '.mat']);
    % load matched int.
    path = [runid '_2P/Cleaned_data/'];
    load([path  'clean_' runid '.mat']);
    % sort order and plot results
    % [x,sortorder]=sort(condiSeq);
    [x,sortorder]=sort(conditions);
    for j = 1:12
        % For unmatched int.
        % clear BPL BPR sampL sampR intL intR intL_dmean intR_dmean
	    % BPL = [BP(1).BP{sortorder(j)}];
	    % BPR = [BP(2).BP{sortorder(j)}];
        % sampL = find(BPL==1);
        % sampR = find(BPR==1);
        % intL=diff(sampL);
        % intR=diff(sampR);
        % For matched int.
        clear intL intR intL_dmean intR_dmean
        intL = intervals{sortorder(j)}(:,1);
        intR = intervals{sortorder(j)}(:,2);
        % remove the mean
        intL_dmean=intL-mean(intL);
        intR_dmean=intR-mean(intR);
        [~,H(1,j,r)]=DFA_main(intL_dmean);
        [~,H(2,j,r)]=DFA_main(intR_dmean);
    end
end
toc % 86 s

%% PLOT 6 DFA H values for all sessions (matched/unmatched int)(sorted order)
% print all H to a table
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
H_L_synch5=H(1,:,9)';H_R_synch5=H(2,:,9)';
H_L_synco5=H(1,:,10)';H_R_synco5=H(2,:,10)';
H_L_synch6=H(1,:,11)';H_R_synch6=H(2,:,11)';
H_L_synco6=H(1,:,12)';H_R_synco6=H(2,:,12)';
Table_H=table(Trials, H_L_synch1, H_R_synch1, H_L_synco1, H_R_synco1,...
    H_L_synch2, H_R_synch2, H_L_synco2, H_R_synco2,...
    H_L_synch3, H_R_synch3, H_L_synco3, H_R_synco3,...
    H_L_synch4, H_R_synch4, H_L_synco4, H_R_synco4,...
    H_L_synch5, H_R_synch5, H_L_synco5, H_R_synco5,...
    H_L_synch6, H_R_synch6, H_L_synco6, H_R_synco6);
% plot all
canvas(0.5, 0.5);
for r=1:numSes
    subplot(numSes,1,r);
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
sgtitle('H from DFA in all sessions (matched int) ^{ *PLOT 6}')
% sgtitle('H from DFA in all sessions (unmatched int) ^{ *PLOT 6-1}')

% average all trials of the same condition in synch/0 (2 subplots)
H_L_uncouple_synch = reshape(H(1,[1:3],[1:2:11]),1,[])  ;% L,uncouple,synch
H_L_Llead_synch    = reshape(H(1,[4:6],[1:2:11]),1,[])  ;% L,L lead,synch
H_L_Rlead_synch    = reshape(H(1,[7:9],[1:2:11]),1,[])  ;% L,R lead,synch
H_L_mutual_synch   = reshape(H(1,[10:12],[1:2:11]),1,[]);% L, mutual,synch
H_R_uncouple_synch = reshape(H(2,[1:3],[1:2:11]),1,[])  ;% R,uncouple,synch
H_R_Llead_synch    = reshape(H(2,[4:6],[1:2:11]),1,[])  ;% R,L lead,synch
H_R_Rlead_synch    = reshape(H(2,[7:9],[1:2:11]),1,[])  ;% R,R lead,synch
H_R_mutual_synch   = reshape(H(2,[10:12],[1:2:11]),1,[]);% R, mutual,synch
H_LR_synch_mean = [mean(H_L_uncouple_synch) mean(H_R_uncouple_synch);...
    mean(H_L_Llead_synch) mean(H_R_Llead_synch);...
    mean(H_L_Rlead_synch) mean(H_R_Rlead_synch);...
    mean(H_L_mutual_synch) mean(H_R_mutual_synch)];

H_L_uncouple_synco = reshape(H(1,[1:3],[2:2:12]),1,[])  ;% R,uncouple,synco
H_L_Llead_synco    = reshape(H(1,[4:6],[2:2:12]),1,[])  ;% R,L lead,synco
H_L_Rlead_synco    = reshape(H(1,[7:9],[2:2:12]),1,[])  ;% R,R lead,synco
H_L_mutual_synco   = reshape(H(1,[10:12],[2:2:12]),1,[]);% R, mutual,synco
H_R_uncouple_synco = reshape(H(2,[1:3],[2:2:12]),1,[])  ;% R,uncouple,synco
H_R_Llead_synco    = reshape(H(2,[4:6],[2:2:12]),1,[])  ;% R,L lead,synco
H_R_Rlead_synco    = reshape(H(2,[7:9],[2:2:12]),1,[])  ;% R,R lead,synco
H_R_mutual_synco   = reshape(H(2,[10:12],[2:2:12]),1,[]);% R, mutual,synco
H_LR_synco_mean = [mean(H_L_uncouple_synco) mean(H_R_uncouple_synco);...
    mean(H_L_Llead_synco) mean(H_R_Llead_synco);...
    mean(H_L_Rlead_synco) mean(H_R_Rlead_synco);...
    mean(H_L_mutual_synco) mean(H_R_mutual_synco)];
% plot with bar groups for L/R
figure;
subplot(2,1,1);
b=bar(H_LR_synch_mean);b(1).FaceColor = [1 0 0]; b(2).FaceColor = [0 0 1];
xticks(1:4);xticklabels({ 'Uncouple','L-lead','R-lead','Mutual'});
ylabel('H');
hold on;
yline(1,'color',deepyellow);
% yline(0.5,'color',deepyellow);
ylim([0.5 0.9]);
title('mean H in all synch sessions (matched int) ^{ *PLOT 6}')
% title('mean H in all synco sessions (unmatched int) ^{ *PLOT 6-1}')
subplot(2,1,2);
b=bar(H_LR_synco_mean);b(1).FaceColor = [1 0 0]; b(2).FaceColor = [0 0 1];
xticks(1:4);xticklabels({ 'Uncouple','L-lead','R-lead','Mutual'});
ylabel('H');
hold on;
yline(1,'color',deepyellow);
% yline(0.5,'color',deepyellow);
ylim([0.5 0.9]);
title('mean H in all synco sessions (matched int) ^{ *PLOT 6}')
% title('mean H in all synco sessions (unmatched int) ^{ *PLOT 6-1}')

% combine synch/o 
H_LR_syn_mean=[mean([H_L_uncouple_synch H_L_uncouple_synco H_R_uncouple_synch H_R_uncouple_synco]);...
    mean([H_L_Llead_synch H_L_Llead_synco H_R_Rlead_synch H_R_Rlead_synco]);... 
    mean([H_R_Llead_synch H_R_Llead_synco H_L_Rlead_synch H_L_Rlead_synco]);...
    mean([H_L_mutual_synch H_L_mutual_synco H_R_mutual_synch H_R_mutual_synco])];
figure;
b=bar(H_LR_syn_mean);
xticks(1:4);xticklabels({'Uncouple','Leading','Following','Mutual'});
ylabel('H');
hold on;
ylim([0.6 0.8])
% yline(0.5,'color',deepyellow);
set(gcf,'color','w'); % set background white for copying in ubuntu
title('mean H in all sessions (matched int) ^{ *PLOT 6}');
% title('mean H in all sessions (unmatched int) ^{ *PLOT 6-1}');

% separate synch and synco with error bar
H_LR_synch_mean=[mean([H_L_uncouple_synch H_R_uncouple_synch]);...
    mean([H_L_Llead_synch H_R_Rlead_synch]);... 
    mean([H_R_Llead_synch H_L_Rlead_synch]);...
    mean([H_L_mutual_synch H_R_mutual_synch])];
H_LR_synco_mean=[mean([H_L_uncouple_synco H_R_uncouple_synco]);...
    mean([H_L_Llead_synco H_R_Rlead_synco]);... 
    mean([H_R_Llead_synco H_L_Rlead_synco]);...
    mean([H_L_mutual_synco H_R_mutual_synco])];
H_LR_synch_std=[std([H_L_uncouple_synch H_R_uncouple_synch]);...
    std([H_L_Llead_synch H_R_Rlead_synch]);... 
    std([H_R_Llead_synch H_L_Rlead_synch]);...
    std([H_L_mutual_synch H_R_mutual_synch])];
H_LR_synco_std=[std([H_L_uncouple_synco H_R_uncouple_synco]);...
    std([H_L_Llead_synco H_R_Rlead_synco]);... 
    std([H_R_Llead_synco H_L_Rlead_synco]);...
    std([H_L_mutual_synco H_R_mutual_synco])];
% grouped barplot with errorbar (4 subplots) and statistical tests
% (Publised)
% https://www.mathworks.com/matlabcentral/answers/102220-how-do-i-place-errorbars-on-my-grouped-bar-graph-using-function-errorbar-in-matlab
canvas(0.23, 0.4);
model_series = [H_LR_synch_mean H_LR_synco_mean];
model_error = [H_LR_synch_std H_LR_synco_std];
b = bar(model_series, 'grouped');
b(1).FaceColor=darkgreen;b(2).FaceColor=pink;
hold on;
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(model_series);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars 
errorbar(x',model_series,model_error,'k','linestyle','none');
hold off
xticks(1:4);xticklabels({'Independent','Leader','Follower','Bidirectional'});
xl = get(gca,'XTickLabel');  
set(gca,'XTickLabel',xl,'fontsize',17,'FontWeight','bold')
xlim([0.5 4.5]);
ylabel('DFA Hurst Exponent','FontSize',17); ylim([0.4 1]);
set(gcf,'color','w'); % set background white for copying in ubuntu
delete(findall(gcf,'type','annotation'))
sg=annotation('textbox',[0.05 0.01 0.5 0.07],'string',...
    {['mean H(matched int) ^{ *PLOT 6}' char(datetime('now'))]})
sg.Rotation=90
yline(0.5,'color', deepyellow)
lg=legend({'Synch','Synco'},'location','north');lg.FontSize-17;
% statistical test (oneway-ANOVA) 
[p,tbl,stats] = anova1([...
    [H_L_uncouple_synch H_R_uncouple_synch]' ...
    [H_L_Llead_synch H_R_Rlead_synch]' ... 
    [H_R_Llead_synch H_L_Rlead_synch]' ...
    [H_L_mutual_synch H_R_mutual_synch]' ...
    [H_L_uncouple_synco H_R_uncouple_synco]' ...
    [H_L_Llead_synco H_R_Rlead_synco]' ... 
    [H_R_Llead_synco H_L_Rlead_synco]' ...
    [H_L_mutual_synco H_R_mutual_synco]'...
    ]); % one-way ANOVA (work)
[p,tbl,stats] = anova2([...
    [[H_L_uncouple_synch H_R_uncouple_synch]' [H_L_uncouple_synco H_R_uncouple_synco]'];...
    [[H_L_Llead_synch H_R_Rlead_synch]' [H_L_Llead_synco H_R_Rlead_synco]'];... 
    [[H_R_Llead_synch H_L_Rlead_synch]' [H_R_Llead_synco H_L_Rlead_synco]'];...
    [[H_L_mutual_synch H_R_mutual_synch]' [H_L_mutual_synco H_R_mutual_synco]']...
    ],36); % two-way ANOVA (work)

% Bonferroni correction
[results,means,~,gnames] = multcompare(stats,"CriticalValueType","bonferroni"); % multi-comparison
tbl = array2table([results],"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])
tb2 = array2table([means],"VariableNames", ["Mean","Standard Error"])
p_bonf=results(:,6);
% statistical test (Two sample t test)
[h,p(1),ci,stats] = ttest2([H_L_uncouple_synch H_R_uncouple_synch]',[H_L_uncouple_synco H_R_uncouple_synco]') % Uncouple
[h,p(2),ci,stats] = ttest2([H_L_Llead_synch H_R_Rlead_synch]',[H_L_Llead_synco H_R_Rlead_synco]') % Leading
[h,p(3),ci,stats] = ttest2([H_R_Llead_synch H_L_Rlead_synch]',[H_R_Llead_synco H_L_Rlead_synco]','Alpha',0.01) % Following
[h,p(4),ci,stats] = ttest2([H_L_mutual_synch H_R_mutual_synch]',[H_L_mutual_synco H_R_mutual_synco]') % Mutual
sprintf('P values for each group: %0.3f\n',p)
[h,p,ci,stats] = ttest2([H_L_mutual_synch H_R_mutual_synch H_L_Llead_synch H_R_Rlead_synch ...
    H_R_Llead_synch H_L_Rlead_synch H_L_mutual_synch H_R_mutual_synch]',...
    [H_L_mutual_synco H_R_mutual_synco H_L_Llead_synco H_R_Rlead_synco ...
    H_R_Llead_synco H_L_Rlead_synco H_L_mutual_synco H_R_mutual_synco]') % Synch vs Synco
% plot lines with errorbar 
canvas(0.23, 0.4);
for i=1:2;
    hold on;
    b=[];
    if i==1; b=errorbar(1:4,H_LR_synch_mean,H_LR_synch_std),'go'; %title('mean H in synch');
    else; b=errorbar(1:4,H_LR_synco_mean,H_LR_synch_std),'ro'; %title('mean H in synco');
    end
    xticks(1:4);xticklabels({'Uncouple','Leading','Following','Mutual'});xlim([0 5])
    ylabel('H');
    ylim([0.4 1])
end
legend('Synch','Synco')
set(gcf,'color','w'); % set background white for copying in ubuntu
% sg=sgtitle({['mean H(matched int) ^{ *PLOT 6}'],char(datetime('now'))});
% title('mean H in all sessions (unmatched int) ^{ *PLOT 6-1}');
delete(findall(gcf,'type','annotation'))
sg=annotation('textbox',[0.05 0.01 0.5 0.07],'string',...
    {['mean H(matched int) ^{ *PLOT 6}' char(datetime('now'))]})
sg.Rotation=90

% Two violin plots
addpath /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/violin
canvas(0.4,0.3);
for violin2=1;
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
title('H synch (L - red; R - blue) (matched int) ^{ *PLOT 6}');
title('H synch (L - red; R - blue) (unmatched int) ^{ *PLOT 6-1}');
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
title('H synco (L - red; R - blue) (matched int) ^{ *PLOT 6}');
title('H synco (L - red; R - blue) (unmatched int) ^{ *PLOT 6-1}');
grid on;
end

% combine synch and synco in violin plot
figure;
for violin1=1;
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
title('mean H in all sessions (matched int) ^{ *PLOT 6}');
% title('mean H in all sessions (unmatched int) ^{ *PLOT 6-1}');
grid on;
end

%% SECT 7 Autocorr summary in all sessions (unmatched/matched int) (sorted order)
clear
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);
sessions={'synch','synco','synch','synco','synch','synco','synch','synco','synch','synco',...
    'synch','synco'};

Autocorrs20lags=nan(2,12,numSes,20);%save the values of all 20 lags
Autocorrs=zeros(2,12,numSes);% sum up the 20 lags
cd /ssd/zhibin/1overf/
tic
for r=1:numSes;
    clear intervals conditions sortorders
    runid = num2str(seeds(r,:));
    % load unmatched int.
    % path = [runid '_2P/Segmented_data/'];
    % load([path  'data_variables' runid '.mat']);
    % load matched int.
    path = [runid '_2P/Cleaned_data/'];
    load([path  'clean_' runid '.mat']);
    load([path  'd_removal' runid '.mat']);
    % sort order and plot results
    % [x,sortorder]=sort(condiSeq);
    [x,sortorder]=sort(conditions);
    for j = 1:12 
        % For unmatched int
        % clear BPL BPR sampL sampR intL intR intL_dmean intR_dmean
	    % BPL = [BP(1).BP{sortorder(j)}];
	    % BPR = [BP(2).BP{sortorder(j)}];
        % sampL = find(BPL==1);
        % sampR = find(BPR==1);
        % intL=diff(sampL);
        % intR=diff(sampR);
        % For matched int.
        clear intL intR intL_dmean intR_dmean
        intL = intervals{sortorder(j)}(:,1);
        intR = intervals{sortorder(j)}(:,2);
        % remove the mean
        intL_dmean=intL-mean(intL);
        intR_dmean=intR-mean(intR);
        acf=[]; [acf] = autocorr(intL_dmean,20); 
        Autocorrs(1,j,r)=sum(acf(2:end));
        Autocorrs20lags(1,j,r,1:20)=acf(2:end);
        acf=[]; [acf] = autocorr(intR_dmean,20); 
        Autocorrs(2,j,r)=sum(acf(2:end));
        Autocorrs20lags(2,j,r,1:20)=acf(2:end);
    end
end
toc
% 88 sec


% print Autocorrs(sum) to a table (show each session)
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
Autocorrs_L_synch5=Autocorrs(1,:,9)';Autocorrs_R_synch5=Autocorrs(2,:,9)';
Autocorrs_L_synco5=Autocorrs(1,:,10)';Autocorrs_R_synco5=Autocorrs(2,:,10)';
Autocorrs_L_synch6=Autocorrs(1,:,11)';Autocorrs_R_synch6=Autocorrs(2,:,11)';
Autocorrs_L_synco6=Autocorrs(1,:,12)';Autocorrs_R_synco6=Autocorrs(2,:,12)';
Table_Autocorrs=table(Trials, Autocorrs_L_synch1, Autocorrs_R_synch1, Autocorrs_L_synco1, Autocorrs_R_synco1,...
    Autocorrs_L_synch2, Autocorrs_R_synch2, Autocorrs_L_synco2, Autocorrs_R_synco2,...
    Autocorrs_L_synch3, Autocorrs_R_synch3, Autocorrs_L_synco3, Autocorrs_R_synco3,...
    Autocorrs_L_synch4, Autocorrs_R_synch4, Autocorrs_L_synco4, Autocorrs_R_synco4,...
    Autocorrs_L_synch5, Autocorrs_R_synch5, Autocorrs_L_synco5, Autocorrs_R_synco5,...
    Autocorrs_L_synch6, Autocorrs_R_synch6, Autocorrs_L_synco6, Autocorrs_R_synco6);

%% PLOT 7 Autocorr summary in all sessions (unmatched/matched int) (sorted order)

% run similar plots in previous section (PLOT 6)
% plot each session
canvas(0.5, 0.5);
for r=1:numSes
    subplot(numSes,1,r);
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
% sgtitle('Sum of Autocorr (20 lags) in all sessions (unmatched int) ^{ *PLOT7}')
sgtitle('Sum of Autocorr (20 lags) in all sessions (matched int) ^{ *PLOT7}')

% average all trials of the same condition (synch vs synco)
% Autocorrs (sum)
Autocorrs_L_uncouple_synch = reshape(Autocorrs(1,[1:3],[1:2:11]),1,[])  % L,uncouple,synch
Autocorrs_L_Llead_synch    = reshape(Autocorrs(1,[4:6],[1:2:11]),1,[])  % L,L lead,synch
Autocorrs_L_Rlead_synch    = reshape(Autocorrs(1,[7:9],[1:2:11]),1,[])  % L,R lead,synch
Autocorrs_L_mutual_synch   = reshape(Autocorrs(1,[10:12],[1:2:11]),1,[])% L, mutual,synch
Autocorrs_R_uncouple_synch = reshape(Autocorrs(2,[1:3],[1:2:11]),1,[])  % R,uncouple,synch
Autocorrs_R_Llead_synch    = reshape(Autocorrs(2,[4:6],[1:2:11]),1,[])  % R,L lead,synch
Autocorrs_R_Rlead_synch    = reshape(Autocorrs(2,[7:9],[1:2:11]),1,[])  % R,R lead,synch
Autocorrs_R_mutual_synch   = reshape(Autocorrs(2,[10:12],[1:2:11]),1,[])% R, mutual,synch
Autocorrs_LR_synch_mean = [mean(Autocorrs_L_uncouple_synch) mean(Autocorrs_R_uncouple_synch);...
    mean(Autocorrs_L_Llead_synch) mean(Autocorrs_R_Llead_synch);...
    mean(Autocorrs_L_Rlead_synch) mean(Autocorrs_R_Rlead_synch);...
    mean(Autocorrs_L_mutual_synch) mean(Autocorrs_R_mutual_synch)];

Autocorrs_L_uncouple_synco = reshape(Autocorrs(1,[1:3],[2:2:12]),1,[])  % R,uncouple,synco
Autocorrs_L_Llead_synco    = reshape(Autocorrs(1,[4:6],[2:2:12]),1,[])  % R,L lead,synco
Autocorrs_L_Rlead_synco    = reshape(Autocorrs(1,[7:9],[2:2:12]),1,[])  % R,R lead,synco
Autocorrs_L_mutual_synco   = reshape(Autocorrs(1,[10:12],[2:2:12]),1,[])% R, mutual,synco
Autocorrs_R_uncouple_synco = reshape(Autocorrs(2,[1:3],[2:2:12]),1,[])  % R,uncouple,synco
Autocorrs_R_Llead_synco    = reshape(Autocorrs(2,[4:6],[2:2:12]),1,[])  % R,L lead,synco
Autocorrs_R_Rlead_synco    = reshape(Autocorrs(2,[7:9],[2:2:12]),1,[]) % R,R lead,synco
Autocorrs_R_mutual_synco   = reshape(Autocorrs(2,[10:12],[2:2:12]),1,[])% R, mutual,synco
Autocorrs_LR_synco_mean = [mean(Autocorrs_L_uncouple_synco) mean(Autocorrs_R_uncouple_synco);...
    mean(Autocorrs_L_Llead_synco) mean(Autocorrs_R_Llead_synco);...
    mean(Autocorrs_L_Rlead_synco) mean(Autocorrs_R_Rlead_synco);...
    mean(Autocorrs_L_mutual_synco) mean(Autocorrs_R_mutual_synco)];

figure;
subplot(2,1,1);
b=bar(Autocorrs_LR_synch_mean);b(1).FaceColor = [1 0 0]; b(2).FaceColor = [0 0 1];
xticks(1:4);xticklabels({ 'Uncouple','L-lead','R-lead','Mutual'});
ylabel('\rho(k)');
% title('mean autocorr sum for 20 lags in all synch sessions (unmatched int) ^{ *PLOT7}')
title('Synchronization')
subplot(2,1,2);
b=bar(Autocorrs_LR_synco_mean);b(1).FaceColor = [1 0 0]; b(2).FaceColor = [0 0 1];
xticks(1:4);xticklabels({ 'Uncouple','L-lead','R-lead','Mutual'});
ylabel('\rho(k)');
% title('mean autocorr sum for 20 lags in all synco sessions (unmatched int) ^{ *PLOT7}')
title('Syncopation')
sgtitle('mean autocorr sum for 20 lags in all synco sessions (matched int) ^{ *PLOT7}')

% average all trials of the same condition (synch vs synco)
% Autocorrs20lags (all 20 lags)
Autocorrs20lags_L_uncouple_synch = squeeze(mean(mean(squeeze(Autocorrs20lags(1,[1:3],[1:2:11],:)),2),1));  % L,uncouple,synch
Autocorrs20lags_L_Llead_synch    = squeeze(mean(mean(squeeze(Autocorrs20lags(1,[4:6],[1:2:11],:)),2),1));   % L,L lead,synch
Autocorrs20lags_L_Rlead_synch    = squeeze(mean(mean(squeeze(Autocorrs20lags(1,[7:9],[1:2:11],:)),2),1));   % L,R lead,synch
Autocorrs20lags_L_mutual_synch   = squeeze(mean(mean(squeeze(Autocorrs20lags(1,[10:12],[1:2:11],:)),2),1)); % L, mutual,synch
Autocorrs20lags_R_uncouple_synch = squeeze(mean(mean(squeeze(Autocorrs20lags(2,[1:3],[1:2:11],:)),2),1));   % R,uncouple,synch
Autocorrs20lags_R_Llead_synch    = squeeze(mean(mean(squeeze(Autocorrs20lags(2,[4:6],[1:2:11],:)),2),1));   % R,L lead,synch
Autocorrs20lags_R_Rlead_synch    = squeeze(mean(mean(squeeze(Autocorrs20lags(2,[7:9],[1:2:11],:)),2),1));   % R,R lead,synch
Autocorrs20lags_R_mutual_synch   = squeeze(mean(mean(squeeze(Autocorrs20lags(2,[10:12],[1:2:11],:)),2),1)); % R, mutual,synch
Autocorrs20lags_LR_synch_mean = [Autocorrs20lags_L_uncouple_synch  Autocorrs20lags_R_uncouple_synch ...
    Autocorrs20lags_L_Llead_synch  Autocorrs20lags_R_Llead_synch ...
    Autocorrs20lags_L_Rlead_synch  Autocorrs20lags_R_Rlead_synch ...
    Autocorrs20lags_L_mutual_synch  Autocorrs20lags_R_mutual_synch];

Autocorrs20lags_L_uncouple_synco = squeeze(mean(mean(squeeze(Autocorrs20lags(1,[1:3],[2:2:12],:)),2),1));   % R,uncouple,synco
Autocorrs20lags_L_Llead_synco    = squeeze(mean(mean(squeeze(Autocorrs20lags(1,[4:6],[2:2:12],:)),2),1));   % R,L lead,synco
Autocorrs20lags_L_Rlead_synco    = squeeze(mean(mean(squeeze(Autocorrs20lags(1,[7:9],[2:2:12],:)),2),1));   % R,R lead,synco
Autocorrs20lags_L_mutual_synco   = squeeze(mean(mean(squeeze(Autocorrs20lags(1,[10:12],[2:2:12],:)),2),1)); % R, mutual,synco
Autocorrs20lags_R_uncouple_synco = squeeze(mean(mean(squeeze(Autocorrs20lags(2,[1:3],[2:2:12],:)),2),1));   % R,uncouple,synco
Autocorrs20lags_R_Llead_synco    = squeeze(mean(mean(squeeze(Autocorrs20lags(2,[4:6],[2:2:12],:)),2),1));   % R,L lead,synco
Autocorrs20lags_R_Rlead_synco    = squeeze(mean(mean(squeeze(Autocorrs20lags(2,[7:9],[2:2:12],:)),2),1));  % R,R lead,synco
Autocorrs20lags_R_mutual_synco   = squeeze(mean(mean(squeeze(Autocorrs20lags(2,[10:12],[2:2:12],:)),2),1)); % R, mutual,synco
Autocorrs20lags_LR_synco_mean = [Autocorrs20lags_L_uncouple_synco Autocorrs20lags_R_uncouple_synco ...
    Autocorrs20lags_L_Llead_synco Autocorrs20lags_R_Llead_synco ...
    Autocorrs20lags_L_Rlead_synco Autocorrs20lags_R_Rlead_synco ...
    Autocorrs20lags_L_mutual_synco Autocorrs20lags_R_mutual_synco];

canvas(0.4,0.5);
condi4names={ 'Uncouple','L-lead','R-lead','Mutual'};
for i=1:4
    subplot(2,4,i);
    plot(Autocorrs20lags_LR_synch_mean(:,2*i-1),'r');hold on;
    plot(Autocorrs20lags_LR_synch_mean(:,2*i),'b');legend('L','R');
    title(condi4names{i},'color',condicolors(i,:));
    subtitle('autocorr synch');
    ylabel('\rho(k)');xlabel('lag');ylim([-0.1 0.5])
    subplot(2,4,4+i);
    plot(Autocorrs20lags_LR_synco_mean(:,2*i-1),'r');hold on;
    plot(Autocorrs20lags_LR_synco_mean(:,2*i),'b');legend('L','R');
    title(condi4names{i},'color',condicolors(i,:));
    subtitle('autocorr synco');
    ylabel('\rho(k)');xlabel('lag'); ylim([-0.1 0.5])
end
% suptitle('mean autocorr for 20 lags from all sessions (unmatched int) ^{ *PLOT7}')
suptitle('mean autocorr for 20 lags from all sessions (matched int) ^{ *PLOT7}')

% combine synch and synco of Autocorrs (sum) in bar graph
Autocorrs_LR_syn_mean=[mean([Autocorrs_L_uncouple_synch Autocorrs_L_uncouple_synco Autocorrs_R_uncouple_synch Autocorrs_R_uncouple_synco]);...
    mean([Autocorrs_L_Llead_synch Autocorrs_L_Llead_synco Autocorrs_R_Rlead_synch Autocorrs_R_Rlead_synco]);... 
    mean([Autocorrs_R_Llead_synch Autocorrs_R_Llead_synco Autocorrs_L_Rlead_synch Autocorrs_L_Rlead_synco]);...
    mean([Autocorrs_L_mutual_synch Autocorrs_L_mutual_synco Autocorrs_R_mutual_synch Autocorrs_R_mutual_synco])];
figure;
b=bar(Autocorrs_LR_syn_mean);
xticks(1:4);xticklabels({'uncouple','Leading','Following','mutual'});
ylabel('\rho(k)');
hold on;
% title('mean autocorr sum for 20 lags in all sessions (unmatched int) ^{ *PLOT7}');
title('mean autocorr sum for 20 lags in all sessions (matched int) ^{ *PLOT7}');

% combine synch and synco of Autocorrs (sum) in violin plot
addpath /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/violin
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
xticks([1:4]);
set(gca,'XtickLabel',{'uncouple','Leading','Following','mutual'});
ylabel('\rho(k)');
% title('mean autocorr sum for 20 lags in all sessions (unmatched int) ^{ *PLOT 7}');
title('mean autocorr sum for 20 lags in all sessions (matched int) ^{ *PLOT 7}');
grid on;

% combine synch and synco of Autocorrs20lags in 3 directions
direction3names={'Independent','Unidirectional','Bidirectional'};
Autocorrs20lags_LR_synch_3mean = [Autocorrs20lags_L_uncouple_synch  Autocorrs20lags_R_uncouple_synch ...
    mean([Autocorrs20lags_L_Llead_synch Autocorrs20lags_R_Rlead_synch],2) ...
    mean([Autocorrs20lags_R_Llead_synch Autocorrs20lags_L_Rlead_synch],2) ...
    Autocorrs20lags_L_mutual_synch  Autocorrs20lags_R_mutual_synch];
Autocorrs20lags_LR_synco_3mean = [Autocorrs20lags_L_uncouple_synco Autocorrs20lags_R_uncouple_synco ...
    mean([Autocorrs20lags_L_Llead_synco Autocorrs20lags_R_Rlead_synco],2) ...
    mean([Autocorrs20lags_R_Llead_synco Autocorrs20lags_L_Rlead_synco],2) ...
    Autocorrs20lags_L_mutual_synco Autocorrs20lags_R_mutual_synco];
canvas(0.4, 0.5);
for i=1:3
subplot(2,3,i)
    plot(Autocorrs20lags_LR_synch_3mean(:,2*i-1),'r');hold on;
    plot(Autocorrs20lags_LR_synch_3mean(:,2*i),'b');legend('L','R');
    if i==2; legend('Leader','follower');end
    title(direction3names{i},'color',condicolors(i,:));
    subtitle('synch autocorr');
    ylabel('\rho(k)');xlabel('lag');ylim([0 0.5])
subplot(2,3,3+i)
    plot(Autocorrs20lags_LR_synco_3mean(:,2*i-1),'r');hold on;
    plot(Autocorrs20lags_LR_synco_3mean(:,2*i),'b');legend('L','R');
    if i==2; legend('Leader','follower');end
    title(direction3names{i},'color',condicolors(i,:));
    subtitle('synco autocorr');
    ylabel('\rho(k)');xlabel('lag');ylim([0 0.5])
end
% sgtitle('mean autocorr sum for 20 lags in all sessions (unmatched int) ^{ *PLOT7}');
sgtitle('mean autocorr sum for 20 lags in all sessions (matched int) ^{ *PLOT7}');

Autocorrs20lags_LR_uncouple_mean=...
    [mean([Autocorrs20lags_L_uncouple_synch Autocorrs20lags_L_uncouple_synco],2) ...
    mean([Autocorrs20lags_R_uncouple_synch Autocorrs20lags_R_uncouple_synco],2)];
Autocorrs20lags_LR_unidirection_mean=...
    [mean([Autocorrs20lags_L_Llead_synch Autocorrs20lags_L_Llead_synco ...
    Autocorrs20lags_R_Rlead_synch Autocorrs20lags_R_Rlead_synco],2) ...
    mean([Autocorrs20lags_R_Llead_synch Autocorrs20lags_R_Llead_synco ...
    Autocorrs20lags_L_Rlead_synch Autocorrs20lags_L_Rlead_synco],2)];
Autocorrs20lags_LR_mutual_mean=...
    [mean([Autocorrs20lags_L_mutual_synch Autocorrs20lags_L_mutual_synco],2) ...
    mean([Autocorrs20lags_R_mutual_synch Autocorrs20lags_R_mutual_synco],2)];
Autocorrs20lags_3dir_LR_syn_mean=...
    [Autocorrs20lags_LR_uncouple_mean Autocorrs20lags_LR_unidirection_mean Autocorrs20lags_LR_mutual_mean];
canvas(0.2, 0.25);
for i=1:3
subplot(1,3,i)
    plot(Autocorrs20lags_3dir_LR_syn_mean(:,2*i-1),'r');hold on;
    plot(Autocorrs20lags_3dir_LR_syn_mean(:,2*i),'b');legend('L','R');
    if i==2; legend('Leader','follower');end
    title(direction3names{i});
    subtitle('autocorr');
    ylabel('\rho(k)');xlabel('lag');ylim([0 0.5])
end
% sgtitle('mean autocorr sum for 20 lags in all sessions (unmatched int) ^{ *PLOT7}');
% sgtitle('mean autocorr sum for 20 lags in all sessions (matched int) ^{ *PLOT7}');
sg=annotation('textbox',[0.3 0.01 0.4 0.05],'string',...
    'Autocorr ^{* PLOT 7}')
set(gcf,'color','w'); % set background white for copying in ubuntu

%% PLOT 8 DFA after d removed (sorted order)
clear
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);
sessions={'synch','synco','synch','synco','synch','synco','synch','synco','synch','synco',...
    'synch','synco'};

H=zeros(2,12,numSes);
cd /ssd/zhibin/1overf/
for r=1:numSes;
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
H_est1=zeros(2,12,numSes);
H_est2=zeros(2,12,numSes);
int_dmean_drm1=cell(2,12,numSes);
cd /ssd/zhibin/1overf/
for r=1:numSes;
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

% % to plot H_est1
% H_est2=H_est1;

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
H_est2_L_synch5=H_est2(1,:,9)';H_est2_R_synch5=H_est2(2,:,9)';
H_est2_L_synco5=H_est2(1,:,10)';H_est2_R_synco5=H_est2(2,:,10)';
H_est2_L_synch6=H_est2(1,:,11)';H_est2_R_synch6=H_est2(2,:,11)';
H_est2_L_synco6=H_est2(1,:,12)';H_est2_R_synco6=H_est2(2,:,12)';
Table_H_est2=table(Trials, H_est2_L_synch1, H_est2_R_synch1, H_est2_L_synco1, H_est2_R_synco1,...
    H_est2_L_synch2, H_est2_R_synch2, H_est2_L_synco2, H_est2_R_synco2,...
    H_est2_L_synch3, H_est2_R_synch3, H_est2_L_synco3, H_est2_R_synco3,...
    H_est2_L_synch4, H_est2_R_synch4, H_est2_L_synco4, H_est2_R_synco4);

% plot H_est2 after
canvas(0.5, 0.5);
for r=1:numSes
    subplot(numSes,1,r);
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
sgtitle('H from DFA after d removal in all sessions (using Leites d estimation) ^{* PLOT 8}');
sgtitle('H from DFA after d removal in all sessions (using d=H-0.5) ^{* PLOT 8}');


%% SECT 9 Xcorr before/after d removal and save new intervals (matched int) (sorted order) 
Trials={'uncouple','uncouple','uncouple','L-lead','L-lead','L-lead',...
    'R-lead','R-lead','R-lead','mutual','mutual','mutual'}';
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);
sessions={'synch','synco','synch','synco','synch','synco','synch','synco','synch','synco',...
    'synch','synco'};

% (Before d removal)********
% d estimation from H (DFA method) based on the good matched intervals
H=zeros(2,12,numSes);XcorrPeakLag=zeros(12,numSes);XcorrPeak=zeros(12,numSes);
Xcorr10Lag=nan(numSes,12,21);
cd /ssd/zhibin/1overf/
tic
for r=1:numSes
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
        % Xcorr based on int_dmean_drm (before d removal)********
        r12=[];lags12=[];
        [r12,lags12]=xcorr(intL_good_dmean,intR_good_dmean,10,'normalized');
        XcorrPeakLag(j,r)=lags12(find(r12==max(r12)));
        XcorrPeak(j,r)=max(r12);
        Xcorr10Lag(r,j,1:21)=r12;
    end
end
toc
% 80 sec
 
% eval(['save ' '/ssd/zhibin/1overf/all_session20220713_0816/H.mat H'])
eval(['save ' '/ssd/zhibin/1overf/all_session20220713_1005/H.mat H']);
d1=H-0.5;% estimate d using DFA method (method1: based on DFA, d=h-0.5)
% load ('/ssd/zhibin/1overf/all_session20220713_1005/H.mat H');

% ( Continue to PLOT 9 )#########################
        
% d removal
H_est2=zeros(2,12,numSes);
int_dmean_drm=cell(2,12,numSes);
cd /ssd/zhibin/1overf/
tic
for r=1:numSes;
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
toc

% save the "int_dmean_drm" under [/ssd/zhibin/1overf/all_session20220713_0816]
% eval(['save ' '/ssd/zhibin/1overf/all_session20220713_0816/int_dmean_drm.mat int_dmean_drm'])
eval(['save ' '/ssd/zhibin/1overf/all_session20220713_1005/int_dmean_drm.mat int_dmean_drm'])
% plot H all session and compare like in PLOT 8

% (After d removal)********
% Xcorr based on int_dmean_drm
% load /ssd/zhibin/1overf/all_session20220713_0816/int_dmean_drm.mat
% load /ssd/zhibin/1overf/all_session20220713_0816/H.mat
load /ssd/zhibin/1overf/all_session20220713_1005/int_dmean_drm.mat
load /ssd/zhibin/1overf/all_session20220713_1005/H.mat
XcorrPeakLag=zeros(12,numSes);XcorrPeak=zeros(12,numSes);
for r=1:numSes;
    for j = 1:12
        % Xcorr after d removal ********
        r12=[];lags12=[];
        [r12,lags12]=xcorr(int_dmean_drm{1,j,r},int_dmean_drm{2,j,r},10,'normalized');
        XcorrPeakLag(j,r)=lags12(find(r12==max(r12)));
        XcorrPeak(j,r)=max(r12);
    end
end
% ( Continue to PLOT 9 again) #########################
%% PLOT 9 Xcorr before/after d removal (matched int) (sorted order) + ANOVA 
% Plot XcorrPeakLag in all sessions
canvas(0.3, 0.8);
for r=1:12
    subplot(numSes,1,r);
    plot(1:12,XcorrPeakLag(:,r),'ko')
    ylabel('lag');ylim([-10 10]);
    yline(0,'color',deepyellow);set(gca, 'YDir','reverse')
    xticks(1:12);xticklabels(Trials);
    title([num2str(seeds(r,:)) '-' sessions{r}]);
    grid on; grid minor;
end
sgtitle({['peak lag from Xcorr in all sessions (-lag = L leading; +lag = R leading)'],...
    ['(matched int) (before d removal) ^{* PLOT 9}']});
sgtitle({['peak lag from Xcorr in all sessions (-lag = L leading; +lag = R leading)'],...
    ['(matched int) (after d removal) ^{* PLOT 9}']});


% print XcorrPeakLag to a table
for table1=1;
XcorrPeakLag_synch1=XcorrPeakLag(:,1);
XcorrPeakLag_synco1=XcorrPeakLag(:,2);
XcorrPeakLag_synch2=XcorrPeakLag(:,3);
XcorrPeakLag_synco2=XcorrPeakLag(:,4);
XcorrPeakLag_synch3=XcorrPeakLag(:,5);
XcorrPeakLag_synco3=XcorrPeakLag(:,6);
XcorrPeakLag_synch4=XcorrPeakLag(:,7);
XcorrPeakLag_synco4=XcorrPeakLag(:,8);
XcorrPeakLag_synch5=XcorrPeakLag(:,9);
XcorrPeakLag_synco5=XcorrPeakLag(:,10);
XcorrPeakLag_synch6=XcorrPeakLag(:,11);
XcorrPeakLag_synco6=XcorrPeakLag(:,12);
end
Table_XcorrPeakLag=table(Trials, XcorrPeakLag_synch1,XcorrPeakLag_synco1,...
    XcorrPeakLag_synch2, XcorrPeakLag_synco2,...
    XcorrPeakLag_synch3, XcorrPeakLag_synco3,...
    XcorrPeakLag_synch4, XcorrPeakLag_synco4,...
    XcorrPeakLag_synch5, XcorrPeakLag_synco5,...
    XcorrPeakLag_synch6, XcorrPeakLag_synco6);
Table_XcorrPeakLag


% organize for PLOT (lag of peak)
% average all trials of the same condition 
% synchronization
XcorrPeakLag_uncouple_synch = reshape(XcorrPeakLag([1:3],[1:2:11]),1,[])'  % uncouple,synch
XcorrPeakLag_Llead_synch    = reshape(XcorrPeakLag([4:6],[1:2:11]),1,[])'  % L lead,synch
XcorrPeakLag_Rlead_synch    = reshape(XcorrPeakLag([7:9],[1:2:11]),1,[])'  % R lead,synch
XcorrPeakLag_mutual_synch   = reshape(XcorrPeakLag([10:12],[1:2:11]),1,[])'% mutual,synch
% print on a table 
Table_XcorrPeakLagSynch_4condi=table(...
    XcorrPeakLag_uncouple_synch,XcorrPeakLag_Llead_synch,...
    XcorrPeakLag_Rlead_synch,XcorrPeakLag_mutual_synch);
Table_XcorrPeakLagSynch_4condi
% ANOVA
% for MATLAB (https://www.mathworks.com/help/stats/one-way-anova.html)
% https://www.mathworks.com/help/stats/multiple-comparisons.html
% (https://www.mathworks.com/matlabcentral/fileexchange/69794-bayesfactor?s_tid=mwa_osa_a)
Table_XcorrPeakLagSynch_4condi=[
    XcorrPeakLag_uncouple_synch XcorrPeakLag_Llead_synch ...
    XcorrPeakLag_Rlead_synch XcorrPeakLag_mutual_synch];
[p,tbl,stats] = anova1(Table_XcorrPeakLagSynch_4condi); % one-way ANOVA (work)
p_synch=p;
[results,means,~,gnames] = multcompare(stats,"CriticalValueType","bonferroni"); % multi-comparison
tbl = array2table([results],"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])
tb2 = array2table([means],"VariableNames", ["Mean","Standard Error"])
p_synch_bonf=results(:,6);
% print out on command window
sprintf(['one-way ANOVA: p = %0.3f\n' ...
    'Bonferroni correction: \n'...
    'control vs L-lead: p = %0.3f\n'...
    'control vs R-lead: p = %0.3f\n'...
    'control vs Mutual: p = %0.3f\n'...
    ' L-lead vs R-lead: p = %0.3f\n'...
    ' L-lead vs Mutual: p = %0.3f\n'...
    ' R-lead vs Mutual: p = %0.3f']...
    , p_synch, p_synch_bonf)
% p_synch_bonf(1), p_synch_bonf(2), p_synch_bonf(3),p_synch_bonf(4), p_synch_bonf(5), p_synch_bonf(6)
% for JASP
Table_XcorrPeakLagSynch_4condi=table(...
    [XcorrPeakLag_uncouple_synch;XcorrPeakLag_Llead_synch;...
    XcorrPeakLag_Rlead_synch;XcorrPeakLag_mutual_synch],...
    [ones(12,1);2*ones(12,1);3*ones(12,1);4*ones(12,1)],'VariableNames',{'XcorrLag','condition'});
% save for JASP
cd /usr/local/MATLAB/R2022a/toolbox/matlab/strfun/
writetable(Table_XcorrPeakLagSynch_4condi,'/ssd/zhibin/1overf/all_session20220713_0816/Table_XcorrPeakLagSynch_4condi.csv');
cd /ssd/zhibin/1overf/all_session20220713_0816
% get mean and std for plot
XcorrPeakLag_LR_synch_mean = [mean(XcorrPeakLag_uncouple_synch);...
    mean(XcorrPeakLag_Llead_synch);...
    mean(XcorrPeakLag_Rlead_synch);...
    mean(XcorrPeakLag_mutual_synch)];
XcorrPeakLag_LR_synch_std = [std(XcorrPeakLag_uncouple_synch);...
    std(XcorrPeakLag_Llead_synch);...
    std(XcorrPeakLag_Rlead_synch);...
    std(XcorrPeakLag_mutual_synch)];

% syncopation
XcorrPeakLag_uncouple_synco = reshape(XcorrPeakLag([1:3],[2:2:12]),1,[])';   % R,uncouple,synco
XcorrPeakLag_Llead_synco    = reshape(XcorrPeakLag([4:6],[2:2:12]),1,[])';   % R,L lead,synco
XcorrPeakLag_Rlead_synco    = reshape(XcorrPeakLag([7:9],[2:2:12]),1,[])';   % R,R lead,synco
XcorrPeakLag_mutual_synco   = reshape(XcorrPeakLag([10:12],[2:2:12]),1,[])'; % R, mutual,synco
% print on a table 
Table_XcorrPeakLagSynco_4condi=table(...
    XcorrPeakLag_uncouple_synco,XcorrPeakLag_Llead_synco,...
    XcorrPeakLag_Rlead_synco,XcorrPeakLag_mutual_synco);
Table_XcorrPeakLagSynco_4condi
% ANOVA
% for MATLAB 
Table_XcorrPeakLagSynco_4condi=[
    XcorrPeakLag_uncouple_synco,XcorrPeakLag_Llead_synco,...
    XcorrPeakLag_Rlead_synco,XcorrPeakLag_mutual_synco];
p=[];tbl=[];stats=[];
[p,tbl,stats] = anova1(Table_XcorrPeakLagSynco_4condi); % one-way ANOVA (work)
p_synco=p;
[results,means,~,gnames] = multcompare(stats,"CriticalValueType","bonferroni"); % multi-comparison
tbl = array2table([results],"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])
tb2 = array2table([means],"VariableNames", ["Mean","Standard Error"])
p_synco_bonf=results(:,6);
% print out on command window
sprintf(['one-way ANOVA: p = %0.3f\n' ...
    'Bonferroni correction: \n'...
    'control vs L-lead: p = %0.3f\n'...
    'control vs R-lead: p = %0.3f\n'...
    'control vs Mutual: p = %0.3f\n'...
    ' L-lead vs R-lead: p = %0.3f\n'...
    ' L-lead vs Mutual: p = %0.3f\n'...
    ' R-lead vs Mutual: p = %0.3f']...
    , p_synco, p_synco_bonf)
% for JASP
Table_XcorrPeakLagSynco_4condi=table(...
    [XcorrPeakLag_uncouple_synco;XcorrPeakLag_Llead_synco;...
    XcorrPeakLag_Rlead_synco;XcorrPeakLag_mutual_synco],...
    [ones(12,1);2*ones(12,1);3*ones(12,1);4*ones(12,1)],'VariableNames',{'XcorrLag','condition'});
% save for JASP
cd /usr/local/MATLAB/R2022a/toolbox/matlab/strfun/
writetable(Table_XcorrPeakLagSynco_4condi,'/ssd/zhibin/1overf/all_session20220713_0816/Table_XcorrPeakLagSynco_4condi.csv');
cd /ssd/zhibin/1overf/all_session20220713_0816
% get mean and std for plot
XcorrPeakLag_LR_synco_mean = [mean(XcorrPeakLag_uncouple_synco);...
    mean(XcorrPeakLag_Llead_synco);...
    mean(XcorrPeakLag_Rlead_synco);...
    mean(XcorrPeakLag_mutual_synco)];
XcorrPeakLag_LR_synco_std = [std(XcorrPeakLag_uncouple_synco);...
    std(XcorrPeakLag_Llead_synco);...
    std(XcorrPeakLag_Rlead_synco);...
    std(XcorrPeakLag_mutual_synco)];
% plot average lag of peak
canvas(0.2, 0.6);
subplot(2,1,1);
errorbar(1:4,XcorrPeakLag_LR_synch_mean,XcorrPeakLag_LR_synch_std,'ko');
xticks(1:4);xticklabels({ 'Uncouple','L-lead','R-lead','Mutual'});
ylabel('average lag of the peak'); % ylabel('\rho(k)');
ylim([-10 10]);yline(0,'color',red);set(gca, 'YDir','reverse');
xlim([0 5]);grid on;  grid minor;
title('average lag of xcorr peak in all synch sessions (-lag: L leading; +lag: R leading) ^{* PLOT 9}');
% subtitle('(matched int) (before d removal)');
subtitle('(matched int) (after d removal)');
subplot(2,1,2);
errorbar(1:4,XcorrPeakLag_LR_synco_mean,XcorrPeakLag_LR_synco_std,'ko');
xticks(1:4);xticklabels({ 'Uncouple','L-lead','R-lead','Mutual'});
ylabel('average lag of the peak'); % ylabel('\rho(k)');
ylim([-10 10]);yline(0,'color',red);set(gca, 'YDir','reverse')
xlim([0 5]);grid on; grid minor;
title('average lag of xcorr peak in all synco sessions (-lag: L leading; +lag: R leading) ^{* PLOT 9}');
% subtitle('(matched int) (before d removal)');
subtitle('(matched int) (after d removal)');

% organize for PLOT (peak)
% for synchronization
% average all trials of the same condition 
XcorrPeak_uncouple_synch = reshape(XcorrPeak([1:3],[1:2:11]),1,[])'  % uncouple,synch
XcorrPeak_Llead_synch    = reshape(XcorrPeak([4:6],[1:2:11]),1,[])'  % L lead,synch
XcorrPeak_Rlead_synch    = reshape(XcorrPeak([7:9],[1:2:11]),1,[])'  % R lead,synch
XcorrPeak_mutual_synch   = reshape(XcorrPeak([10:12],[1:2:11]),1,[])'% mutual,synch
% print on a table
Table_XcorrPeakSynch_4condi=table(...
    XcorrPeak_uncouple_synch,XcorrPeak_Llead_synch,...
    XcorrPeak_Rlead_synch,XcorrPeak_mutual_synch);
Table_XcorrPeakSynch_4condi
% for MATLAB 
Table_XcorrPeakSynch_4condi=[
    XcorrPeak_uncouple_synch XcorrPeak_Llead_synch ...
    XcorrPeak_Rlead_synch XcorrPeak_mutual_synch];
clear p tbl stats
[p,tbl,stats] = anova1(Table_XcorrPeakSynch_4condi); % one-way ANOVA (work)
p_synch=p;
[results,means,~,gnames] = multcompare(stats,"CriticalValueType","bonferroni"); % multi-comparison
tbl = array2table([results],"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])
tb2 = array2table([means],"VariableNames", ["Mean","Standard Error"])
p_synch_bonf=results(:,6);
% print out on command window
sprintf(['one-way ANOVA: p = %0.3f\n' ...
    'Bonferroni correction: \n'...
    'control vs L-lead: p = %0.3f\n'...
    'control vs R-lead: p = %0.3f\n'...
    'control vs Mutual: p = %0.3f\n'...
    ' L-lead vs R-lead: p = %0.3f\n'...
    ' L-lead vs Mutual: p = %0.3f\n'...
    ' R-lead vs Mutual: p = %0.3f']...
    , p_synch, p_synch_bonf)
% for jasp
Table_XcorrPeakSynch_4condi=table(...
    [XcorrPeak_uncouple_synch;XcorrPeak_Llead_synch;...
    XcorrPeak_Rlead_synch;XcorrPeak_mutual_synch],...
    [ones(12,1);2*ones(12,1);3*ones(12,1);4*ones(12,1)],'VariableNames',{'XcorrPeak','condition'});
% save for JASP
cd /usr/local/MATLAB/R2022a/toolbox/matlab/strfun/
writetable(Table_XcorrPeakSynch_4condi,'/ssd/zhibin/1overf/all_session20220713_0816/Table_XcorrPeakSynch_4condi.csv');
cd /ssd/zhibin/1overf/all_session20220713_0816
% get mean and std for plot
XcorrPeak_LR_synch_mean = [mean(XcorrPeak_uncouple_synch);...
    mean(XcorrPeak_Llead_synch);...
    mean(XcorrPeak_Rlead_synch);...
    mean(XcorrPeak_mutual_synch)];
XcorrPeak_LR_synch_std = [std(XcorrPeak_uncouple_synch);...
    std(XcorrPeak_Llead_synch);...
    std(XcorrPeak_Rlead_synch);...
    std(XcorrPeak_mutual_synch)];
% for syncopation
XcorrPeak_uncouple_synco = reshape(XcorrPeak([1:3],[2:2:12]),1,[])'  % R,uncouple,synco
XcorrPeak_Llead_synco    = reshape(XcorrPeak([4:6],[2:2:12]),1,[])'  % R,L lead,synco
XcorrPeak_Rlead_synco    = reshape(XcorrPeak([7:9],[2:2:12]),1,[])'  % R,R lead,synco
XcorrPeak_mutual_synco   = reshape(XcorrPeak([10:12],[2:2:12]),1,[])'% R, mutual,synco
% print on a table 
Table_XcorrPeakSynco_4condi=table(...
    XcorrPeak_uncouple_synco,XcorrPeak_Llead_synco,...
    XcorrPeak_Rlead_synco,XcorrPeak_mutual_synco);
Table_XcorrPeakSynco_4condi
% for MATLAB 
Table_XcorrPeakSynco_4condi=[
    XcorrPeak_uncouple_synco XcorrPeak_Llead_synco ...
    XcorrPeak_Rlead_synco XcorrPeak_mutual_synco];
clear p tbl stats
[p,tbl,stats] = anova1(Table_XcorrPeakSynco_4condi); % one-way ANOVA (work)
p_synch=p;
[results,means,~,gnames] = multcompare(stats,"CriticalValueType","bonferroni"); % multi-comparison
tbl = array2table([results],"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])
tb2 = array2table([means],"VariableNames", ["Mean","Standard Error"])
p_synch_bonf=results(:,6);
% print out on command window
sprintf(['one-way ANOVA: p = %0.3f\n' ...
    'Bonferroni correction: \n'...
    'control vs L-lead: p = %0.3f\n'...
    'control vs R-lead: p = %0.3f\n'...
    'control vs Mutual: p = %0.3f\n'...
    ' L-lead vs R-lead: p = %0.3f\n'...
    ' L-lead vs Mutual: p = %0.3f\n'...
    ' R-lead vs Mutual: p = %0.3f']...
    , p_synch, p_synch_bonf)
% for jasp
Table_XcorrPeakSynco_4condi=table(...
    [XcorrPeak_uncouple_synco;XcorrPeak_Llead_synco;...
    XcorrPeak_Rlead_synco;XcorrPeak_mutual_synco],...
    [ones(12,1);2*ones(12,1);3*ones(12,1);4*ones(12,1)],'VariableNames',{'XcorrPeak','condition'});
% save for JASP
cd /usr/local/MATLAB/R2022a/toolbox/matlab/strfun/
writetable(Table_XcorrPeakSynco_4condi,'/ssd/zhibin/1overf/all_session20220713_0816/Table_XcorrPeakSynco_4condi.csv');
cd /ssd/zhibin/1overf/all_session20220713_0816
% get mean and std for plot
XcorrPeak_LR_synco_mean = [mean(XcorrPeak_uncouple_synco);...
    mean(XcorrPeak_Llead_synco);...
    mean(XcorrPeak_Rlead_synco);...
    mean(XcorrPeak_mutual_synco)];
XcorrPeak_LR_synco_std = [std(XcorrPeak_uncouple_synco);...
    std(XcorrPeak_Llead_synco);...
    std(XcorrPeak_Rlead_synco);...
    std(XcorrPeak_mutual_synco)];
% plot average peak
% condi4names={ 'Uncouple','L-lead','R-lead','Mutual'};
canvas(0.1, 0.5);
subplot(2,1,1);
errorbar(1:4,XcorrPeak_LR_synch_mean,XcorrPeak_LR_synch_std,'ko');
% for i=1:4
%     errorbar(i,XcorrPeak_LR_synch_mean(i),XcorrPeak_LR_synch_std(i),'o','Color',condicolors(i,:));
%     hold on;
% end
xticks(1:4);xticklabels({ 'Uncouple','L-lead','R-lead','Mutual'});
ylabel('\rho(k)'); grid on;  grid minor;
ylim([-0.1 1]);yline(0,'color',red);xlim([0 5]);
title('average peak of xcorr in all synch sessions ^{* PLOT 9}');
% subtitle('(matched int) (before d removal)');
subtitle('(matched int) (after d removal)');
subplot(2,1,2);
errorbar(1:4,XcorrPeak_LR_synco_mean,XcorrPeak_LR_synco_std,'ko');
xticks(1:4);xticklabels({ 'Uncouple','L-lead','R-lead','Mutual'});
ylabel('\rho(k)'); grid on;  grid minor;
ylim([-0.1 1]);yline(0,'color',red);xlim([0 5]);
title('average peak of xcorr in all synco sessions ^{* PLOT 9}');
% subtitle('(matched int) (before d removal)');
subtitle('(matched int) (after d removal)');
%

%% PLOT 9-1 Xcorr10Lag before d removal (matched int) (sorted order)
% average all trials of the same condition (synch vs synco)
Xcorr10Lag; % sorted order (12 ses x 12 blocks x 21 lags) from SECT 9
sorted4inds=[1:3; 4:6; 7:9; 10:12];
% Examine all Xcorr10Lag
canvas(1,1);
for s=1:12
    for b=1:12
        subplot(12,12,12*(s-1)+b)
        plot(-10:1:10,squeeze(Xcorr10Lag(s,b,:)),'k');
        i=[];
        if any([1:3]==b);
            i=1;
        elseif any([4:6]==b);
            i=2;
        elseif any([7:9]==b);
            i=3;
        else any([10:12]==b);
            i=4;
        end;
        title(condi4names{i},'color',condicolors(i,:));
        ylabel('\rho(k)');xlabel('lag');ylim([-0.1 0.8])
        yline(0,'color',[1 0.8 0.2]);xline(0,'color',[1 0.8 0.2]);
        xline(-1,'color',[1 0.8 0.2]);xline(1,'color',[1 0.8 0.2]);
    end
end

% Examine Xcorr10Lag in Mutual Synch, to see how many trials have two peak
figure
for s=[1 3 5 7 9 11]
    for b=10:12
        subplot(6,3,3*((s+1)/2-1)+b-9)
        plot(-10:1:10,squeeze(Xcorr10Lag(s,b,:)),'k');
        title([syn2names{1} '   ' condi4names{4}],'color',condicolors(i,:));
        ylabel('\rho(k)');xlabel('lag');ylim([-0.2 1])
        yline(0,'color',[1 0.8 0.2]);xline(0,'color',[1 0.8 0.2]);
        xline(-1,'color',[1 0.8 0.2]);xline(1,'color',[1 0.8 0.2]);
    end
end

% organized into 2x4 conditions (now 1x4 with diff colors)
canvas(0.3,0.2);
condi4names={ 'Uncouple','L-lead','R-lead','Mutual'};
for i=1:4
    subplot(1,4,i);
    plot(-10:1:10,mean(squeeze(mean(Xcorr10Lag(1:2:11,sorted4inds(i,:),:),2)),1),'b');hold on;
    plot(-10:1:10,mean(squeeze(mean(Xcorr10Lag(2:2:12,sorted4inds(i,:),:),2)),1),'r');
    title(condi4names{i},'color',condicolors(i,:));
    subtitle('xcorr synco');
    ylabel('\rho(k)');xlabel('lag');ylim([-0.1 0.8])
    yline(0,'color',[1 0.8 0.2]);xline(0,'color',[1 0.8 0.2]);
end
lg=legend('Synch','Synco','location','eastoutside');
lg.Position = [0.9475 0.25 0.01 0.25];
sgtitle('mean xcorr for 10 lags from all sessions (matched int) ^{ *PLOT 9-1}')
set(gcf,'color','w'); % set background white for copying in ubuntu

% combine Xcorr10Lag in 3 directions (2x3suplots)
Xcorr10Lag; % sorted order
direction3names={'Independent','Unidirectional','Bidirectional'};
sorted3inds={[1:3],[4:9],[10:12]};
canvas(0.45, 0.6);
for i=1:3
subplot(2,3,i)
    plot(-10:1:10,mean(squeeze(mean(Xcorr10Lag(1:2:11,sorted3inds{i},:),2)),1),'k');
    title(direction3names{i},'color',condicolors(i,:));
    subtitle('xcorr synch');ylabel('\rho');ylim([-0.1 0.8])
    xlabel('lag (L-leading <- 0 -> R-leading)');
    yline(0,'color',[1 0.8 0.2]);xline(0,'color',[1 0.8 0.2]);
    subplot(2,3,3+i)
    plot(-10:1:10,mean(squeeze(mean(Xcorr10Lag(2:2:12,sorted3inds{i},:),2)),1),'k');
    title(direction3names{i},'color',condicolors(i,:));
    subtitle('xcorr synco');ylabel('\rho');ylim([-0.1 0.8])
    xlabel('lag (L-leading <- 0 -> R-leading)');
    yline(0,'color',[1 0.8 0.2]);xline(0,'color',[1 0.8 0.2]);
end
sgtitle('mean xcorr sum for 10 lags in all sessions (matched int) ^{ *PLOT 9-1}');
% in 1x3 with different colors and 

% combine synch and synco in 1x3 
Xcorr10Lag; % sorted order
direction3names={'Independent','Unidirectional','Bidirectional'};
sorted3inds={[1:3],[4:9],[10:12]};
sorted4inds=[1:3; 4:6; 7:9; 10:12];
canvas(0.45, 0.35);
for i=1:3
subplot(1,3,i)
    if i==2;
        % L-lead
        L_mat=squeeze(mean(Xcorr10Lag(:,sorted4inds(2,:),:),2));
        % R-lead
        R_mat=squeeze(mean(Xcorr10Lag(:,sorted4inds(3,:),:),2));
        % combine
        LR_mean=mean([L_mat;fliplr(R_mat)]);
        % plot Leader vs Follower
        plot(-10:1:10,LR_mean,'k');
        title(direction3names{i});
        ylabel('\rho(k)');xlabel({'lag'; '(Leader-leading <- 0 -> Follower-leading)'});ylim([-0.1 0.8])
    elseif i==1 | i==3;
    plot(-10:1:10,mean(squeeze(mean(Xcorr10Lag(:,sorted3inds{i},:),2)),1),'k');
    title(direction3names{i});
    ylabel('\rho(k)');xlabel({'lag'; '(L-leading <- 0 -> R-leading)'});ylim([-0.1 0.8])
    end
    yline(0,'color',[1 0.8 0.2]);xline(0,'color',[1 0.8 0.2]);
    grid on; grid minor
end
% sgtitle('mean xcorr sum for 10 lags in all sessions (matched int) ^{ *PLOT 9-1}');
delete(findall(gcf,'type','annotation'))
set([sg], 'fitboxtotext','on',...
    'edgecolor','none')
sg=annotation('textbox',[0.01 0.01 0.07 0.05],'string',...
    'Xcorr ^{* PLOT 9-1}')
set(gcf,'color','w'); % set background white for copying in ubuntu

% separate synch and synco with diff color (1 x 3 direction) (published)
Xcorr10Lag; % sorted order
direction3names={'Independent','Unidirectional','Bidirectional'};
sorted3inds={[1:3],[4:9],[10:12]};
sorted4inds=[1:3; 4:6; 7:9; 10:12];
syn2Ind;
canvas(0.45, 0.37); % figure;
for i=1:3
    for syn=1:2
    subplot(1,3,i)
    if i==2;
        % L-lead
        L_mat=squeeze(mean(Xcorr10Lag(syn2Ind{syn},sorted4inds(2,:),:),2));
        % R-lead
        R_mat=squeeze(mean(Xcorr10Lag(syn2Ind{syn},sorted4inds(3,:),:),2));
        % combine
        LR_mean=mean([L_mat;fliplr(R_mat)]);
        % plot Leader vs Follower
        hold on;
        plot(-10:1:10,LR_mean,'linewidth',4,'color',[syn2colors(syn,:) 1]);
        ax = gca;
        ax.FontSize = 15; 
        title(direction3names{i},'FontSize',15);
        ylabel('\rho','FontSize',20,'FontWeight','bold');
        xlabel({'Lag','    Leader -leading <- 0 -> Follower - following'},'FontSize',15);ylim([-0.1 0.8])
    elseif i==1 | i==3;
        hold on;
        plot(-10:1:10,mean(squeeze(mean(Xcorr10Lag(syn2Ind{syn},sorted3inds{i},:),2)),1),'linewidth',4,'color',[syn2colors(syn,:) 1]);
        ax = gca;
        ax.FontSize = 15; 
        title(direction3names{i},'FontSize',17);
        ylabel('\rho','FontSize',20,'FontWeight','bold');xlabel({'Lag', 'A -leading <- 0 -> B -leading'},'FontSize',17);ylim([-0.1 0.8])
    end    
    hold off;
    end
    yline(0,'color',[1 0.8 0.2]);xline(0,'color',[1 0.8 0.2]);
    xline(-1,'color',[1 0.8 0.2]);xline(1,'color',[1 0.8 0.2]);
    lg=legend({'synch','synco','','','',''});
    lg.FontSize=17;
end
% lg=legend('synch','synco','location','eastoutside');
% sgtitle('mean xcorr sum for 10 lags in synch/o sessions (matched int) ^{ *PLOT 9-1}');
set(gcf,'color','w'); % set background white for copying in ubuntu
%% SECT 10 redo H and intervals_H_removed and save them separately (matched int)
% in original order
clear
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);
sessions={'synch','synco','synch','synco','synch','synco','synch','synco','synch','synco',...
    'synch','synco'};
cd /ssd/zhibin/1overf/
H=[];intervals_H_removed=[];
tic
for s=1:numSes
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
toc % 98 sec
%% SECT 10-1 extract H and condition from all sessions & 4 states (matched int)
% original order
% clear
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);
sessions={'synch','synco','synch','synco','synch','synco','synch','synco','synch','synco',...
    'synch','synco'};
cd /ssd/zhibin/1overf/
H_all=[];
for s=1:numSes
    runid=num2str(seeds(s,:));
    clear H
    load(['/ssd/zhibin/1overf/' runid '_2P/Cleaned_data/d_removal' runid '.mat' ])   
    for b=1:12
        H_all(:,s,b)=H(:,b); % 2 subject x 12 session x 12 trials
    end
end

% organize conditions in all sessions into a vector of 192
condition_all=[];
tic
for s=1:numSes
    clear conditions
    runid=num2str(seeds(s,:));
    load(['/ssd/zhibin/1overf/' runid '_2P/Cleaned_data/clean_' runid '.mat'],'conditions');
    condition_all(s,:)=conditions; % 12 session x 12 trials
end
toc % 78 sec
clear dataR dataL
% reshape into a vector in time sequence
condition_all=reshape(condition_all',[],1); % 144 x 1 

% indices for 4 states
states4names={'Uncoupled','Leading','Following','Mutual'};
% find the indices for each condition in L R conbined sequence (4 states)
% and organize for PLS
uncoupleInd_LR=[find(condition_all==1);12*numSes+find(condition_all==1)];
leadingInd_LR=[find(condition_all==2);12*numSes+find(condition_all==3)];
followingInd_LR=[find(condition_all==3);12*numSes+find(condition_all==2)];
mutualInd_LR=[find(condition_all==4);12*numSes+find(condition_all==4)];
Inds4_LR=[uncoupleInd_LR leadingInd_LR followingInd_LR mutualInd_LR];

% H-int matching for L and R as in Marmelat 2012
% organize H_all 
H_all; % 2x12numSesx12trials
H_Lall=squeeze(H_all(1,:,:)); % 12session by 12trials
H_Lall=reshape(H_Lall',[],1); % in time sequence
H_Rall=squeeze(H_all(2,:,:)); % 12session by 12trials
H_Rall=reshape(H_Rall',[],1); % in time sequence
H_all_LR=[H_Lall; H_Rall]
% Inds4(:,c);
% uncoupleInd(1:24);
% leadingInd(1:24); % for L
% followingInd(1:24); % for L
% mutualInd(1:24);
%% SECT 10-2 Construct testing dataset
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);

% organize conditions in all sessions into a vector of 288
condition_all=[]; % 12 session x 12 trials
% initialize the testing dataset
testing_data=[]; % 12 session x 12 trials
tic
for s=1:numSes
    clear conditions
    runid=num2str(seeds(s,:));
    load(['/ssd/zhibin/1overf/' runid '_2P/Cleaned_data/clean_' runid '.mat' ]);
    condition_all(s,:)=conditions; 
    % tesing_data has the following structure: L&R assinged same value
    % synco are negaitive and synch are positive
    % sessions marked with 0.01 to 0.12
    if mod(s,2)==0;
        % assigned negative to synco condi
        testing_data(s,:)=-1*conditions-s/100; 
    else mod(s,2)==1;
        % append sequential integer as two decimals for each pair of subject in each session
        testing_data(s,:)=conditions+s/100; 
    end
end
toc % took 78 sec

% reshape into a vector in time sequence
condition_all=reshape(condition_all',[],1); % 144 x 1 
testing_data_L=reshape(testing_data',[],1); % 144 x 1 
testing_data_R=testing_data_L; % 144 x 1 
testing_data_LR=[testing_data_L; testing_data_R];
% testing EEG dataset in 32channels
% same structure as testing data but mark each channel with 0.0001-0.0032
channel_mark1=nan(12,32);
for s=1:numSes
    if mod(s,2)==0;
        % assigned negative to synco condi
        channel_mark1(s,:)=-[1:32]/10000; 
    else mod(s,2)==1;
        % append sequential integer as two decimals for each pair of subject in each session
        channel_mark1(s,:)=[1:32]/10000; 
    end
end
channel_mark2=nan(144,32);
for s=1:numSes
    channel_mark2([1:12]+12*(s-1),:)=repmat(channel_mark1(s,:),12,1);
end
testing_EEG_L=repmat(testing_data_L,1,32) + channel_mark2;
testing_EEG_R=repmat(testing_data_R,1,32) + channel_mark2;
testing_EEG_LR=[testing_EEG_L;testing_EEG_R];
%% PLOT 10-2 H matching as in Marmelat 2012 (matched int)
% H-int matching for L and R as in Marmelat 2012
% indicies in the synch and synco time sequence
synchind=[1:3 7:9 13:15 19:21 25:27 31:33]; % 3 trials x 6 sessions
syncoind=[4:6 10:12 16:18 22:24 28:30 34:36]; % 3 trials x 6 sessions
synind=[synchind;syncoind];

%********************** 2 subplots (uncouple and mutual)
for  subplots2 = 1;
canvas(0.2,0.25);
subplot(1,2,1); %##############
plot(H_Lall(uncoupleInd_LR(synchind)),H_Rall(uncoupleInd_LR(synchind)),'.','MarkerSize',30,'color',darkgreen);
hold on;
plot(H_Lall(uncoupleInd_LR(syncoind)),H_Rall(uncoupleInd_LR(syncoind)),'.','MarkerSize',30,'color',pink);
xlabel('DFA exponent, Participant L');ylabel('DFA exponent, Participant R');
title('Uncouple');
A=[];Alpha1=[];FitValues=[];
A=polyfit(H_Lall(uncoupleInd_LR(1:36)),H_Rall(uncoupleInd_LR(1:36)),1);
Alpha1=A(1); % the slope, the first order polynomial coefficient from polyfit (Hurst Componenet >1 ?)
FitValues=polyval(A,H_Lall(uncoupleInd_LR(1:36)));
hold on; plot(H_Lall(uncoupleInd_LR(1:36)),FitValues,'k-');
xlim([0.2 1.4]);ylim([0.2 1.4]);plot([0 1.4], [0 1.4],'m--');hold off;
grid on;
legend('synch','synco','','','location','northwest');
subplot(1,2,2); %##############
plot(H_Lall(mutualInd_LR(synchind)),H_Rall(mutualInd_LR(synchind)),'.','MarkerSize',30,'color',darkgreen);
hold on;
plot(H_Lall(mutualInd_LR(syncoind)),H_Rall(mutualInd_LR(syncoind)),'.','MarkerSize',30,'color',pink);
xlabel('DFA exponent, Participant L');ylabel('DFA exponent, Participant R');
title('Mutual');
A=[];Alpha1=[];FitValues=[];
A=polyfit(H_Lall(mutualInd_LR(1:36)),H_Rall(mutualInd_LR(1:36)),1);
Alpha1=A(1); 
FitValues=polyval(A,H_Lall(mutualInd_LR(1:36)));
hold on; plot(H_Lall(mutualInd_LR(1:36)),FitValues,'k-');
xlim([0.2 1.4]);ylim([0.2 1.4]);plot([0 1.4], [0 1.4],'m--');hold off;
grid on;
legend('synch','synco','','','location','northwest');
end

%********************** 4 subplots (4 states)
for subplots4 = 1;
canvas(0.2,0.5);
subplot(2,2,1); %##############
plot(H_Lall(uncoupleInd_LR(synchind)),H_Rall(uncoupleInd_LR(synchind)),'.','MarkerSize',30,'color',darkgreen);
hold on;
plot(H_Lall(uncoupleInd_LR(syncoind)),H_Rall(uncoupleInd_LR(syncoind)),'.','MarkerSize',30,'color',pink);
xlabel('DFA exponent, Participant L');ylabel('DFA exponent, Participant R');
title('Uncouple');
A=[];Alpha1=[];FitValues=[];
A=polyfit(H_Lall(uncoupleInd_LR(1:36)),H_Rall(uncoupleInd_LR(1:36)),1);
Alpha1=A(1); % the slope, the first order polynomial coefficient from polyfit (Hurst Componenet >1 ?)
FitValues=polyval(A,H_Lall(uncoupleInd_LR(1:36)));
hold on; plot(H_Lall(uncoupleInd_LR(1:36)),FitValues,'k-');
xlim([0.2 1.4]);ylim([0.2 1.4]);plot([0 1.4], [0 1.4],'m--');hold off;
legend('synch','synco','','','location','northwest');
grid on;
subplot(2,2,2); %##############
plot(H_Lall(leadingInd_LR(synchind)),H_Rall(leadingInd_LR(synchind)),'.','MarkerSize',30,'color',darkgreen);
hold on;
plot(H_Lall(leadingInd_LR(syncoind)),H_Rall(leadingInd_LR(syncoind)),'.','MarkerSize',30,'color',pink);
xlabel('DFA exponent, Participant L');ylabel('DFA exponent, Participant R');
title('L-lead');
A=[];Alpha1=[];FitValues=[];
A=polyfit(H_Lall(leadingInd_LR(1:36)),H_Rall(leadingInd_LR(1:36)),1);
Alpha1=A(1); % the slope, the first order polynomial coefficient from polyfit (Hurst Componenet >1 ?)
FitValues=polyval(A,H_Lall(leadingInd_LR(1:36)));
hold on; plot(H_Lall(leadingInd_LR(1:36)),FitValues,'k-');
xlim([0.2 1.4]);ylim([0.2 1.4]);plot([0 1.4], [0 1.4],'m--');hold off;
legend('synch','synco','','','location','northwest');
grid on;
subplot(2,2,3); %##############
plot(H_Lall(followingInd_LR(synchind)),H_Rall(followingInd_LR(synchind)),'.','MarkerSize',30,'color',darkgreen);
hold on;
plot(H_Lall(followingInd_LR(syncoind)),H_Rall(followingInd_LR(syncoind)),'.','MarkerSize',30,'color',pink);
xlabel('DFA exponent, Participant L');ylabel('DFA exponent, Participant R');
title('R-lead');
A=[];Alpha1=[];FitValues=[];
A=polyfit(H_Lall(followingInd_LR(1:36)),H_Rall(followingInd_LR(1:36)),1);
Alpha1=A(1); % the slope, the first order polynomial coefficient from polyfit (Hurst Componenet >1 ?)
FitValues=polyval(A,H_Lall(followingInd_LR(1:36)));
hold on; plot(H_Lall(followingInd_LR(1:36)),FitValues,'k-');
xlim([0.2 1.4]);ylim([0.2 1.4]);plot([0 1.4], [0 1.4],'m--');hold off;
legend('synch','synco','','','location','northwest');
grid on;
subplot(2,2,4); %##############
plot(H_Lall(mutualInd_LR(synchind)),H_Rall(mutualInd_LR(synchind)),'.','MarkerSize',30,'color',darkgreen);
hold on;
plot(H_Lall(mutualInd_LR(syncoind)),H_Rall(mutualInd_LR(syncoind)),'.','MarkerSize',30,'color',pink);
xlabel('DFA exponent, Participant L');ylabel('DFA exponent, Participant R');
title('mutual');
A=[];Alpha1=[];FitValues=[];
A=polyfit(H_Lall(mutualInd_LR(1:36)),H_Rall(mutualInd_LR(1:36)),1);
Alpha1=A(1); 
FitValues=polyval(A,H_Lall(mutualInd_LR(1:36)));
hold on; plot(H_Lall(mutualInd_LR(1:36)),FitValues,'k-');
xlim([0.2 1.4]);ylim([0.2 1.4]);plot([0 1.4], [0 1.4],'m--');hold off;
legend('synch','synco','','','location','northwest');
grid on;
end

%********************** 3 subplots
for subplots3 =1;
figure; % canvas(0.5,0.35);
subplot(1,3,1); %##############
for subplot1 =1;
plot(H_Lall(uncoupleInd_LR(synchind(1:3))),H_Rall(uncoupleInd_LR(synchind(1:3))),'.','MarkerSize',30,'color',darkgreen);hold on;
plot(H_Lall(uncoupleInd_LR(synchind(4:6))),H_Rall(uncoupleInd_LR(synchind(4:6))),'square','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(uncoupleInd_LR(synchind(7:9))),H_Rall(uncoupleInd_LR(synchind(7:9))),'^','MarkerSize',8,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(uncoupleInd_LR(synchind(10:12))),H_Rall(uncoupleInd_LR(synchind(10:12))),'pentagram','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(uncoupleInd_LR(synchind(13:15))),H_Rall(uncoupleInd_LR(synchind(13:15))),'*','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(uncoupleInd_LR(synchind(16:18))),H_Rall(uncoupleInd_LR(synchind(16:18))),'diamond','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(uncoupleInd_LR(syncoind(1:3))),H_Rall(uncoupleInd_LR(syncoind(1:3))),'.','MarkerSize',30,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(uncoupleInd_LR(syncoind(4:6))),H_Rall(uncoupleInd_LR(syncoind(4:6))),'square','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(uncoupleInd_LR(syncoind(7:9))),H_Rall(uncoupleInd_LR(syncoind(7:9))),'^','MarkerSize',8,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(uncoupleInd_LR(syncoind(10:12))),H_Rall(uncoupleInd_LR(syncoind(10:12))),'pentagram','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(uncoupleInd_LR(syncoind(13:15))),H_Rall(uncoupleInd_LR(syncoind(13:15))),'*','MarkerSize',8,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(uncoupleInd_LR(syncoind(16:18))),H_Rall(uncoupleInd_LR(syncoind(16:18))),'diamond','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
xlabel('DFA exponent, Participant L','FontSize',15);
ylabel('DFA exponent, Participant R','FontSize',15);
title('Independent','FontSize',15);
% A=[];S=[]; Alpha1=[];FitValues=[];
% A=polyfit(H_Lall(uncoupleInd_LR(1:36)),H_Rall(uncoupleInd_LR_LR(1:36)),1);
% Alpha1=A(1); % the slope, the first order polynomial coefficient from polyfit (Hurst Componenet >1 ?)
% FitValues=polyval(A,H_Lall(uncoupleInd_LR(1:36)));
% hold on; plot(H_Lall(uncoupleInd_LR(1:36)),FitValues,'k-','LineWidth',3);

% fit the synch data
A=[];Alpha1=[];FitValues=[];RHO=[];xxx=[];yyy=[];
[xxx,yyy]=keepvalues(H_Lall(uncoupleInd_LR(1:2:35)),H_Rall(uncoupleInd_LR(1:2:35)),0.5,1);
A=polyfit(xxx,yyy,1);
[RHO,~]=corr(xxx,yyy);
Alpha1=A(1); % the slope, the first order polynomial coefficient from polyfit (Hurst Componenet >1 ?)
FitValues=polyval(A,xxx);
hold on; plot(xxx,FitValues,'-','color',darkgreen,'LineWidth',3);
text(min(xxx),max(FitValues)+0.025,sprintf('\\rho=%.2f',RHO),'Color',darkgreen,'FontSize', 15)

% fit the synco data
A=[];Alpha1=[];FitValues=[];RHO=[];xxx=[];yyy=[];
[xxx,yyy]=keepvalues(H_Lall(uncoupleInd_LR(2:2:36)),H_Rall(uncoupleInd_LR(2:2:36)),0.5,1);
A=polyfit(xxx,yyy,1);
[RHO,~]=corr(xxx,yyy);
Alpha1=A(1); % the slope, the first order polynomial coefficient from polyfit (Hurst Componenet >1 ?)
FitValues=polyval(A,xxx);
hold on; plot(xxx,FitValues,'-','color',pink,'LineWidth',3);
xlim([0.5 1]);ylim([0.5 1]);
plot([0.5 1], [0.5 1],'m--');
hold off;
text(min(xxx),min(FitValues)-0.025,sprintf('\\rho=%.2f',RHO),'Color',pink,'FontSize', 15)
% legend('synch','','','','synco','','','','location','northwest');
% subtitle(['(uncouple: linear regression slope = ' num2str(Alpha1,'%.3f') ')'])
grid on;
end
subplot(1,3,2); %##############
for subplot2 =1;
plot(H_Lall(leadingInd_LR(synchind(1:3))),H_Rall(leadingInd_LR(synchind(1:3))),'.','MarkerSize',30,'color',darkgreen); hold on;
plot(H_Lall(leadingInd_LR(synchind(4:6))),H_Rall(leadingInd_LR(synchind(4:6))),'square','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(leadingInd_LR(synchind(7:9))),H_Rall(leadingInd_LR(synchind(7:9))),'^','MarkerSize',8,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(leadingInd_LR(synchind(10:12))),H_Rall(leadingInd_LR(synchind(10:12))),'pentagram','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(leadingInd_LR(synchind(13:15))),H_Rall(leadingInd_LR(synchind(13:15))),'*','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(leadingInd_LR(synchind(16:18))),H_Rall(leadingInd_LR(synchind(16:18))),'diamond','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(leadingInd_LR(syncoind(1:3))),H_Rall(leadingInd_LR(syncoind(1:3))),'.','MarkerSize',30,'color',pink);
plot(H_Lall(leadingInd_LR(syncoind(4:6))),H_Rall(leadingInd_LR(syncoind(4:6))),'square','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(leadingInd_LR(syncoind(7:9))),H_Rall(leadingInd_LR(syncoind(7:9))),'^','MarkerSize',8,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(leadingInd_LR(syncoind(10:12))),H_Rall(leadingInd_LR(syncoind(10:12))),'pentagram','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(leadingInd_LR(syncoind(13:15))),H_Rall(leadingInd_LR(syncoind(13:15))),'*','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(leadingInd_LR(syncoind(16:18))),H_Rall(leadingInd_LR(syncoind(16:18))),'diamond','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Rall(followingInd_LR(synchind(1:3))),H_Lall(followingInd_LR(synchind(1:3))),'.','MarkerSize',30,'color',darkgreen);
plot(H_Rall(followingInd_LR(synchind(4:6))),H_Lall(followingInd_LR(synchind(4:6))),'square','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Rall(followingInd_LR(synchind(7:9))),H_Lall(followingInd_LR(synchind(7:9))),'^','MarkerSize',8,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Rall(followingInd_LR(synchind(10:12))),H_Lall(followingInd_LR(synchind(10:12))),'pentagram','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Rall(followingInd_LR(synchind(13:15))),H_Lall(followingInd_LR(synchind(13:15))),'*','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Rall(followingInd_LR(synchind(16:18))),H_Lall(followingInd_LR(synchind(16:18))),'diamond','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Rall(followingInd_LR(syncoind(1:3))),H_Lall(followingInd_LR(syncoind(1:3))),'.','MarkerSize',30,'color',pink);
plot(H_Rall(followingInd_LR(syncoind(4:6))),H_Lall(followingInd_LR(syncoind(4:6))),'square','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Rall(followingInd_LR(syncoind(7:9))),H_Lall(followingInd_LR(syncoind(7:9))),'^','MarkerSize',8,'color',pink,'MarkerFaceColor',pink);
plot(H_Rall(followingInd_LR(syncoind(10:12))),H_Lall(followingInd_LR(syncoind(10:12))),'pentagram','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Rall(followingInd_LR(syncoind(13:15))),H_Lall(followingInd_LR(syncoind(13:15))),'pentagram','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Rall(followingInd_LR(syncoind(16:18))),H_Lall(followingInd_LR(syncoind(16:18))),'pentagram','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
xlabel('DFA exponent, Leader','FontSize',15);
ylabel('DFA exponent, Follower','FontSize',15);
title('Unidirectional','FontSize',15);

% A=[];Alpha1=[];FitValues=[];
% A=polyfit([H_Lall(leadingInd_LR(1:36)); H_Rall(followingInd_LR(1:36))], [H_Rall(leadingInd_LR(1:36)); H_Lall(followingInd_LR(1:36))],1);
% Alpha1=A(1); % the slope, the first order polynomial coefficient from polyfit (Hurst Componenet >1 ?)
% FitValues=polyval(A,[H_Lall(leadingInd_LR(1:36)); H_Rall(followingInd_LR(1:36))]);
% hold on; plot([H_Lall(leadingInd_LR(1:36)); H_Rall(followingInd_LR(1:36))],FitValues,'k-','LineWidth',3);
% fit the synch data
A=[];Alpha1=[];FitValues=[];RHO=[];xxx=[];yyy=[];
[xxx,yyy]=keepvalues([H_Lall(leadingInd_LR(1:2:35)); H_Rall(followingInd_LR(1:2:35))], ...
    [H_Rall(leadingInd_LR(1:2:35)); H_Lall(followingInd_LR(1:2:35))],0.5,1);
A=polyfit(xxx, yyy,1);
[RHO,~]=corr(xxx, yyy);
Alpha1=A(1); % the slope, the first order polynomial coefficient from polyfit (Hurst Componenet >1 ?)
FitValues=polyval(A,xxx);
hold on; 
plot(xxx,FitValues,'-','color',darkgreen,'LineWidth',3);
text(max(xxx)-0.1,max(FitValues)+0.02,sprintf('\\rho=%.2f',RHO),'Color',darkgreen,'FontSize', 15)
% fit the synco data
A=[];Alpha1=[];FitValues=[];RHO=[];xxx=[];yyy=[];
[xxx,yyy]=keepvalues([H_Lall(leadingInd_LR(2:2:36)); H_Rall(followingInd_LR(2:2:36))],...
    [H_Rall(leadingInd_LR(2:2:36)); H_Lall(followingInd_LR(2:2:36))],0.5,1);
A=polyfit(xxx, yyy,1);
[RHO,~]=corr(xxx, yyy);
Alpha1=A(1); % the slope, the first order polynomial coefficient from polyfit (Hurst Componenet >1 ?)
FitValues=polyval(A,xxx);
hold on; plot(xxx,FitValues,'-','color',pink,'LineWidth',3);
text(max(xxx)-0.1,max(FitValues)-0.03,sprintf('\\rho=%.2f',RHO),'Color',pink,'FontSize', 15)
xlim([0.5 1]);ylim([0.5 1]);
plot([0.5 1], [0.5 1],'m--');
hold off;
% legend('synch','','','','synco','','','','','','','','','','','','location','northwest');
% subtitle(['(leading or following: linear regression slope = ' num2str(Alpha1,'%.3f') ')'])
grid on;
end
subplot(1,3,3); %##############
for subplot3 =1;
plot(H_Lall(mutualInd_LR(synchind(1:3))),H_Rall(mutualInd_LR(synchind(1:3))),'.','MarkerSize',30,'color',darkgreen);hold on;
plot(H_Lall(mutualInd_LR(synchind(4:6))),H_Rall(mutualInd_LR(synchind(4:6))),'square','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(mutualInd_LR(synchind(7:9))),H_Rall(mutualInd_LR(synchind(7:9))),'^','MarkerSize',8,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(mutualInd_LR(synchind(10:12))),H_Rall(mutualInd_LR(synchind(10:12))),'pentagram','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(mutualInd_LR(synchind(13:15))),H_Rall(mutualInd_LR(synchind(13:15))),'*','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(mutualInd_LR(synchind(16:18))),H_Rall(mutualInd_LR(synchind(16:18))),'diamond','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(mutualInd_LR(syncoind(1:3))),H_Rall(mutualInd_LR(syncoind(1:3))),'.','MarkerSize',30,'color',pink);
plot(H_Lall(mutualInd_LR(syncoind(4:6))),H_Rall(mutualInd_LR(syncoind(4:6))),'square','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(mutualInd_LR(syncoind(7:9))),H_Rall(mutualInd_LR(syncoind(7:9))),'^','MarkerSize',8,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(mutualInd_LR(syncoind(10:12))),H_Rall(mutualInd_LR(syncoind(10:12))),'pentagram','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(mutualInd_LR(syncoind(13:15))),H_Rall(mutualInd_LR(syncoind(13:15))),'*','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(mutualInd_LR(syncoind(16:18))),H_Rall(mutualInd_LR(syncoind(16:18))),'diamond','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
xlabel('DFA exponent, Participant L','FontSize',15);
ylabel('DFA exponent, Participant R','FontSize',15);
title('Bidirectional','FontSize',15); 
xlim([0.5 1]);ylim([0.5 1]);

% fit all the data
A=[];Alpha1=[];FitValues=[];RHO=[];xxx=[];yyy=[];
[xxx,yyy]=keepvalues(H_Lall(mutualInd_LR(1:36)),...
    H_Rall(mutualInd_LR(1:36)),0.5,1);
A=polyfit(xxx,yyy,1);
[RHO,~]=corr(xxx,yyy);
Alpha1=A(1); 
FitValues=polyval(A,xxx);
hold on; plot(xxx,FitValues,'k-','LineWidth',3);
text(min(xxx),min(FitValues)-0.025,sprintf('\\rho=%.2f',RHO),'Color',[0 0 0],'FontSize', 15)
% fit the synch data
A=[];Alpha1=[];FitValues=[];RHO=[];xxx=[];yyy=[];
[xxx,yyy]=keepvalues(H_Lall(mutualInd_LR(1:2:35)),...
    H_Rall(mutualInd_LR(1:2:35)),0.5,1);
A=polyfit(xxx,yyy,1);
[RHO,~]=corr(xxx,yyy);
Alpha1=A(1); % the slope, the first order polynomial coefficient from polyfit (Hurst Componenet >1 ?)
FitValues=polyval(A,xxx);
hold on; plot(xxx,FitValues,'-','color',darkgreen,'LineWidth',3);
text(min(xxx),min(FitValues)+0.07,sprintf('\\rho=%.2f',RHO),'Color',darkgreen,'FontSize', 15)
% fit the synco data
A=[];Alpha1=[];FitValues=[];RHO=[];xxx=[];yyy=[];
[xxx,yyy]=keepvalues(H_Lall(mutualInd_LR(2:2:36)),...
    H_Rall(mutualInd_LR(2:2:36)),0.5,1);
A=polyfit(xxx,yyy,1);
[RHO,~]=corr(xxx,yyy);
Alpha1=A(1); % the slope, the first order polynomial coefficient from polyfit (Hurst Componenet >1 ?)
FitValues=polyval(A,xxx);
hold on; plot(xxx,FitValues,'-','color',pink,'LineWidth',3);
text(max(xxx)-0.15,max(FitValues)-0.025,sprintf('\\rho=%.2f',RHO),'Color',pink,'FontSize', 15)
xlim([0.5 1]);ylim([0.5 1]);
plot([0.5 1], [0.5 1],'m--');
hold off;
% subtitle(['(mutual): linear regression slope = ' num2str(Alpha1,'%.3f') ')'])
grid on;
end
lg=legend('synch - subj pair 1','synch - subj pair 2','synch - subj pair 3','synch - subj pair 4', ...
    'synch - subj pair 5','synch - subj pair 6', ...
    'synco - subj pair 1','synco - subj pair 2','synco - subj pair 3','synco - subj pair 4',...
    'synco - subj pair 5','synco - subj pair 6','location','eastoutside');
lg.Position = [0.9475 0.15 0.01 0.75];
% sgtitle(['H matching as in Marmelat 2012  ^{* PLOT 10-2}']);
delete(findall(gcf,'type','annotation'))
sg=annotation('textbox',[0.07 0.01 0.22 0.05],'string',...
    'H matching as in Marmelat 2012  ^{* PLOT 10-2}')
sg.Rotation=90
set(gcf,'color','w'); % set background white for copying in ubuntu
end

%% PLOT 10-2 H matching as in Marmelat 2012 (matched int) Reproduced the plot in poster
for subplots3 =1 
figure; % canvas(0.5,0.35);
subplot(1,3,1); %##############
for subplot1 =1;
plot(H_Lall(uncoupleInd_LR(synchind(1:3))),H_Rall(uncoupleInd_LR(synchind(1:3))),'.','MarkerSize',30,'color',darkgreen);hold on;
plot(H_Lall(uncoupleInd_LR(synchind(4:6))),H_Rall(uncoupleInd_LR(synchind(4:6))),'square','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(uncoupleInd_LR(synchind(7:9))),H_Rall(uncoupleInd_LR(synchind(7:9))),'^','MarkerSize',8,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(uncoupleInd_LR(synchind(10:12))),H_Rall(uncoupleInd_LR(synchind(10:12))),'pentagram','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(uncoupleInd_LR(synchind(13:15))),H_Rall(uncoupleInd_LR(synchind(13:15))),'*','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(uncoupleInd_LR(synchind(16:18))),H_Rall(uncoupleInd_LR(synchind(16:18))),'diamond','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(uncoupleInd_LR(syncoind(1:3))),H_Rall(uncoupleInd_LR(syncoind(1:3))),'.','MarkerSize',30,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(uncoupleInd_LR(syncoind(4:6))),H_Rall(uncoupleInd_LR(syncoind(4:6))),'square','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(uncoupleInd_LR(syncoind(7:9))),H_Rall(uncoupleInd_LR(syncoind(7:9))),'^','MarkerSize',8,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(uncoupleInd_LR(syncoind(10:12))),H_Rall(uncoupleInd_LR(syncoind(10:12))),'pentagram','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(uncoupleInd_LR(syncoind(13:15))),H_Rall(uncoupleInd_LR(syncoind(13:15))),'*','MarkerSize',8,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(uncoupleInd_LR(syncoind(16:18))),H_Rall(uncoupleInd_LR(syncoind(16:18))),'diamond','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
ax = gca;
ax.FontSize = 13; 
xlabel('DFA Hurst Exponent, Participant A','FontSize',15);
ylabel('DFA Husrt Exponent, Participant B','FontSize',15);
title('Independent','FontSize',15);


% A=[];S=[]; Alpha1=[];FitValues=[];
% A=polyfit(H_Lall(uncoupleInd_LR(1:36)),H_Rall(uncoupleInd_LR_LR(1:36)),1);
% Alpha1=A(1); % the slope, the first order polynomial coefficient from polyfit (Hurst Componenet >1 ?)
% FitValues=polyval(A,H_Lall(uncoupleInd_LR(1:36)));
% hold on; plot(H_Lall(uncoupleInd_LR(1:36)),FitValues,'k-','LineWidth',3);

% fit the synch data
A=[];Alpha1=[];FitValues=[];RHO=[];xxx=[];yyy=[];
[xxx,yyy]=keepvalues(H_Lall(uncoupleInd_LR(1:2:35)),H_Rall(uncoupleInd_LR(1:2:35)),0.2,1.4);
A=polyfit(xxx,yyy,1);
[RHO,~]=corr(xxx,yyy);
Alpha1=A(1); % the slope, the first order polynomial coefficient from polyfit (Hurst Componenet >1 ?)
FitValues=polyval(A,xxx);
hold on; plot(xxx,FitValues,'-','color',darkgreen,'LineWidth',3);
text(max(xxx),min(FitValues),sprintf('\\rho=%.2f',RHO),'Color',darkgreen,'FontSize', 15)

% fit the synco data
A=[];Alpha1=[];FitValues=[];RHO=[];xxx=[];yyy=[];
[xxx,yyy]=keepvalues(H_Lall(uncoupleInd_LR(2:2:36)),H_Rall(uncoupleInd_LR(2:2:36)),0.2,1.4);
A=polyfit(xxx,yyy,1);
[RHO,~]=corr(xxx,yyy);
Alpha1=A(1); % the slope, the first order polynomial coefficient from polyfit (Hurst Componenet >1 ?)
FitValues=polyval(A,xxx);
hold on; plot(xxx,FitValues,'-','color',pink,'LineWidth',3);
% xlim([0.5 1]);ylim([0.5 1]);
xlim([0.2 1.4]);ylim([0.2 1.4]);
plot([0.5 1], [0.5 1],'m--');
hold off;
text(max(xxx),max(FitValues),sprintf('\\rho=%.2f',RHO),'Color',pink,'FontSize', 15)
% legend('synch','','','','synco','','','','location','northwest');
% subtitle(['(uncouple: linear regression slope = ' num2str(Alpha1,'%.3f') ')'])
grid on;
end
subplot(1,3,2); %##############
for subplot2 =1;
plot(H_Lall(leadingInd_LR(synchind(1:3))),H_Rall(leadingInd_LR(synchind(1:3))),'.','MarkerSize',30,'color',darkgreen); hold on;
plot(H_Lall(leadingInd_LR(synchind(4:6))),H_Rall(leadingInd_LR(synchind(4:6))),'square','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(leadingInd_LR(synchind(7:9))),H_Rall(leadingInd_LR(synchind(7:9))),'^','MarkerSize',8,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(leadingInd_LR(synchind(10:12))),H_Rall(leadingInd_LR(synchind(10:12))),'pentagram','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(leadingInd_LR(synchind(13:15))),H_Rall(leadingInd_LR(synchind(13:15))),'*','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(leadingInd_LR(synchind(16:18))),H_Rall(leadingInd_LR(synchind(16:18))),'diamond','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(leadingInd_LR(syncoind(1:3))),H_Rall(leadingInd_LR(syncoind(1:3))),'.','MarkerSize',30,'color',pink);
plot(H_Lall(leadingInd_LR(syncoind(4:6))),H_Rall(leadingInd_LR(syncoind(4:6))),'square','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(leadingInd_LR(syncoind(7:9))),H_Rall(leadingInd_LR(syncoind(7:9))),'^','MarkerSize',8,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(leadingInd_LR(syncoind(10:12))),H_Rall(leadingInd_LR(syncoind(10:12))),'pentagram','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(leadingInd_LR(syncoind(13:15))),H_Rall(leadingInd_LR(syncoind(13:15))),'*','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(leadingInd_LR(syncoind(16:18))),H_Rall(leadingInd_LR(syncoind(16:18))),'diamond','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Rall(followingInd_LR(synchind(1:3))),H_Lall(followingInd_LR(synchind(1:3))),'.','MarkerSize',30,'color',darkgreen);
plot(H_Rall(followingInd_LR(synchind(4:6))),H_Lall(followingInd_LR(synchind(4:6))),'square','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Rall(followingInd_LR(synchind(7:9))),H_Lall(followingInd_LR(synchind(7:9))),'^','MarkerSize',8,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Rall(followingInd_LR(synchind(10:12))),H_Lall(followingInd_LR(synchind(10:12))),'pentagram','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Rall(followingInd_LR(synchind(13:15))),H_Lall(followingInd_LR(synchind(13:15))),'*','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Rall(followingInd_LR(synchind(16:18))),H_Lall(followingInd_LR(synchind(16:18))),'diamond','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Rall(followingInd_LR(syncoind(1:3))),H_Lall(followingInd_LR(syncoind(1:3))),'.','MarkerSize',30,'color',pink);
plot(H_Rall(followingInd_LR(syncoind(4:6))),H_Lall(followingInd_LR(syncoind(4:6))),'square','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Rall(followingInd_LR(syncoind(7:9))),H_Lall(followingInd_LR(syncoind(7:9))),'^','MarkerSize',8,'color',pink,'MarkerFaceColor',pink);
plot(H_Rall(followingInd_LR(syncoind(10:12))),H_Lall(followingInd_LR(syncoind(10:12))),'pentagram','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Rall(followingInd_LR(syncoind(13:15))),H_Lall(followingInd_LR(syncoind(13:15))),'pentagram','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Rall(followingInd_LR(syncoind(16:18))),H_Lall(followingInd_LR(syncoind(16:18))),'pentagram','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
ax = gca;
ax.FontSize = 13; 
xlabel('DFA Hurst Exponent, Leader','FontSize',15);
ylabel('DFA Hurst Exponent, Follower','FontSize',15);
title('Unidirectional','FontSize',15);

% A=[];Alpha1=[];FitValues=[];
% A=polyfit([H_Lall(leadingInd_LR(1:36)); H_Rall(followingInd_LR(1:36))], [H_Rall(leadingInd_LR(1:36)); H_Lall(followingInd_LR(1:36))],1);
% Alpha1=A(1); % the slope, the first order polynomial coefficient from polyfit (Hurst Componenet >1 ?)
% FitValues=polyval(A,[H_Lall(leadingInd_LR(1:36)); H_Rall(followingInd_LR(1:36))]);
% hold on; plot([H_Lall(leadingInd_LR(1:36)); H_Rall(followingInd_LR(1:36))],FitValues,'k-','LineWidth',3);

% fit the synch data
A=[];Alpha1=[];FitValues=[];RHO=[];xxx=[];yyy=[];
[xxx,yyy]=keepvalues([H_Lall(leadingInd_LR(1:2:35)); H_Rall(followingInd_LR(1:2:35))], ...
    [H_Rall(leadingInd_LR(1:2:35)); H_Lall(followingInd_LR(1:2:35))],0.2,1.4);
A=polyfit(xxx, yyy,1);
[RHO,~]=corr(xxx, yyy);
Alpha1=A(1); % the slope, the first order polynomial coefficient from polyfit (Hurst Componenet >1 ?)
FitValues=polyval(A,xxx);
hold on; 
plot(xxx,FitValues,'-','color',darkgreen,'LineWidth',3);
text(max(xxx)-0.1,max(FitValues)-0.03,sprintf('\\rho=%.2f',RHO),'Color',darkgreen,'FontSize', 15)

% fit the synco data
A=[];Alpha1=[];FitValues=[];RHO=[];xxx=[];yyy=[];
[xxx,yyy]=keepvalues([H_Lall(leadingInd_LR(2:2:36)); H_Rall(followingInd_LR(2:2:36))],...
    [H_Rall(leadingInd_LR(2:2:36)); H_Lall(followingInd_LR(2:2:36))],0.2,1.4);
A=polyfit(xxx, yyy,1);
[RHO,~]=corr(xxx, yyy);
Alpha1=A(1); % the slope, the first order polynomial coefficient from polyfit (Hurst Componenet >1 ?)
FitValues=polyval(A,xxx);
hold on; plot(xxx,FitValues,'-','color',pink,'LineWidth',3);
text(max(xxx)-0.1,max(FitValues)+0.1,sprintf('\\rho=%.2f',RHO),'Color',pink,'FontSize', 15)
% xlim([0.5 1]);ylim([0.5 1]);
xlim([0.2 1.4]);ylim([0.2 1.4]);
plot([0.5 1], [0.5 1],'m--');
plot([0.5 1], [0.5 1],'m--');
hold off;
% legend('synch','','','','synco','','','','','','','','','','','','location','northwest');
% subtitle(['(leading or following: linear regression slope = ' num2str(Alpha1,'%.3f') ')'])
grid on;
end
subplot(1,3,3); %##############
for subplot3 =1;
plot(H_Lall(mutualInd_LR(synchind(1:3))),H_Rall(mutualInd_LR(synchind(1:3))),'.','MarkerSize',30,'color',darkgreen);hold on;
plot(H_Lall(mutualInd_LR(synchind(4:6))),H_Rall(mutualInd_LR(synchind(4:6))),'square','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(mutualInd_LR(synchind(7:9))),H_Rall(mutualInd_LR(synchind(7:9))),'^','MarkerSize',8,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(mutualInd_LR(synchind(10:12))),H_Rall(mutualInd_LR(synchind(10:12))),'pentagram','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(mutualInd_LR(synchind(13:15))),H_Rall(mutualInd_LR(synchind(13:15))),'*','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(mutualInd_LR(synchind(16:18))),H_Rall(mutualInd_LR(synchind(16:18))),'diamond','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(mutualInd_LR(syncoind(1:3))),H_Rall(mutualInd_LR(syncoind(1:3))),'.','MarkerSize',30,'color',pink);
plot(H_Lall(mutualInd_LR(syncoind(4:6))),H_Rall(mutualInd_LR(syncoind(4:6))),'square','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(mutualInd_LR(syncoind(7:9))),H_Rall(mutualInd_LR(syncoind(7:9))),'^','MarkerSize',8,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(mutualInd_LR(syncoind(10:12))),H_Rall(mutualInd_LR(syncoind(10:12))),'pentagram','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(mutualInd_LR(syncoind(13:15))),H_Rall(mutualInd_LR(syncoind(13:15))),'*','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(mutualInd_LR(syncoind(16:18))),H_Rall(mutualInd_LR(syncoind(16:18))),'diamond','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
ax = gca;
ax.FontSize = 13; 
xlabel('DFA Hurst Exponent, Participant A','FontSize',15);
ylabel('DFA Husrt Exponent, Participant B','FontSize',15);
title('Bidirectional','FontSize',15);
% xlim([0.5 1]);ylim([0.5 1]);
xlim([0.2 1.4]);ylim([0.2 1.4]);
% fit all the data
A=[];Alpha1=[];FitValues=[];RHO=[];xxx=[];yyy=[];
[xxx,yyy]=keepvalues(H_Lall(mutualInd_LR(1:36)),...
    H_Rall(mutualInd_LR(1:36)),0.2,1.4);
A=polyfit(xxx,yyy,1);
[RHO,~]=corr(xxx,yyy);
Alpha1=A(1); 
FitValues=polyval(A,xxx);
hold on; plot(xxx,FitValues,'k-','LineWidth',3);
text(max(xxx)-0.2,max(FitValues)-0.2,sprintf('\\rho=%.2f',RHO),'Color',[0 0 0],'FontSize', 15)

% fit the synch data
A=[];Alpha1=[];FitValues=[];RHO=[];xxx=[];yyy=[];
[xxx,yyy]=keepvalues(H_Lall(mutualInd_LR(1:2:35)),...
    H_Rall(mutualInd_LR(1:2:35)),0.2,1.4);
A=polyfit(xxx,yyy,1);
[RHO,~]=corr(xxx,yyy);
Alpha1=A(1); % the slope, the first order polynomial coefficient from polyfit (Hurst Componenet >1 ?)
FitValues=polyval(A,xxx);
hold on; plot(xxx,FitValues,'-','color',darkgreen,'LineWidth',3);
text(max(xxx)-0.1,max(FitValues)-0.1,sprintf('\\rho=%.2f',RHO),'Color',darkgreen,'FontSize', 15)

% fit the synco data
A=[];Alpha1=[];FitValues=[];RHO=[];xxx=[];yyy=[];
[xxx,yyy]=keepvalues(H_Lall(mutualInd_LR(2:2:36)),...
    H_Rall(mutualInd_LR(2:2:36)),0.2,1.4);
A=polyfit(xxx,yyy,1);
[RHO,~]=corr(xxx,yyy);
Alpha1=A(1); % the slope, the first order polynomial coefficient from polyfit (Hurst Componenet >1 ?)
FitValues=polyval(A,xxx);
hold on; plot(xxx,FitValues,'-','color',pink,'LineWidth',3);
text(max(xxx)-0.5,max(FitValues)-0.25,sprintf('\\rho=%.2f',RHO),'Color',pink,'FontSize', 15)
% xlim([0.5 1]);ylim([0.5 1]);
xlim([0.2 1.4]);ylim([0.2 1.4]);
plot([0.5 1], [0.5 1],'m--');
hold off;
% subtitle(['(mutual): linear regression slope = ' num2str(Alpha1,'%.3f') ')'])
grid on;
end
lg=legend('synch - subj pair 1','synch - subj pair 2','synch - subj pair 3','synch - subj pair 4', ...
    'synch - subj pair 5','synch - subj pair 6', ...
    'synco - subj pair 1','synco - subj pair 2','synco - subj pair 3','synco - subj pair 4',...
    'synco - subj pair 5','synco - subj pair 6','location','eastoutside');
lg.Position = [0.9475 0.15 0.01 0.75];
% sgtitle(['H matching as in Marmelat 2012  ^{* PLOT 10-2}']);
delete(findall(gcf,'type','annotation'))
sg=annotation('textbox',[0.07 0.01 0.22 0.05],'string',...
    'H matching as in Marmelat 2012  ^{* PLOT 10-2}')
sg.Rotation=90
set(gcf,'color','w'); % set background white for copying in ubuntu
end


%% PLOT 10-3 EEG matching for theta only in mutual condition
% select 6 channels  ---  Left: 5-F3,9-FC5,14-T7; Right: 7-F4,12-FC6,18-T6
chan6=[5 7 9 12 14 18];
chan6names={'F3','F4','FC5','FC6','T7','T6'};

% select the EEG data
theta_LR_chan; % from PLOT13
theta_LR_chan_copy=theta_LR_chan; % make a copy for testing
theta_LR_chan=theta_LR_chan_copy; % restore;
theta_LR_chan=testing_EEG_LR;

%********************** 6 subplots separating synch/o
cmin=0;cmax=7;
% cmin=-4.1232;cmax=-4.0101; % testing synco
% cmin=4.0101;cmax=4.1232; % testing synch
% cmin=-4.1232;cmax=4.1232; % testing both
for subplots6 =1;
canvas(0.25,0.7);
for i=1:6
subplot(3,2,i); %##############
theta_L_chan=theta_LR_chan(1:144,chan6(i));
theta_R_chan=theta_LR_chan(145:288,chan6(i));
plot(theta_L_chan(mutualInd_LR(synchind(1:3))),theta_R_chan(mutualInd_LR(synchind(1:3))),'.','MarkerSize',30,'color',darkgreen);hold on;
plot(theta_L_chan(mutualInd_LR(synchind(4:6))),theta_R_chan(mutualInd_LR(synchind(4:6))),'square','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(theta_L_chan(mutualInd_LR(synchind(7:9))),theta_R_chan(mutualInd_LR(synchind(7:9))),'^','MarkerSize',8,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(theta_L_chan(mutualInd_LR(synchind(10:12))),theta_R_chan(mutualInd_LR(synchind(10:12))),'pentagram','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(theta_L_chan(mutualInd_LR(synchind(13:15))),theta_R_chan(mutualInd_LR(synchind(13:15))),'*','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(theta_L_chan(mutualInd_LR(synchind(16:18))),theta_R_chan(mutualInd_LR(synchind(16:18))),'diamond','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(theta_L_chan(mutualInd_LR(syncoind(1:3))),theta_R_chan(mutualInd_LR(syncoind(1:3))),'.','MarkerSize',30,'color',pink);
plot(theta_L_chan(mutualInd_LR(syncoind(4:6))),theta_R_chan(mutualInd_LR(syncoind(4:6))),'square','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(theta_L_chan(mutualInd_LR(syncoind(7:9))),theta_R_chan(mutualInd_LR(syncoind(7:9))),'^','MarkerSize',8,'color',pink,'MarkerFaceColor',pink);
plot(theta_L_chan(mutualInd_LR(syncoind(10:12))),theta_R_chan(mutualInd_LR(syncoind(10:12))),'pentagram','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(theta_L_chan(mutualInd_LR(syncoind(13:15))),theta_R_chan(mutualInd_LR(syncoind(13:15))),'*','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(theta_L_chan(mutualInd_LR(syncoind(16:18))),theta_R_chan(mutualInd_LR(syncoind(16:18))),'diamond','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
xlabel('Theta power, Participant L');ylabel('Theta power, Participant R');
title('bidirectional'); 
% fit the synch data
A=[];Alpha1=[];FitValues=[];RHO=[];
A=polyfit(theta_L_chan(mutualInd_LR(1:2:35)),theta_R_chan(mutualInd_LR(1:2:35)),1);
[RHO,~]=corr(theta_L_chan(mutualInd_LR(1:2:35)),theta_R_chan(mutualInd_LR(1:2:35)));
Alpha1=A(1); % the slope, the first order polynomial coefficient from polyfit (Hurst Componenet >1 ?)
FitValues=polyval(A,theta_L_chan(mutualInd_LR(1:2:35)));
hold on; plot(theta_L_chan(mutualInd_LR(1:2:35)),FitValues,'-','color',darkgreen,'LineWidth',3);
text(max(theta_L_chan(mutualInd_LR(1:2:35))),max(FitValues),sprintf('\\rho=%.2f',RHO),'Color',[0 0 0])
% fit the synco data
A=[];Alpha1=[];FitValues=[];RHO=[];
A=polyfit(theta_L_chan(mutualInd_LR(2:2:36)),theta_R_chan(mutualInd_LR(2:2:36)),1);
[RHO,~]=corr(theta_L_chan(mutualInd_LR(2:2:36)),theta_R_chan(mutualInd_LR(2:2:36)));
Alpha1=A(1); % the slope, the first order polynomial coefficient from polyfit (Hurst Componenet >1 ?)
FitValues=polyval(A,theta_L_chan(mutualInd_LR(2:2:36)));
hold on; plot(theta_L_chan(mutualInd_LR(2:2:36)),FitValues,'-','color',pink,'LineWidth',3);
text(max(theta_L_chan(mutualInd_LR(2:2:36))),max(FitValues),sprintf('\\rho=%.2f',RHO),'Color',[0 0 0])
xlim([cmin cmax]);ylim([cmin cmax]);
plot([cmin cmax], [cmin cmax],'m--');
hold off;
% subtitle(['(mutual): linear regression slope = ' num2str(Alpha1,'%.3f') ')'])
grid on;
title(chan6names{i});
end;
lg=legend('synch - subj pair 1','synch - subj pair 2','synch - subj pair 3','synch - subj pair 4', ...
    'synch - subj pair 5','synch - subj pair 6', ...
    'synco - subj pair 1','synco - subj pair 2','synco - subj pair 3','synco - subj pair 4',...
    'synco - subj pair 5','synco - subj pair 6','location','eastoutside');
lg.Position = [0.9475 0.15 0.01 0.025];
% sgtitle(['H matching as in Marmelat 2012  ^{* PLOT 10-2}']);
delete(findall(gcf,'type','annotation'))
% sg=annotation('textbox',[0.07 0.01 0.22 0.05],'string',...
%     'Sum of theta power for chanel 7-F4  ^{* PLOT 10-2}')
% sg=annotation('textbox',[0.07 0.01 0.22 0.05],'string',...
%     'Sum of theta power for chanel 14-T7  ^{* PLOT 10-2}')
sg=annotation('textbox',[0.07 0.01 0.5 0.05],'string',...
    'Sum of theta power for 6 channels in bidirectional state ^{* PLOT 10-2}')
sg.Rotation=90
set(gcf,'color','w'); % set background white for copying in ubuntu
end

% single regression for all 6 chan combining synch/o
cmin=0;cmax=7;
for subplots6 =1;
canvas(0.25,0.75);
for i=1:2
subplot(1,2,i); %##############
theta_L_chan=theta_LR_chan(1:144,chan6(i));
theta_R_chan=theta_LR_chan(145:288,chan6(i));
plot(theta_L_chan(mutualInd_LR(synchind(1:3))),theta_R_chan(mutualInd_LR(synchind(1:3))),'.','MarkerSize',30,'color',darkgreen);hold on;
plot(theta_L_chan(mutualInd_LR(synchind(4:6))),theta_R_chan(mutualInd_LR(synchind(4:6))),'square','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(theta_L_chan(mutualInd_LR(synchind(7:9))),theta_R_chan(mutualInd_LR(synchind(7:9))),'^','MarkerSize',8,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(theta_L_chan(mutualInd_LR(synchind(10:12))),theta_R_chan(mutualInd_LR(synchind(10:12))),'pentagram','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(theta_L_chan(mutualInd_LR(synchind(13:15))),theta_R_chan(mutualInd_LR(synchind(13:15))),'*','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(theta_L_chan(mutualInd_LR(synchind(16:18))),theta_R_chan(mutualInd_LR(synchind(16:18))),'diamond','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(theta_L_chan(mutualInd_LR(syncoind(1:3))),theta_R_chan(mutualInd_LR(syncoind(1:3))),'.','MarkerSize',30,'color',pink);
plot(theta_L_chan(mutualInd_LR(syncoind(4:6))),theta_R_chan(mutualInd_LR(syncoind(4:6))),'square','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(theta_L_chan(mutualInd_LR(syncoind(7:9))),theta_R_chan(mutualInd_LR(syncoind(7:9))),'^','MarkerSize',8,'color',pink,'MarkerFaceColor',pink);
plot(theta_L_chan(mutualInd_LR(syncoind(10:12))),theta_R_chan(mutualInd_LR(syncoind(10:12))),'pentagram','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(theta_L_chan(mutualInd_LR(syncoind(13:15))),theta_R_chan(mutualInd_LR(syncoind(13:15))),'*','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(theta_L_chan(mutualInd_LR(syncoind(16:18))),theta_R_chan(mutualInd_LR(syncoind(16:18))),'diamond','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
xlabel('Theta power, Participant L','FontSize',15);
ylabel('Theta power, Participant R','FontSize',15);
title('bidirectional'); 
% fit the synch/o combined data
A=[];Alpha1=[];FitValues=[];RHO=[];
A=polyfit(theta_L_chan(mutualInd_LR(1:36)),theta_R_chan(mutualInd_LR(1:36)),1);
[RHO,~]=corr(theta_L_chan(mutualInd_LR(1:36)),theta_R_chan(mutualInd_LR(1:36)));
Alpha1=A(1); % the slope, the first order polynomial coefficient from polyfit (Hurst Componenet >1 ?)
FitValues=polyval(A,theta_L_chan(mutualInd_LR(1:36)));
hold on; plot(theta_L_chan(mutualInd_LR(1:36)),FitValues,'-','color',[0 0 0],'LineWidth',3);
text(max(theta_L_chan(mutualInd_LR(1:36))),max(FitValues),sprintf('\\rho=%.2f',RHO),'Color',[0 0 0],'FontSize', 15)
xlim([cmin cmax]);ylim([cmin cmax]);
plot([cmin cmax], [cmin cmax],'m--');
hold off;
% subtitle(['(mutual): linear regression slope = ' num2str(Alpha1,'%.3f') ')'])
grid on;
title(['Channel ' chan6names{i}] ,'FontSize', 14)
end;
lg=legend('synch - subj pair 1','synch - subj pair 2','synch - subj pair 3','synch - subj pair 4', ...
    'synch - subj pair 5','synch - subj pair 6', ...
    'synco - subj pair 1','synco - subj pair 2','synco - subj pair 3','synco - subj pair 4',...
    'synco - subj pair 5','synco - subj pair 6','location','eastoutside');
lg.Position = [0.9475 0.15 0.01 0.025];
% sgtitle(['H matching as in Marmelat 2012  ^{* PLOT 10-2}']);
delete(findall(gcf,'type','annotation'))
% sg=annotation('textbox',[0.07 0.01 0.22 0.05],'string',...
%     'Sum of theta power for chanel 7-F4  ^{* PLOT 10-2}')
% sg=annotation('textbox',[0.07 0.01 0.22 0.05],'string',...
%     'Sum of theta power for chanel 14-T7  ^{* PLOT 10-2}')
sg=annotation('textbox',[0.07 0.01 0.5 0.05],'string',...
    'Sum of theta power for 6 channels in bidirectional state ^{* PLOT 10-2}')
sg.Rotation=90
set(gcf,'color','w'); % set background white for copying in ubuntu
end

% all six channels in one plot
theta_L_chan=mean(theta_LR_chan(1:144,[5 7 9 12 14 18]),2);
theta_R_chan=mean(theta_LR_chan(145:288,[5 7 9 12 14 18]),2);
%********************** 6 subplots in 1
cmin=0;cmax=7;
for subplots6to1 =1;
canvas(0.25,0.25);
plot(theta_L_chan(mutualInd_LR(synchind(1:3))),theta_R_chan(mutualInd_LR(synchind(1:3))),'.','MarkerSize',30,'color',darkgreen);hold on;
plot(theta_L_chan(mutualInd_LR(synchind(4:6))),theta_R_chan(mutualInd_LR(synchind(4:6))),'square','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(theta_L_chan(mutualInd_LR(synchind(7:9))),theta_R_chan(mutualInd_LR(synchind(7:9))),'^','MarkerSize',8,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(theta_L_chan(mutualInd_LR(synchind(10:12))),theta_R_chan(mutualInd_LR(synchind(10:12))),'pentagram','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(theta_L_chan(mutualInd_LR(synchind(13:15))),theta_R_chan(mutualInd_LR(synchind(13:15))),'*','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(theta_L_chan(mutualInd_LR(synchind(16:18))),theta_R_chan(mutualInd_LR(synchind(16:18))),'diamond','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(theta_L_chan(mutualInd_LR(syncoind(1:3))),theta_R_chan(mutualInd_LR(syncoind(1:3))),'.','MarkerSize',30,'color',pink);
plot(theta_L_chan(mutualInd_LR(syncoind(4:6))),theta_R_chan(mutualInd_LR(syncoind(4:6))),'square','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(theta_L_chan(mutualInd_LR(syncoind(7:9))),theta_R_chan(mutualInd_LR(syncoind(7:9))),'^','MarkerSize',8,'color',pink,'MarkerFaceColor',pink);
plot(theta_L_chan(mutualInd_LR(syncoind(10:12))),theta_R_chan(mutualInd_LR(syncoind(10:12))),'pentagram','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(theta_L_chan(mutualInd_LR(syncoind(13:15))),theta_R_chan(mutualInd_LR(syncoind(13:15))),'*','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(theta_L_chan(mutualInd_LR(syncoind(16:18))),theta_R_chan(mutualInd_LR(syncoind(16:18))),'diamond','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
xlabel('Theta power, Participant L','FontSize',15);ylabel('Theta power, Participant R','FontSize',15);
title('Sum of theta power from F3,F4,FC5,FC6,T7,T6 in bidrectional state'); 
% fit the synch data
A=[];Alpha1=[];FitValues=[];RHO=[]
A=polyfit(theta_L_chan(mutualInd_LR(1:2:35)),theta_R_chan(mutualInd_LR(1:2:35)),1);
[RHO,~]=corr(theta_L_chan(mutualInd_LR(1:2:35)),theta_R_chan(mutualInd_LR(1:2:35)));
Alpha1=A(1); % the slope, the first order polynomial coefficient from polyfit (Hurst Componenet >1 ?)
FitValues=polyval(A,theta_L_chan(mutualInd_LR(1:2:35)));
hold on; plot(theta_L_chan(mutualInd_LR(1:2:35)),FitValues,'-','color',darkgreen,'LineWidth',3);
text(max(theta_L_chan(mutualInd_LR(1:2:35))),max(FitValues),sprintf('\\rho=%.2f',RHO),'Color',[0 0 0])
% fit the synco data
A=[];Alpha1=[];FitValues=[];RHO=[]
A=polyfit(theta_L_chan(mutualInd_LR(2:2:36)),theta_R_chan(mutualInd_LR(2:2:36)),1);
[RHO,~]=corr(theta_L_chan(mutualInd_LR(2:2:36)),theta_R_chan(mutualInd_LR(2:2:36)));
Alpha1=A(1); % the slope, the first order polynomial coefficient from polyfit (Hurst Componenet >1 ?)
FitValues=polyval(A,theta_L_chan(mutualInd_LR(2:2:36)));
hold on; plot(theta_L_chan(mutualInd_LR(2:2:36)),FitValues,'-','color',pink,'LineWidth',3);
text(max(theta_L_chan(mutualInd_LR(2:2:36))),max(FitValues),sprintf('\\rho=%.2f',RHO),'Color',[0 0 0])
xlim([cmin cmax]);ylim([cmin cmax]);
plot([cmin cmax], [cmin cmax],'m--');
hold off;
% subtitle(['(mutual): linear regression slope = ' num2str(Alpha1,'%.3f') ')'])
grid on;
lg=legend('synch - subj pair 1','synch - subj pair 2','synch - subj pair 3','synch - subj pair 4', ...
    'synch - subj pair 5','synch - subj pair 6', ...
    'synco - subj pair 1','synco - subj pair 2','synco - subj pair 3','synco - subj pair 4',...
    'synco - subj pair 5','synco - subj pair 6','location','eastoutside');
lg.Position = [0.9475 0.15 0.01 0.025];
% sgtitle(['H matching as in Marmelat 2012  ^{* PLOT 10-2}']);
delete(findall(gcf,'type','annotation'))
sg=annotation('textbox',[0.07 0.01 0.5 0.05],'string',...
    'Sum of theta power for 6 channels in bidirectional state ^{* PLOT 10-3}')
sg.Rotation=90
set(gcf,'color','w'); % set background white for copying in ubuntu
end
%% PLOT 11 Xcorr(0) with slinding win of 20 intervals (matched int) in one subj
% clear
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);
sessions={'synch','synco','synch','synco','synch','synco','synch','synco','synch','synco',...
    'synch','synco'};
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
sgtitle(['session ' num2str(seed) ' ' sessionTypes{session} ' BPint corr xcorr diff ^{* PLOT 11}']);
%% PLOT 11-1 Xcorr(0) with slinding win of 20 intervals (matched int) in all subj
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);
sessions={'synch','synco','synch','synco','synch','synco','synch','synco','synch','synco',...
    'synch','synco'};
% sliding win
BPwin=20;
% organize conditions in all sessions into a matrix of 8 x 12
condition_all_mat=NaN(numSes,12); % (Refer to SECT 10-1)
% Compute Xcorr(0) for all trials
intNum=NaN(numSes,12);BPint_xcorrSeries=cell(numSes,12);
for s=1:numSes;
    clear intervals conditions
    runid=num2str(seeds(s,:));
    load(['/ssd/zhibin/1overf/' runid '_2P/Cleaned_data/clean_' runid '.mat' ]);
    condition_all_mat(s,:)=conditions;
    intL_good_dmean=[];intR_good_dmean=[];
    for b=1:12
        % remove the mean (no need for DFA)
        intL_good_dmean=intervals{b}(:,1)-mean(intervals{b}(:,1));
        intR_good_dmean=intervals{b}(:,2)-mean(intervals{b}(:,2));
        % sliding window of 20 to compute xcorr at lag 0
        intNum(s,b)=size(intervals{b},1);
        for k=1:(intNum(s,b)-BPwin+1)
            % xcorr
            r12=[];lags12=[];
            [r12,lags12]=xcorr(intL_good_dmean(k:k+BPwin-1),intR_good_dmean(k:k+BPwin-1),10,'normalized');
            BPint_xcorrSeries{s,b}(k)=r12(11);
        end
    end
end
% find indicies for all mutual conditions (Refer to SECT 10-1)
find(condition_all_mat==4)
% plot xcorr(0) for synch/synco (4 subplots for 4 pairs of subj)
canvas(0.4,0.4)
for p=1:6
    subplot(2,3,p)
        % synch
        inds=[];inds=find(condition_all_mat((2*p-1),:)==4);
        for t=1:3
            plot(BPint_xcorrSeries{(2*p-1),inds(t)},'color',darkgreen);hold on;
        end
        % synco
        inds=[];inds=find(condition_all_mat((2*p),:)==4);
        for t=1:3
            plot(BPint_xcorrSeries{(2*p),inds(t)},'color',pink);hold on;
        end
        hold off;
        ylabel('xcorr coef'); xlabel(['taps']);
        ylim([-1 1]); yline(0,'color',deepyellow);
        legend({'synch','','','synco','',''},'location','southeast');
        title(['pair of subjects - ' num2str(p)]);
end
sgtitle({['Mutual: xcorr(0) smoo win ' num2str(BPwin) ' ^{* PLOT 11-1}']},'Color',condicolors(4,:));
%% SECT 11-2 EEG (-500ms)-> Xcorr(0)/(-1/+1)(matched int) -Corr -PLS in all subj
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);
sessions={'synch','synco','synch','synco','synch','synco','synch','synco','synch','synco',...
    'synch','synco'};
% organize conditions in all sessions into a matrix of 8 x 12
condition_all_mat=NaN(numSes,12); % (Refer to SECT 10-1)
% Compute Xcorr(0)/Xcorr(-1/+1) for all trials
BPint_xcorrSeries=nan(2,numSes,12);
BPint_xcorrSeries_synch=nan(2,numSes,12);
BPint_xcorrSeries_synco=nan(2,numSes,12);
tic
for s=1:numSes;
    clear intervals conditions
    runid=num2str(seeds(s,:));
    load(['/ssd/zhibin/1overf/' runid '_2P/Cleaned_data/clean_' runid '.mat' ]);
    condition_all_mat(s,:)=conditions;
    intL_good_dmean=[];intR_good_dmean=[];
    for b=1:12
        % remove the mean (no need for DFA)
        intL_good_dmean=intervals{b}(:,1)-mean(intervals{b}(:,1));
        intR_good_dmean=intervals{b}(:,2)-mean(intervals{b}(:,2));
        [r12,lags12]=xcorr(intL_good_dmean,intR_good_dmean,10,'normalized');
        % for synch
        BPint_xcorrSeries_synch(1,s,b)=r12(10);% xcorr(-1)
        BPint_xcorrSeries_synch(2,s,b)=r12(12);% xcorr(+1)
        % for synco
        BPint_xcorrSeries_synco(1,s,b)=r12(11);% xcorr(0)
        BPint_xcorrSeries_synco(2,s,b)=r12(11);% xcorr(0)
        % for both syn types
        % BPint_xcorrSeries(1,s,b)=r12(11);% xcorr(0)
        % BPint_xcorrSeries(2,s,b)=r12(11);% xcorr(0)
    end
end
toc % 77 sec

% make a copy of the results (skip)
% BPint_xcorrSeries_synch=BPint_xcorrSeries;
% BPint_xcorrSeries_synco=BPint_xcorrSeries;

% select the result to analyze
% BPint_xcorrSeries=BPint_xcorrSeries_synch;
% BPint_xcorrSeries=BPint_xcorrSeries_synco;

% Organize BPint_xcorrSeries for PLS
% same way as we did with H in PLOT 13
BPint_xcorrSeries; % (2xnumSesx12) for all sessions (matched int) (original order)
BPint_xcorrSeries_L=squeeze(BPint_xcorrSeries(1,:,:));
BPint_xcorrSeries_R=squeeze(BPint_xcorrSeries(2,:,:));
% squeeze into 1 vector from the 96 blocks for each subject, for corr with pow in each chan
BPint_xcorrSeries_L=reshape(BPint_xcorrSeries_L',[],1);% 144x1 (each element from one block in time sequence) 
BPint_xcorrSeries_R=reshape(BPint_xcorrSeries_R',[],1);
% Combine L and R
BPint_xcorrSeries_LR=[BPint_xcorrSeries_L;BPint_xcorrSeries_R]; % to be used for corr and PLS

% indices for synch or synco
Inds_synch=[1:12 25:36 49:60 73:84 97:108 121:132]';
Inds_synch_LR=[Inds_synch;12*numSes+Inds_synch];
Inds_synco=[13:24 37:48 61:72 85:96 109:120 133:144]';
Inds_synco_LR=[Inds_synco;12*numSes+Inds_synco];
% ******
% select_ind_LR=select_ind_LR;
% select_ind_LR=Inds_synch_LR;
% select_ind_LR=Inds_synco_LR;
select_ind_LR; % *******

% indicies for 4 states
Inds4_LR; % 1:288 for all sessions % from SECT 10-1 
% indices for synch or synco
Inds_synch4=Inds4_LR([1:3 7:9 13:15 19:21 25:27 31:33],:); % synchind
Inds_synch4_LR=[Inds_synch4;12*numSes+Inds_synch4];
Inds_synco4=Inds4_LR([4:6 10:12 16:18 22:24 28:30 34:36],:); % syncoind
Inds_synco4_LR=[Inds_synco4;12*numSes+Inds_synco4];
% *******
% select_Inds4_LR=Inds4_LR; 
% select_Inds4_LR=Inds_synch4_LR; 
% select_Inds4_LR=Inds_synco4_LR; 
select_Inds4_LR; % *******

% run SECT 12 and PLOT 13 to organize EEG first (-500ms)
    % Corr
    % Compute the correlation between sum-EEG pow and BPint_xcorrSeries_LR
    % for c=1:32
    %     delta_L_Xcorr_L_corr(c)=corr(delta_L_chan(:,c),BPint_xcorrSeries_L);
    %     theta_L_Xcorr_L_corr(c)=corr(theta_L_chan(:,c),BPint_xcorrSeries_L);
    %     alpha_L_Xcorr_L_corr(c)=corr(alpha_L_chan(:,c),BPint_xcorrSeries_L);
    %     beta_L_Xcorr_L_corr(c)=corr(beta_L_chan(:,c),BPint_xcorrSeries_L);
    %     gamma_L_Xcorr_L_corr(c)=corr(gamma_L_chan(:,c),BPint_xcorrSeries_L);
    %     delta_R_Xcorr_R_corr(c)=corr(delta_R_chan(:,c),BPint_xcorrSeries_R);
    %     theta_R_Xcorr_R_corr(c)=corr(theta_R_chan(:,c),BPint_xcorrSeries_R);
    %     alpha_R_Xcorr_R_corr(c)=corr(alpha_R_chan(:,c),BPint_xcorrSeries_R);
    %     beta_R_Xcorr_R_corr(c)=corr(beta_R_chan(:,c),BPint_xcorrSeries_R);
    %     gamma_R_Xcorr_R_corr(c)=corr(gamma_R_chan(:,c),BPint_xcorrSeries_R);
    % end
% indices for 4 states from PLOT 10-1
% fix the scale in the data (as in PLOT 16)
% delta_LR_chan = delta_LR_chan./(ones(288,1)*std(delta_LR_chan)); %organized in SECT 13
% theta_LR_chan = theta_LR_chan./(ones(288,1)*std(theta_LR_chan));
% alpha_LR_chan = alpha_LR_chan./(ones(288,1)*std(alpha_LR_chan));
% beta_LR_chan = beta_LR_chan./(ones(288,1)*std(beta_LR_chan));
% gamma_LR_chan = gamma_LR_chan./(ones(288,1)*std(gamma_LR_chan));
delta_LR_chan;theta_LR_chan;alpha_LR_chan;beta_LR_chan;gamma_LR_chan; %organized from PLOT 13

% Correlation combining all states (all states x 5 freq)
delta_LR_Xcorr_LR_corr=zeros(32,1);
theta_LR_Xcorr_LR_corr=zeros(32,1);
alpha_LR_Xcorr_LR_corr=zeros(32,1);
beta_LR_Xcorr_LR_corr=zeros(32,1);
gamma_LR_Xcorr_LR_corr=zeros(32,1);
for c=1:32
    delta_LR_Xcorr_LR_corr(c)=corr(delta_LR_chan(select_ind_LR,c),BPint_xcorrSeries_LR(select_ind_LR));
    theta_LR_Xcorr_LR_corr(c)=corr(theta_LR_chan(select_ind_LR,c),BPint_xcorrSeries_LR(select_ind_LR));
    alpha_LR_Xcorr_LR_corr(c)=corr(alpha_LR_chan(select_ind_LR,c),BPint_xcorrSeries_LR(select_ind_LR));
    beta_LR_Xcorr_LR_corr(c)=corr(beta_LR_chan(select_ind_LR,c),BPint_xcorrSeries_LR(select_ind_LR));
    gamma_LR_Xcorr_LR_corr(c)=corr(gamma_LR_chan(select_ind_LR,c),BPint_xcorrSeries_LR(select_ind_LR));
end
freq5_LR_Xcorr_LR_corr=[delta_LR_Xcorr_LR_corr theta_LR_Xcorr_LR_corr alpha_LR_Xcorr_LR_corr ...
    beta_LR_Xcorr_LR_corr gamma_LR_Xcorr_LR_corr];

% Plot the 5 freq (1 x 5 freq)
for plot_1by5=1;
canvas(0.3,0.2)
cmin=-0.6;cmax=0.6;
    for i=1:5
    subplot(1,5,i)
    selected_LR_Xcorr_LR_4corr=[];
    selected_LR_Xcorr_LR_4corr=freq5_LR_Xcorr_LR_corr(:,i);
    topoplot(selected_LR_Xcorr_LR_4corr,channels,'nosedir','+X','style','map');
    clim([cmin cmax]);
    end
colormap(hnc)
cb=colorbar;
cb.AxisLocation = 'out';
cb.Position = [0.92 0.15 0.01 0.75];
set(gcf,'color','w'); % set background white for copying in ubuntu
if true %annotations
delete(findall(gcf,'type','annotation'))
h0=annotation('textbox',[0.17 0.95 0.05 0.03],'string','Delta','color',[0 0 0])
h1=annotation('textbox',[0.33 0.95 0.05 0.03],'string','Theta','color',[0 0 0])
h2=annotation('textbox',[0.5 0.95 0.05 0.03],'string','Alpha','color',[0 0 0])
h3=annotation('textbox',[0.66 0.95 0.05 0.03],'string','Beta','color',[0 0 0])
h4=annotation('textbox',[0.82 0.95 0.05 0.03],'string','Gamma','color',[0 0 0])
v0=annotation('textbox',[0.1 0.4 0.05 0.03],'string','All conditions');
set(v0,'Rotation',90);
set([h0 h1 h2 h3 h4 v0], 'fitboxtotext','on',...
    'edgecolor','none')
end
%{
sg=annotation('textbox',[0.3 0.01 0.4 0.05],'string',...
    'Correlation of sum-EEG (-500ms) and Xcorr(0) ^{* PLOT 11-3}')
sg=annotation('textbox',[0.3 0.01 0.4 0.15],'string',...
     {['Correlation of sum-EEG (-500ms) and Xcorr(-1/+1) in Synch ^{* PLOT 11-2} '], char(datetime('now'))},'Color',syn2colors(1,:));
sg=annotation('textbox',[0.3 0.01 0.4 0.18],'string',...
    {['Correlation of sum-EEG (-500ms) and Xcorr(0) in Synco ^{* PLOT 11-2} '], char(datetime('now'))},'Color',syn2colors(2,:));
%} 

end

% Correlation in 4 states (4 states x 5 freq)
delta_LR_Xcorr_LR_4corr=zeros(4,32);
theta_LR_Xcorr_LR_4corr=zeros(4,32);
alpha_LR_Xcorr_LR_4corr=zeros(4,32);
beta_LR_Xcorr_LR_4corr=zeros(4,32);
gamma_LR_Xcorr_LR_4corr=zeros(4,32);
for s=1:4
    for c=1:32
        delta_LR_Xcorr_LR_4corr(s,c)=corr(delta_LR_chan(select_Inds4_LR(:,s),c),BPint_xcorrSeries_LR(select_Inds4_LR(:,s)));
        theta_LR_Xcorr_LR_4corr(s,c)=corr(theta_LR_chan(select_Inds4_LR(:,s),c),BPint_xcorrSeries_LR(select_Inds4_LR(:,s)));
        alpha_LR_Xcorr_LR_4corr(s,c)=corr(alpha_LR_chan(select_Inds4_LR(:,s),c),BPint_xcorrSeries_LR(select_Inds4_LR(:,s)));
        beta_LR_Xcorr_LR_4corr(s,c)=corr(beta_LR_chan(select_Inds4_LR(:,s),c),BPint_xcorrSeries_LR(select_Inds4_LR(:,s)));
        gamma_LR_Xcorr_LR_4corr(s,c)=corr(gamma_LR_chan(select_Inds4_LR(:,s),c),BPint_xcorrSeries_LR(select_Inds4_LR(:,s)));
    end
end

% Examine the corrcoef values
%{ 
clf
EEG_test=beta_LR_chan(select_Inds4_LR(:,4),5);
Xcorr_test=BPint_xcorrSeries_LR(select_Inds4_LR(:,4));
plot(Xcorr_test,EEG_test,'b.');
corrcoef(Xcorr_test,EEG_test) 
%}

% Plot 4 states(4x5)
for plot_4by5=1;
canvas(0.3,0.5);
cmin=-0.7;cmax=0.7;
for s=1:4
    subplot(4,5,5*(s-1)+1);
    topoplot(delta_LR_Xcorr_LR_4corr(s,:),channels,'nosedir','+X','style','fill','numcontour',12);
    clim([cmin cmax]);
    subplot(4,5,5*(s-1)+2);
    topoplot(theta_LR_Xcorr_LR_4corr(s,:),channels,'nosedir','+X','style','fill','numcontour',12);
    clim([cmin cmax]);
    subplot(4,5,5*(s-1)+3);
    topoplot(alpha_LR_Xcorr_LR_4corr(s,:),channels,'nosedir','+X','style','fill','numcontour',12);
    clim([cmin cmax]);
    subplot(4,5,5*(s-1)+4);
    topoplot(beta_LR_Xcorr_LR_4corr(s,:),channels,'nosedir','+X','style','fill','numcontour',12);
    clim([cmin cmax]);
    subplot(4,5,5*(s-1)+5);
    topoplot(gamma_LR_Xcorr_LR_4corr(s,:),channels,'nosedir','+X','style','fill','numcontour',12);
    clim([cmin cmax]);
end
colormap(hotncold(12))
% https://www.mathworks.com/help/matlab/ref/matlab.graphics.illustration.colorbar-properties.html
cb=colorbar;
cb.AxisLocation = 'out';
cb.Position = [0.92 0.15 0.01 0.75];
% annotations
delete(findall(gcf,'type','annotation'))
h0=annotation('textbox',[0.17 0.95 0.05 0.03],'string','Delta','color',[0 0 0])
h1=annotation('textbox',[0.33 0.95 0.05 0.03],'string','Theta','color',[0 0 0])
h2=annotation('textbox',[0.5 0.95 0.05 0.03],'string','Alpha','color',[0 0 0])
h3=annotation('textbox',[0.66 0.95 0.05 0.03],'string','Beta','color',[0 0 0])
h4=annotation('textbox',[0.81 0.95 0.05 0.03],'string','Gamma','color',[0 0 0])
v0=annotation('textbox',[0.14 0.15 0.05 0.03],'string','Mutual','color',condicolors(4,:))
v1=annotation('textbox',[0.14 0.37 0.05 0.03],'string','Following','color',condicolors(3,:))
v2=annotation('textbox',[0.14 0.59 0.05 0.03],'string','Leading','color',condicolors(2,:))
v3=annotation('textbox',[0.14 0.81 0.05 0.03],'string','Uncoupled','color',condicolors(1,:))
set(v0,'Rotation',90);set(v1,'Rotation',90);set(v2,'Rotation',90);set(v3,'Rotation',90);
set([h0 h1 h2 h3 h4 v0 v1 v2 v3], 'fitboxtotext','on',...
    'edgecolor','none')
end
delete(sg)
%{
sg=annotation('textbox',[0.3 0.01 0.4 0.05],'string',...
    'Correlation of sum-EEG (-500ms) and Xcorr(0) ^{* PLOT 11-3}')
sg=annotation('textbox',[0.3 0.01 0.45 0.07],'string',...
     {['Correlation of sum-EEG (-500ms) and Xcorr(-1/+1) in Synch ^{* PLOT 11-2} '], char(datetime('now'))},'Color',syn2colors(1,:));
sg=annotation('textbox',[0.3 0.01 0.4 0.07],'string',...
    {['Correlation of sum-EEG (-500ms) and Xcorr(0) in Synco ^{* PLOT 11-2} '], char(datetime('now'))},'Color',syn2colors(2,:));
%} 
set(gcf,'color','w'); % set background white for copying in ubuntu


% PLS
addpath /home/zhibinz2/Documents/GitHub/matlab/ramesh/plsmodel
addpath(genpath('/home/zhibinz2/Documents/GitHub/matlab/external/')); 
addpath /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/PLS
% get pow5forpls3 from PLOT 16
pow5forpls3; % 288x160
% (ALL states: 5freq x 32chan = 160 predictors x 288 trials)
% -updated to use mynpls_pred function
reg=[];ypred_fit=[];X_MC=[];Y_MC=[];
clear plsmodel;
% [reg,ypred_fit,X_MC,Y_MC] = mynpls_pred(pow5forpls3,BPint_xcorrSeries_LR);
[reg,ypred_fit,X_MC,Y_MC] = mynpls_pred(pow5forpls3(select_ind_LR,:),BPint_xcorrSeries_LR(select_ind_LR));
ssOr=sum(sum((Y_MC-mean(Y_MC)).^2))
AIC=nan(1,10);R2s=nan(1,10);ssERs=nan(1,10);
for Fac=1:10;
    R2=[];ssEr=[];
    ssEr=sum(sum((Y_MC-ypred_fit(:,Fac)).^2))
    ssERs(Fac)=ssEr;
    R2=100*(1-ssEr/ssOr)
    R2s(Fac)=R2;
    % Cal AIC
    AIC(Fac)=log(ssEr)+2*(Fac);
end
canvas(0.4, 0.3);
subplot(2,1,1);
plot(1:10,AIC,'.','MarkerSize',20);xlabel('nFac');ylabel('AIC');
xlim([0 11]);
subtitle('All states');
subplot(2,1,2);
plot(1:10,R2s,'r.','MarkerSize',20);xlabel('nFac');ylabel('R^2');
xlim([0 11]);
subtitle('All states');
set(gcf,'color','w'); % set background white for copying in ubuntu
% sgtitle({['AIC & R^2 for PLS model in all statues: sum-EEG(-500ms) -> Xcorr(0) ^{ *PLOT 11-2}'], char(datetime('now'))});
sgtitle({['AIC & R^2 for PLS model in all statues: sum-EEG(-500ms) -> Xcorr(-1/+1) in synch ^{ *PLOT 11-2}'], char(datetime('now'))},'Color',syn2colors(1,:));
% sgtitle({['AIC & R^2 for PLS model in all statues: sum-EEG(-500ms) -> Xcorr(0) in synco ^{ *PLOT 11-2}'], char(datetime('now'))},'Color',syn2colors(2,:));


% Try cross validation (myxvalidation)
% select the Fac with minimal AIC
Fac=1
% Choose #of elements to remove
% nSet=6; % for synch and synco combined
% nSet=3; % for synch or synco only
% do the 8-fold cross validation
ssEr=[];R2=[];
% X=pow5forpls3;Y=BPint_xcorrSeries_LR;
X=pow5forpls3(select_ind_LR,:);Y=BPint_xcorrSeries_LR(select_ind_LR);
nSet=1;
[ssEr,R2]  = myxvalidation(X,Y,Fac,nSet);
% reshape regs{1:4}{Fac} 
clear plsmodel
% reshape from 160 x 1 back to 5 x 32 (freq x chan)
plsmodel.weights=reshape(reg{Fac},5,32); 
canvas(0.25,0.2)
cmin=-0.02;cmax=0.02;
imagesc(plsmodel.weights);colorbar; % caxis([-2 2]*10E-7); % by default, imagesc reverse the Y 
yticks([1:5]);yticklabels({'Delta','Theta','Alpha','Beta','Gamma'});
colormap('jet'); clim([cmin cmax]);
set(gca, 'YDir','normal');
xticks([1:32]);xticklabels([labels]);xtickangle(90);grid on;
set(gcf,'color','w'); % set background white for copying in ubuntu
% title({['PLS model (R2= ' num2str(round(R2,1)) ') in all statues: ' ...
%     'sum-EEG(-500ms) -> Xcorr(0)  ^{* PLOT 11-2}'], char(datetime('now'))});
title({['PLS model (R2= ' num2str(round(R2,1)) ') in all statues: ' ...
    'sum-EEG(-500ms) -> Xcorr(-1/+1) in synch ^{* PLOT 11-2}'], char(datetime('now'))},'Color',syn2colors(1,:));
title({['PLS model (R2= ' num2str(round(R2,1)) ') in all statues: ' ...
    'sum-EEG(-500ms) -> Xcorr(0) in synco ^{* PLOT 11-2}'], char(datetime('now'))},'Color',syn2colors(2,:));

% (Each of the 4 states: 5freq x 32 chan = 160 predictors x 72
% trials)-updated to use mynpls_pred function
cd /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/PLS
regs=[];AIC=nan(2,10);R2s=nan(4,10);ssERs=nan(4,10);
for c=1:4 % four states
    R2=[];reg=[];ypred=[];
    [reg,ypred_fit,X_MC,Y_MC] = mynpls_pred(pow5forpls3(select_Inds4_LR(:,c),:),BPint_xcorrSeries_LR(select_Inds4_LR(:,c)));
    regs{c}=reg;
    ssOr=[];
    ssOr=sum(sum((Y_MC-mean(Y_MC)).^2))
    for Fac=1:10;
        R2=[];ssEr=[];
        ssEr=sum(sum((Y_MC-ypred_fit(:,Fac)).^2))
        ssERs(c,Fac)=ssEr;
        R2=100*(1-ssEr/ssOr)
        R2s(c,Fac)=R2;
        % Cal AIC
        AIC(c,Fac)=log(ssEr)+2*(Fac);
    end
end
canvas(0.4,0.4);
for c=1:4
    subplot(2,4,c);
    plot(1:10,AIC(c,:),'.','MarkerSize',20);xlabel('nFac');ylabel('AIC');
    ylim([0 17]);xlim([0 11]);
    subtitle(states4names{c},'Color',condicolors(c,:));
    subplot(2,4,4+c);
    plot(1:10,R2s(c,:),'r.','MarkerSize',20);xlabel('nFac');ylabel('R^2');
    ylim([10 110]);xlim([0 11]);
    subtitle(states4names{c},'Color',condicolors(c,:));
end
set(gcf,'color','w'); % set background white for copying in ubuntu
% sgtitle({['AIC & R^2 for PLS model in all 4 statues: sum-EEG(-500ms) -> Xcorr(0) ^{ *PLOT 11-2}'], char(datetime('now'))});
sgtitle({['AIC & R^2 for PLS model in all 4 statues: sum-EEG(-500ms) -> Xcorr(-1/+1) in synch^{ *PLOT 11-2}'], char(datetime('now'))});
% sgtitle({['AIC & R^2 for PLS model in all 4 statues: sum-EEG(-500ms) -> Xcorr(0) in synco^{ *PLOT 11-2}'], char(datetime('now'))});
% Try cross validation (myxvalidation)
% select the Fac with minimal AIC
Fac=1
ssErs=[];R2s=[];
% Choose #of elements to remove
% nSet=6; % for synch and synco combined
% nSet=3; % for synch or synco only
% do the cross validation
for c=1:4
    ssEr=[];R2=[];
    X=pow5forpls3(select_Inds4_LR(:,c),:);Y=BPint_xcorrSeries_LR(select_Inds4_LR(:,c));
    nSet=1;
    [ssEr,R2]  = myxvalidation(X,Y,Fac,nSet);
    ssErs(c)=ssEr;
    R2s(c)=R2;
end
% reshape regs{1:4}{Fac} 
clear plsmodel
for c=1:4
    plsmodel(c).weights=reshape(regs{c}{Fac},5,32); 
end
% Plot the 4 states
canvas(0.35,0.4);
cmin=-0.02;cmax=0.02;
for c=1:4
    subplot(2,2,c);
    imagesc(plsmodel(c).weights);colorbar; 
    yticks([1:5]);yticklabels({'Delta','Theta','Alpha','Beta','Gamma'});
    set(gca, 'YDir','normal');
    xticks([1:32]);xticklabels([labels]);xtickangle(90);
    colormap('jet'); clim([cmin cmax]);
    % caxis([-0.1 0.1]);
    title([states4names{c} ': PLS model (R^2= ' num2str(round(R2s(c),1)) '  Fac= ' num2str(Fac) ') '], ...
        'Color',condicolors(c,:));
    grid on;
end
set(gcf,'color','w'); % set background white for copying in ubuntu
colormap(hnc)
% sgtitle({['PLS model: sum-EEG(-500ms) -> Xcorr(0) ^{* PLOT 11-2}'], char(datetime('now'))});
sgtitle({['PLS model: sum-EEG(-500ms) -> Xcorr(-1/+1) in synch ^{* PLOT 11-2}'],char(datetime('now'))});
sgtitle({['PLS model: sum-EEG(-500ms) -> Xcorr(0) in synco ^{* PLOT 11-2}'],char(datetime('now'))});
% topoplot for uncouple and mutual (now in 4 states)
addpath /home/zhibinz2/Documents/GitHub/matlab-archive/hnlcode/common/gen_code/color
hnc = hotncold(100);
band5names={'Delta','Theta','Alpha','Beta','Gamma'};
states2names={'Uncouple','Mutual'};
states4names;
% figure;
canvas(0.3,0.8)
for states4=1
% c=[1 4];
c=[1 2 3 4];
for s=1:4 % 2
    for b=1:5
        subplot(4,5,(s-1)*5+b)
        topoplot(plsmodel(c(s)).weights(b,:),channels,'nosedir','+X','style','map');
        clim([cmin cmax]);
    end
end
colormap(hnc)
% https://www.mathworks.com/help/matlab/ref/matlab.graphics.illustration.colorbar-properties.html
cb=colorbar;
cb.AxisLocation = 'out';
cb.Position = [0.92 0.15 0.01 0.75];
set(gcf,'color','w'); % set background white for copying in ubuntu
delete(findall(gcf,'type','annotation'))
h0=annotation('textbox',[0.17 0.95 0.05 0.03],'string','Delta','color',[0 0 0])
h1=annotation('textbox',[0.33 0.95 0.05 0.03],'string','Theta','color',[0 0 0])
h2=annotation('textbox',[0.5 0.95 0.05 0.03],'string','Alpha','color',[0 0 0])
h3=annotation('textbox',[0.66 0.95 0.05 0.03],'string','Beta','color',[0 0 0])
h4=annotation('textbox',[0.82 0.95 0.05 0.03],'string','Gamma','color',[0 0 0])
v0=annotation('textbox',[0.1 0.15 0.05 0.03],'string',...
    {[states4names{c(4)}], [' (R^2= ' num2str(round(R2s(c(4)),1)) ...
            '  Fac= ' num2str(Fac) ') ']}...
            ,'color',condicolors(4,:));
set(v0,'Rotation',90);
v1=annotation('textbox',[0.1 0.35 0.05 0.03],'string',...
    {[states4names{c(3)}], [' (R^2= ' num2str(round(R2s(c(3)),1)) ...
            '  Fac= ' num2str(Fac) ') ']} ...
            ,'color',condicolors(3,:));
set(v1,'Rotation',90);
v2=annotation('textbox',[0.1 0.56 0.05 0.03],'string',...
    {[states4names{c(2)}], [' (R^2= ' num2str(round(R2s(c(2)),1)) ...
            '  Fac= ' num2str(Fac) ') ']} ...
            ,'color',condicolors(2,:));
set(v2,'Rotation',90);
v3=annotation('textbox',[0.1 0.79 0.05 0.03],'string',...
    {[states4names{c(1)}], [' (R^2= ' num2str(round(R2s(c(1)),1)) ...
            '  Fac= ' num2str(Fac) ') ']} ...
            ,'color',condicolors(1,:));
set(v3,'Rotation',90);
end
delete(sg)
% sg=annotation('textbox',[0.3 0.01 0.4 0.05],'string',...
%     {['PLS model: sum-EEG(-500ms) -> Xcorr(0) ^{* PLOT 11-2}'],char(datetime('now'))});
sg=annotation('textbox',[0.3 0.03 0.4 0.06],'string',...
    {['PLS model: sum-EEG(-500ms) -> Xcorr(-1/+1) in synch ^{* PLOT 11-2}'],char(datetime('now'))});
sg=annotation('textbox',[0.3 0.03 0.4 0.06],'string',...
    {['PLS model: sum-EEG(-500ms) -> Xcorr(0) in synco^{* PLOT 11-2}'],char(datetime('now'))});
set([h0 h1 h2 h3 h4 v0 v1 v2 v3], 'fitboxtotext','on',...
    'edgecolor','none')

% topoplot for following state only in synco
canvas(0.3,0.2)
for states1=1;
cmin=-0.007;cmax=0.007;
c=[3];
for s=1    
    for b=1:5
        subplot(1,5,(s-1)*5+b)
        topoplot(plsmodel(c(s)).weights(b,:),channels,'nosedir','+X','style','map');
        clim([cmin cmax]);
    end
end
colormap(hnc)
cb=colorbar;
cb.AxisLocation = 'out';
cb.Position = [0.92 0.15 0.01 0.75];
set(gcf,'color','w'); % set background white for copying in ubuntu
delete(findall(gcf,'type','annotation'))
h0=annotation('textbox',[0.17 0.95 0.05 0.03],'string','Delta','color',[0 0 0])
h1=annotation('textbox',[0.33 0.95 0.05 0.03],'string','Theta','color',[0 0 0])
h2=annotation('textbox',[0.5 0.95 0.05 0.03],'string','Alpha','color',[0 0 0])
h3=annotation('textbox',[0.66 0.95 0.05 0.03],'string','Beta','color',[0 0 0])
h4=annotation('textbox',[0.82 0.95 0.05 0.03],'string','Gamma','color',[0 0 0])
v0=annotation('textbox',[0.1 0.4 0.05 0.03],'string',...
    {[states4names{c(1)}], [' (R^2= ' num2str(round(R2s(c(1)),1)) ...
            '  Fac= ' num2str(Fac) ') ']}...
            ,'color',condicolors(3,:));
set(v0,'Rotation',90);
end
set([h0 h1 h2 h3 h4 v0 v1 v2 v3], 'fitboxtotext','on',...
    'edgecolor','none')
sg=annotation('textbox',[0.3 0.01 0.35 0.15],'string',...
    {['PLS model: sum-EEG(-500ms) -> Xcorr(0) in synco^{* PLOT 11-2}'],char(datetime('now'))});

% topoplot for leading state only in synch
canvas(0.3,0.2);
for states1=1;
cmin=-0.005;cmax=0.005;
c=[2];
for s=1    
    for b=1:5
        subplot(1,5,(s-1)*5+b)
        topoplot(plsmodel(c(s)).weights(b,:),channels,'nosedir','+X','style','map');
        clim([cmin cmax]);
    end
end
colormap(hnc)
cb=colorbar;
cb.AxisLocation = 'out';
cb.Position = [0.92 0.15 0.01 0.75];
set(gcf,'color','w'); % set background white for copying in ubuntu
delete(findall(gcf,'type','annotation'))
h0=annotation('textbox',[0.17 0.95 0.05 0.03],'string','Delta','color',[0 0 0])
h1=annotation('textbox',[0.33 0.95 0.05 0.03],'string','Theta','color',[0 0 0])
h2=annotation('textbox',[0.5 0.95 0.05 0.03],'string','Alpha','color',[0 0 0])
h3=annotation('textbox',[0.66 0.95 0.05 0.03],'string','Beta','color',[0 0 0])
h4=annotation('textbox',[0.82 0.95 0.05 0.03],'string','Gamma','color',[0 0 0])
v0=annotation('textbox',[0.1 0.4 0.05 0.03],'string',...
    {[states4names{c(1)}], [' (R^2= ' num2str(round(R2s(c(1)),1)) ...
            '  Fac= ' num2str(Fac) ') ']}...
            ,'color',condicolors(2,:));
set(v0,'Rotation',90);
end
set([h0 h1 h2 h3 h4 v0 v1 v2 v3], 'fitboxtotext','on',...
    'edgecolor','none')
sg=annotation('textbox',[0.3 0.01 0.35 0.15],'string',...
    {['PLS model: sum-EEG(-500ms) -> Xcorr(0) in synch^{* PLOT 11-2}'],char(datetime('now'))});
%% SECT 12 Compute EEG power (-/+ 500ms) 5 bands
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);
sessions={'synch','synco','synch','synco','synch','synco','synch','synco','synch','synco',...
    'synch','synco'};
EEGwin=0.5; % second

% EEG -500ms before the matched tap
zEEG500_L=cell(numSes,12);zEEG500_R=cell(numSes,12);
delta_L=cell(numSes,12);theta_L=cell(numSes,12);alpha_L=cell(numSes,12);beta_L=cell(numSes,12);gamma_L=cell(numSes,12);
delta_R=cell(numSes,12);theta_R=cell(numSes,12);alpha_R=cell(numSes,12);beta_R=cell(numSes,12);gamma_R=cell(numSes,12);
cd /ssd/zhibin/1overf/
tic
for s=1:numSes % each session
    runid=num2str(seeds(s,:));
    clear dataL dataR samples
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
toc % about 2 min
% clear working space
clear zEEG_L zEEG_R zEEG500_L zEEG500_R 
clear dataR dataL

% EEG +500ms After the matched tap
delta_LL=cell(numSes,12);theta_LL=cell(numSes,12);alpha_LL=cell(numSes,12);beta_LL=cell(numSes,12);gamma_LL=cell(numSes,12);
delta_RR=cell(numSes,12);theta_RR=cell(numSes,12);alpha_RR=cell(numSes,12);beta_RR=cell(numSes,12);gamma_RR=cell(numSes,12);
zEEG500_LL=cell(numSes,12);zEEG500_RR=cell(numSes,12);
cd /ssd/zhibin/1overf/
tic
for s=1:numSes % each session
    runid=num2str(seeds(s,:));
    clear dataL dataR
    load(['/ssd/zhibin/1overf/' runid '_2P/Cleaned_data/clean_' runid '.mat' ])
    zEEG_L={};zEEG_R={};
    for b=1:12 % each block
        zEEG_L{b}=zscore(dataL{b}(:,1:32),[],1);
        zEEG_R{b}=zscore(dataR{b}(:,1:32),[],1);
        for i=1:(size(samples{b},1)-2) % each matched interval (excluding the last tap or two)          
            zEEG500_LL{s,b}{i}=zEEG_L{b}(samples{b}(i,1):samples{b}(i,1)+999,:); % start from first matched tap 
            zEEG500_RR{s,b}{i}=zEEG_R{b}(samples{b}(i,2):samples{b}(i,2)+999,:); % end at the 2nd or 3rd last matched tap
            [delta_LL{s,b}(i,:), theta_LL{s,b}(i,:), alpha_LL{s,b}(i,:), beta_LL{s,b}(i,:), gamma_LL{s,b}(i,:)]...
                =sum5band(zEEG500_LL{s,b}{i},sr,EEGwin);
            [delta_RR{s,b}(i,:), theta_RR{s,b}(i,:), alpha_RR{s,b}(i,:), beta_RR{s,b}(i,:), gamma_RR{s,b}(i,:)]...
                =sum5band(zEEG500_RR{s,b}{i},sr,EEGwin);
        end
    end
end
toc  % about 2 min

%% SECT 12-1 Exact evoked potential to other's tapping (based on matched BP only)
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);
% erp_L={};erp_R={};
zEEG500_LL={};zEEG500_RR={}; % saving EEG for computing EEG across all subjects
tic
% extract EEG based on matched BP only
for s=1:12; % session 1-12
    runid=num2str(seeds(s,:));
    clear dataL dataR samples conditions
    load(['/ssd/zhibin/1overf/' runid '_2P/Cleaned_data/clean_' runid '.mat' ])
    zEEG_L={};zEEG_R={}; % saving individual blocks
    for b=1:12; % block 1 -12
    zEEG_L{b}=zscore(dataL{b}(:,1:32),[],1);
    zEEG_R{b}=zscore(dataR{b}(:,1:32),[],1);
            for i=1:(size(samples{b},1)-2) % each matched interval (excluding the last tap or two) 
                if conditions(b)==1;
            zEEG500_LL{s,b}(:,:,i)=zEEG_L{b}(samples{b}(i,1):samples{b}(i,1)+999,:); % start from first matched tap
            zEEG500_RR{s,b}(:,:,i)=zEEG_R{b}(samples{b}(i,2):samples{b}(i,2)+999,:); % end at the 2nd or 3rd last matched tap
                elseif conditions(b)==2;
            zEEG500_LL{s,b}(:,:,i)=zEEG_L{b}(samples{b}(i,1):samples{b}(i,1)+999,:); % start from first matched tap
            zEEG500_RR{s,b}(:,:,i)=zEEG_R{b}(samples{b}(i,1):samples{b}(i,1)+999,:); % end at the 2nd or 3rd last matched tap
                elseif conditions(b)==3;
            zEEG500_LL{s,b}(:,:,i)=zEEG_L{b}(samples{b}(i,2):samples{b}(i,2)+999,:); % start from first matched tap
            zEEG500_RR{s,b}(:,:,i)=zEEG_R{b}(samples{b}(i,2):samples{b}(i,2)+999,:); % end at the 2nd or 3rd last matched tap
                else conditions(b)==4;
            zEEG500_LL{s,b}(:,:,i)=zEEG_L{b}(samples{b}(i,2):samples{b}(i,2)+999,:); % start from first matched tap
            zEEG500_RR{s,b}(:,:,i)=zEEG_R{b}(samples{b}(i,1):samples{b}(i,1)+999,:); % end at the 2nd or 3rd last matched tap
                end
            end
            % erp_L{s,b}=mean(zEEG500_LL{s,b},3); % ERP in 32 channels
            % erp_R{s,b}=mean(zEEG500_RR{s,b},3);
    end
end
toc; % 108 second
% clear for workspace
clear dataL dataR erp_L erp_R 
% prepare indicies for conditions and states
condition_all; % from SECT 10-1
Inds4_LR=[uncoupleInd_LR leadingInd_LR followingInd_LR mutualInd_LR]; % from SECT 10-1
% reorganzie cell sequence for zEEG500_LL and zEEG500_RR
zEEG500_LL_all=reshape(zEEG500_LL',[],1);
zEEG500_RR_all=reshape(zEEG500_RR',[],1);
zEEG500_LLRR_all=[zEEG500_LL_all; zEEG500_RR_all]; % combine cell array
% clear for workspace
clear zEEG500_LL_all zEEG500_RR_all
% compute erp of 4 states across all sessions
% zEEG500_LL_4ss={};
% zEEG500_RR_4ss={};
canvas(0.5,1);
% selectChan=7; % F4
tic
for ss=1:4
    cat_stateEEG=[];
    cat_stateEEG=cat(3,zEEG500_LLRR_all{Inds4_LR(:,ss)});
    stateERP=mean(cat_stateEEG,3);
    for selectChan=1:32
        subplot(4,8,selectChan)
        plot([1:1000]/2,stateERP(:,selectChan),'color',condicolors(ss,:));
        hold on;
        if ss==4;
        title(['Chan' labels{selectChan}]); xlabel('time (ms)');
        xlim([0 500]);ylim([-0.2 0.2])
        end
    end
end
lg=legend(states4names,'location','eastoutside');
lg.Position = [0.9475 0.4 0.01 0.25];
set(gcf,'color','w'); % set background white for copying in ubuntu
sgtitle('ERP of 4 states')
toc % 7 sec
% clear for workspace
clear cat_stateEEG
% separateing syn types - all 32 chan
syn2names={'Synch','Synco'};
for t=1:2
    zEEG500_LL_4ss={};
    zEEG500_RR_4ss={};
    canvas(0.5,1);
    for ss=1:4
    cat_stateEEG=[];
    cat_stateEEG=cat(3,zEEG500_LLRR_all{Inds4_LR([synind(t,:) 36+synind(t,:)],ss)});
    stateERP=mean(cat_stateEEG,3);
    for selectChan=1:32
        subplot(4,8,selectChan)
        plot([1:1000]/2,stateERP(:,selectChan),'color',condicolors(ss,:));
        hold on;
        if ss==4;
        title(['Chan' labels{selectChan}]); xlabel('time (ms)');
        xlim([0 500]);ylim([-0.2 0.2])
        end
    end
    end
lg=legend(states4names,'location','eastoutside');
lg.Position = [0.9475 0.4 0.01 0.25];
set(gcf,'color','w'); % set background white for copying in ubuntu
sgtitle(['ERP of 4 states: ' syn2names{t}],'Color',syn2colors(t,:));
end
% clear for workspace
clear cat_stateEEG

% separateing syn types - chan 24-32
syn2names={'Synch','Synco'};
canvas(0.5,0.5);
for t=1:2
    % zEEG500_LL_4ss={};
    % zEEG500_RR_4ss={};
    for ss=1:4
    cat_stateEEG=[];
    cat_stateEEG=cat(3,zEEG500_LLRR_all{Inds4_LR([synind(t,:) 36+synind(t,:)],ss)});
    stateERP=mean(cat_stateEEG,3);
    for selectChan=24:32
        subplot(2,9,(t-1)*9+selectChan-23)
        plot([1:1000]/2,stateERP(:,selectChan),'color',condicolors(ss,:));
        hold on;
        if ss==4;
        title(['Chan' labels{selectChan}],'Color',syn2colors(t,:)); xlabel('time (ms)');
        xlim([0 500]);ylim([-0.2 0.25])
        end
    end
    end
lg=legend(states4names,'location','eastoutside');
lg.Position = [0.9475 0.4 0.01 0.25];
set(gcf,'color','w'); % set background white for copying in ubuntu
sgtitle('ERP of 4 states: Synch (top row); Synco (bottom row)');
end
% clear for workspace
clear cat_stateEEG

%% SECT 12-2 Exact evoked potential to other or one's own history (based on matched BP and pacers)

% use prevous synch session's sample as time 0 for uncouple and leading...

% or extract pacers from synco sessions (skip)
Pacers={}; % 12 sessions x 12 blocks
for s=1:12; % session 1-12
    runid=num2str(seeds(s,:));
    clear dataL dataR 
    load(['/ssd/zhibin/1overf/' runid '_2P/Cleaned_data/clean_' runid '.mat' ])
    % Pacers=?
end

% or extract pacers photocells (optimal, saved)
PacerTime01_all=cell(2,numSes,12);
cd /ssd/zhibin/1overf/
for r=1:numSes;
    clearvars -except PacerTime01_all r seeds numSes
    tic
    runid = num2str(seeds(r,:));
    path = [runid '_2P/'];
    % load([path runid 'workspace.mat']);
    load([path runid 'workspace.mat'],'PacerTimeL01','PacerTimeR01','SegIndL','SegIndR');
    for i=1:12
    PacerTime01_all{1,r,i}=PacerTimeL01(SegIndL(2*i-1):SegIndL(2*i),:); % in cells
    PacerTime01_all{2,r,i}=PacerTimeR01(SegIndR(2*i-1):SegIndR(2*i),:);
    end
    toc
end
clearvars -except PacerTime01_all 
cd /ssd/zhibin/1overf/all_session20220713_1005
save("PacerTime01_all.mat",'PacerTime01_all');
load('PacerTime01_all.mat');

% or extract showframes from synco session (skip)
cd /ssd/zhibin/1overf/
for r=1:numSes;
    clearvars -except PacerTimeL01
    runid = num2str(seeds(r,:));
    path = [runid '_2P/Segmented_data/'];
    load([path  'Showframes' runid '.mat']);
    Showframes(1)
    Showframes(2)
end

% extract FB and samples from all sessions (saved)
tic
FB_all=cell(2*numSes,1);
samples_all=cell(numSes,1);
cd /ssd/zhibin/1overf/
for r=1:numSes;
    clearvars -except FB_all samples_all seeds r numSes
    runid = num2str(seeds(r,:));
    path = [runid '_2P/Segmented_data/'];
    load([path  'data_variables' runid '.mat']);
    FB_all{r*2-1,1}=FB(1,1).FB;
    FB_all{r*2,1}=FB(2,1).FB;
    clearvars -except FB_all samples_all seeds runid r numSes
    load(['/ssd/zhibin/1overf/' runid '_2P/Cleaned_data/clean_' runid '.mat' ]);
    samples_all{r}=samples;
end
clearvars -except FB_all samples_all
cd /ssd/zhibin/1overf/all_session20220713_1005
save("FB_all_samples_all.mat",'FB_all','samples_all');
toc 232 sec

% matched PacerTime01_all to samples_all
PacerTime01_all; % 2 subj x 12 sessions x 12 blocks
samples_all; % 12 sessions (12 blocks; taps x 2 subj)
condition_all;
matched_Psamp=[];
for s=1:12
    for b=1:12
        tic
        BPsampL=samples_all{s}{b}(:,1);
        PsampL=find(PacerTime01_all{1,s,b}==1); 
        % only even# session (synco uncouple and leading) has pacers all the way to end of blocks
        % and mutual might have empty Psampl
        BPsampR=samples_all{s}{b}(:,2);
        PsampR=find(PacerTime01_all{2,s,b}==1);
        % matche the pacers to the already matched BPsamp
        % find the closest pacer before each matched BPsamp
        for i=1:size(samples_all{s}{b},1)
            clear all_positive_PsampL all_positive_PsampR
            if ~isempty(PsampL) % if empty 
            all_positive_PsampL=PsampL((BPsampL(i)-PsampL)>0); % all positive difference in PsampL
            matched_Psamp{s}{b}(i,1)=all_positive_PsampL(end);
            end
            if ~isempty(PsampR) % if only one value, there is 
            all_positive_PsampR=PsampR((BPsampR(i)-PsampR)>0); % all positive difference in PsampR
            matched_Psamp{s}{b}(i,2)=all_positive_PsampR(end);
            end
        end
        toc
    end
end
cd G:\1overf\all_session20220713_1005
save("matched_Psamp.mat",'matched_Psamp');
% clear workspace
clear PacerTime01_all FB_all zEEG500_LLRR_all

cd /ssd/zhibin/1overf/all_session20220713_1005
load('matched_Psamp.mat')

% extract EEG again (this time based on matched BP and Pacers in uncouple and leding in synco
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);
zEEG500_LL={};zEEG500_RR={}; % saving EEG for computing EEG across all subjects
tic
% extract EEG based on matched BP only
for s=1:12; % session 1-12
    runid=num2str(seeds(s,:));
    clear dataL dataR samples conditions
    load(['/ssd/zhibin/1overf/' runid '_2P/Cleaned_data/clean_' runid '.mat' ])
    zEEG_L={};zEEG_R={}; % saving individual blocks
    for b=1:12; % block 1 -12
    zEEG_L{b}=zscore(dataL{b}(:,1:32),[],1);
    zEEG_R{b}=zscore(dataR{b}(:,1:32),[],1);
            for i=1:(size(samples{b},1)-2) % each matched interval (excluding the last tap or two) 
                if conditions(b)==1 & mod(s,2)==1 % synch
            zEEG500_LL{s,b}(:,:,i)=zEEG_L{b}(samples{b}(i,1):samples{b}(i,1)+999,:); % start from first matched tap
            zEEG500_RR{s,b}(:,:,i)=zEEG_R{b}(samples{b}(i,2):samples{b}(i,2)+999,:); % end at the 2nd or 3rd last matched tap
                elseif conditions(b)==1 & mod(s,2)==0 % synco
            zEEG500_LL{s,b}(:,:,i)=zEEG_L{b}(matched_Psamp{s}{b}(i,1):matched_Psamp{s}{b}(i,1)+999,:); % start from first matched tap
            zEEG500_RR{s,b}(:,:,i)=zEEG_R{b}(matched_Psamp{s}{b}(i,2):matched_Psamp{s}{b}(i,2)+999,:); % end at the 2nd or 3rd last matched tap
                elseif conditions(b)==2 & mod(s,2)==1 % synch
            zEEG500_LL{s,b}(:,:,i)=zEEG_L{b}(samples{b}(i,1):samples{b}(i,1)+999,:); % start from first matched tap
            zEEG500_RR{s,b}(:,:,i)=zEEG_R{b}(samples{b}(i,1):samples{b}(i,1)+999,:); % end at the 2nd or 3rd last matched tap
                elseif conditions(b)==2 & mod(s,2)==0 % synco
            zEEG500_LL{s,b}(:,:,i)=zEEG_L{b}(matched_Psamp{s}{b}(i,1):matched_Psamp{s}{b}(i,1)+999,:); % start from first matched tap
            zEEG500_RR{s,b}(:,:,i)=zEEG_R{b}(samples{b}(i,1):samples{b}(i,1)+999,:); % end at the 2nd or 3rd last matched tap
                elseif conditions(b)==3 & mod(s,2)==1 % synch
            zEEG500_LL{s,b}(:,:,i)=zEEG_L{b}(samples{b}(i,2):samples{b}(i,2)+999,:); % start from first matched tap
            zEEG500_RR{s,b}(:,:,i)=zEEG_R{b}(samples{b}(i,2):samples{b}(i,2)+999,:); % end at the 2nd or 3rd last matched tap
                elseif conditions(b)==3 & mod(s,2)==0 % synco
            zEEG500_LL{s,b}(:,:,i)=zEEG_L{b}(samples{b}(i,2):samples{b}(i,2)+999,:); % start from first matched tap
            zEEG500_RR{s,b}(:,:,i)=zEEG_R{b}(matched_Psamp{s}{b}(i,2):matched_Psamp{s}{b}(i,2)+999,:); % end at the 2nd or 3rd last matched tap
                else conditions(b)==4;
            zEEG500_LL{s,b}(:,:,i)=zEEG_L{b}(samples{b}(i,2):samples{b}(i,2)+999,:); % start from first matched tap
            zEEG500_RR{s,b}(:,:,i)=zEEG_R{b}(samples{b}(i,1):samples{b}(i,1)+999,:); % end at the 2nd or 3rd last matched tap
                end
            end
            % erp_L{s,b}=mean(zEEG500_LL{s,b},3); % ERP in 32 channels
            % erp_R{s,b}=mean(zEEG500_RR{s,b},3);
    end
end
toc; %  second
% clear for workspace
clear dataL dataR zEEG_L zEEG_R erp_L erp_R 
% prepare indicies for conditions and states
condition_all; % from SECT 10-1
Inds4_LR=[uncoupleInd_LR leadingInd_LR followingInd_LR mutualInd_LR]; % from SECT 10-1
% reorganzie cell sequence for zEEG500_LL and zEEG500_RR
zEEG500_LL_all=reshape(zEEG500_LL',[],1);
zEEG500_RR_all=reshape(zEEG500_RR',[],1);
zEEG500_LLRR_all=[zEEG500_LL_all; zEEG500_RR_all]; % combine cell array
% clear for workspace
clear zEEG500_LL_all zEEG500_RR_all
% save 
cd /ssd/zhibin/1overf/all_session20220713_1005
tic
save('zEEG500_LLRR_all.mat','zEEG500_LLRR_all');
toc % 271 sec

% separateing syn types - chan 24-32
syn2names={'Synch','Synco'};
figure; % canvas(0.5,0.5);
for t=1:2
    % zEEG500_LL_4ss={};
    % zEEG500_RR_4ss={};
    for ss=1:4
    cat_stateEEG=[];
    cat_stateEEG=cat(3,zEEG500_LLRR_all{Inds4_LR([synind(t,:) 36+synind(t,:)],ss)});
    stateERP=mean(cat_stateEEG,3);
    for selectChan=24:32
        subplot(2,9,(t-1)*9+selectChan-23)
        plot([1:1000]/2,stateERP(:,selectChan),'color',condicolors(ss,:));
        hold on;
        if ss==4;
        title(['Chan' labels{selectChan}],'Color',syn2colors(t,:)); xlabel('time (ms)');
        xlim([0 500]);ylim([-0.2 0.28])
        end
    end
    end
lg=legend(states4names,'location','eastoutside');
lg.Position = [0.9475 0.4 0.01 0.25];
set(gcf,'color','w'); % set background white for copying in ubuntu
sgtitle('ERP of 4 states: Synch (top row); Synco (bottom row)');
end
% clear for workspace
clear cat_stateEEG

%% SECT 13 Examine & Organize sum-EEG power and H for corr and PLS
% organize EEG power
% zEEG500_L; zEEG500_R; % for all sessions from SECT 12
delta_L;theta_L;alpha_L;beta_L;gamma_L;
delta_R;theta_R;alpha_R;beta_R;gamma_R;

% Eaxamine corr of EEG for each tap in trial level
% organize 12 x 12 cell arrays into 144 x 1 cell array (288 x 1 combining LR)
tic
delta_LR_all=[reshape(delta_L',[],1);reshape(delta_R',[],1)]; 
theta_LR_all=[reshape(theta_L',[],1);reshape(theta_R',[],1)]; 
alpha_LR_all=[reshape(alpha_L',[],1);reshape(alpha_R',[],1)]; 
beta_LR_all=[reshape(beta_L',[],1);reshape(beta_R',[],1)]; 
gamma_LR_all=[reshape(gamma_L',[],1);reshape(gamma_R',[],1)]; 
LR_all_chan=cat(2,delta_LR_all,theta_LR_all,alpha_LR_all,beta_LR_all,gamma_LR_all);
toc
clear delta_L theta_L alpha_L beta_L gamma_L delta_R theta_R alpha_R beta_R gamma_R;
clear delta_LR_all theta_LR_all alpha_LR_all beta_LR_all gamma_LR_all
% select trial and compute the corr
tic
Corr4x5_all=cell(2,3,5,18);
for t=1:2 % t=1 % synch/o
for s=[1 2 4] % s=4 % one of 4 states
for f=1:5 % f=2 % frequency band
for subj=1:18 % subj=1; %trials
    Corr4x5=[];
    if s==1;
     Corr4x5=...
    corr( ...
    LR_all_chan{(Inds4_LR(synind(t,subj),s)),f},...
    LR_all_chan{(Inds4_LR(36+synind(t,subj),s)),f}...
    );
     Corr4x5=(Corr4x5+Corr4x5')./2;
    Corr4x5_all{t,1,f,subj}=Corr4x5;
    elseif s==4;
     Corr4x5=...
    corr( ...
    LR_all_chan{(Inds4_LR(synind(t,subj),s)),f},...
    LR_all_chan{(Inds4_LR(36+synind(t,subj),s)),f}...
    );
     Corr4x5=(Corr4x5+Corr4x5')./2;
    Corr4x5_all{t,3,f,subj}=Corr4x5;
    else s==2 | s==3;
     Corr4x5=...
    corr( ...
    [LR_all_chan{(Inds4_LR(synind(t,subj),2)),f}; LR_all_chan{(Inds4_LR(36+synind(t,subj),2)),f}],...
    [LR_all_chan{(Inds4_LR(36+synind(t,subj),3)),f}; LR_all_chan{(Inds4_LR(synind(t,subj),3)),f}]...
    );   
    Corr4x5_all{t,2,f,subj}=Corr4x5;
    end
end
end
end
end
toc
% plot
for t=1:2 % t=1 % synch/o
    cmin=-0.6;cmax=0.6;
    figure
    clf
for s=1:3 % one of 4 states
for f=1:5 % f
    Corr4x5select=[];
    for subj=1:18    
        Corr4x5select=cat(3, Corr4x5select,Corr4x5_all{t,s,f,subj}); 
    end
    % average of corrcoef 
    Corr4x5m=mean(Corr4x5select,3);
    subplot(3,5,5*(s-1)+f)
    imagesc(Corr4x5m);colorbar; title(band5names{f})
    if s==1;  ylabel('Subject A');xlabel('Subject B');subtitle('Independent');
    elseif s==2; ylabel('Leader');xlabel('Follower');subtitle('Unidirectional');
    else s==3; ylabel('Subject A');xlabel('Subject B');subtitle('Bidirectional');
    end
    colormap(hnc)
    clim([cmin cmax]);
end
end
    if t==1; sgtitle(['synch ']); else t==2; sgtitle(['synco']);end;
end
% only select gamma to display
for t=1:2
    cmin=-0.5; cmax=0.5;
end


% sum the band power in each channel in each block (8x12=96)
delta_L_sum=[];theta_L_sum=[];alpha_L_sum=[];beta_L_sum=[];gamma_L_sum=[];
delta_R_sum=[];theta_R_sum=[];alpha_R_sum=[];beta_R_sum=[];gamma_R_sum=[];
for s=1:numSes
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
% clear out some variables for workspace
clear delta_L theta_L alpha_L beta_L gamma_L delta_R theta_R alpha_R beta_R gamma_R
% squeeze into 32 vectors for corr with H
delta_L_chan=[];theta_L_chan=[];alpha_L_chan=[];beta_L_chan=[];gamma_L_chan=[];
delta_R_chan=[];theta_R_chan=[];alpha_R_chan=[];beta_R_chan=[];gamma_R_chan=[];
for c=1:32
    delta_L_chan(:,c)=reshape(delta_L_sum(:,:,c)',[],1); % 144x32 (each element from one block in time sequence) 
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
% clear out variables for workspace
clear delta_L_sum theta_L_sum alpha_L_sum beta_L_sum gamma_L_sum delta_R_sum theta_R_sum alpha_R_sum beta_R_sum gamma_R_sum 
% Combine L and R for 288 predictors in PLS
delta_LR_chan=[delta_L_chan;delta_R_chan]; % (288 x 32)
theta_LR_chan=[theta_L_chan;theta_R_chan];
alpha_LR_chan=[alpha_L_chan;alpha_R_chan];
beta_LR_chan=[beta_L_chan;beta_R_chan];
gamma_LR_chan=[gamma_L_chan;gamma_R_chan];
% fix the scale in the data (as in PLOT 16, now skiped in PLOT 16)
delta_LR_chan = delta_LR_chan./(ones(288,1)*std(delta_LR_chan)); % organized in SECT 13
theta_LR_chan = theta_LR_chan./(ones(288,1)*std(theta_LR_chan));
alpha_LR_chan = alpha_LR_chan./(ones(288,1)*std(alpha_LR_chan));
beta_LR_chan = beta_LR_chan./(ones(288,1)*std(beta_LR_chan));
gamma_LR_chan = gamma_LR_chan./(ones(288,1)*std(gamma_LR_chan));

% Examine the power according to condition and session type
synind; %{ indicies of synind from PLOT 10-2 %}

% Separating L/R
cmin=-4;cmax=4;
for t=1:2 % type of session (synch/o)
    canvas(0.65,0.5);
    for s=1:4 % 4 states
    % Delta
    % Left suj
    subplot(4,10,10*(s-1)+1)
    topoplot(mean(delta_LR_chan(Inds4_LR(synind(t,:),s),:),1),channels,'nosedir','+X','style','map');
    % colorbar;
    clim([cmin cmax]);
    % Right suj
    subplot(4,10,10*(s-1)+2)
    topoplot(mean(delta_LR_chan(Inds4_LR(36+synind(t,:),s),:),1),channels,'nosedir','+X','style','map');
    % colorbar;
    clim([cmin cmax]);
    % Theta
    % Left suj
    subplot(4,10,10*(s-1)+3)
    topoplot(mean(theta_LR_chan(Inds4_LR(synind(t,:),s),:),1),channels,'nosedir','+X','style','map');
    % colorbar;
    clim([cmin cmax]);
    % Right suj
    subplot(4,10,10*(s-1)+4)
    topoplot(mean(theta_LR_chan(Inds4_LR(36+synind(t,:),s),:),1),channels,'nosedir','+X','style','map');
    % colorbar;
    clim([cmin cmax]);
    % Alpha
    % Left suj
    subplot(4,10,10*(s-1)+5)
    topoplot(mean(alpha_LR_chan(Inds4_LR(synind(t,:),s),:),1),channels,'nosedir','+X','style','map');
    % colorbar;
    clim([cmin cmax]);
    % Right suj
    subplot(4,10,10*(s-1)+6)
    topoplot(mean(alpha_LR_chan(Inds4_LR(36+synind(t,:),s),:),1),channels,'nosedir','+X','style','map');
    % colorbar;
    clim([cmin cmax]);
    % Beta
    % Left suj
    subplot(4,10,10*(s-1)+7)
    topoplot(mean(beta_LR_chan(Inds4_LR(synind(t,:),s),:),1),channels,'nosedir','+X','style','map');
    % colorbar;
    clim([cmin cmax]);
    % Right suj
    subplot(4,10,10*(s-1)+8)
    topoplot(mean(beta_LR_chan(Inds4_LR(36+synind(t,:),s),:),1),channels,'nosedir','+X','style','map');
    % colorbar;
    clim([cmin cmax]);
    % Gamma
    % Left suj
    subplot(4,10,10*(s-1)+9)
    topoplot(mean(gamma_LR_chan(Inds4_LR(synind(t,:),s),:),1),channels,'nosedir','+X','style','map');
    % colorbar;
    clim([cmin cmax]);
    % Right suj
    subplot(4,10,10*(s-1)+10)
    topoplot(mean(gamma_LR_chan(Inds4_LR(36+synind(t,:),s),:),1),channels,'nosedir','+X','style','map');
    % colorbar;
    clim([cmin cmax]);
    end
    colormap(hnc)
    if true % annotation
    cb=colorbar;
    cb.AxisLocation = 'out';
    cb.Position = [0.92 0.15 0.01 0.75];
    delete(findall(gcf,'type','annotation'))
    h0=annotation('textbox',[0.15 0.95 0.05 0.03],'string',{'Delta','L'},'color',[0 0 0]);
    h01=annotation('textbox',[0.23 0.95 0.05 0.03],'string',{'Delta','R'},'color',[0 0 0]);
    h1=annotation('textbox',[0.31 0.95 0.05 0.03],'string',{'Theta','L'},'color',[0 0 0]);
    h11=annotation('textbox',[0.39 0.95 0.05 0.03],'string',{'Theta','R'},'color',[0 0 0]);
    h2=annotation('textbox',[0.47 0.95 0.05 0.03],'string',{'Alpha','L'},'color',[0 0 0]);
    h21=annotation('textbox',[0.55 0.95 0.05 0.03],'string',{'Alpha','R'},'color',[0 0 0]);
    h3=annotation('textbox',[0.63 0.95 0.05 0.03],'string',{'Beta','L'},'color',[0 0 0]);
    h31=annotation('textbox',[0.71 0.95 0.05 0.03],'string',{'Beta','R'},'color',[0 0 0]);
    h4=annotation('textbox',[0.79 0.95 0.05 0.03],'string',{'Gamma','L'},'color',[0 0 0]);
    h41=annotation('textbox',[0.86 0.95 0.05 0.03],'string',{'Gamma','R'},'color',[0 0 0]);
    v0=annotation('textbox',[0.12 0.15 0.05 0.03],'string','Mutual','color',condicolors(4,:));
    v1=annotation('textbox',[0.12 0.37 0.05 0.03],'string','Following','color',condicolors(3,:));
    v2=annotation('textbox',[0.12 0.59 0.05 0.03],'string','Leading','color',condicolors(2,:));
    v3=annotation('textbox',[0.12 0.81 0.05 0.03],'string','Uncouple','color',condicolors(1,:));
    set(v0,'Rotation',90);set(v1,'Rotation',90);set(v2,'Rotation',90);set(v3,'Rotation',90);
    set([h0 h1 h2 h3 h4 v0 v1 v2 v3 h01 h11 h21 h31 h41], 'fitboxtotext','on',...
        'edgecolor','none')
    set(gcf,'color','w'); % set background white for copying in ubuntu
    end
    % title annotation
    % delete(sg);
    if t == 1
        sg=annotation('textbox',[0.3 0.01 0.4 0.07],'string',...
        {'Average of sum-EEG power (-500ms) in synch ^{* PLOT 13}',char(datetime('now'))});
    else % t == 2; 
        sg=annotation('textbox',[0.3 0.01 0.4 0.07],'string',...
        {'Average of sum-EEG power (-500ms) in synco ^{* PLOT 13}',char(datetime('now'))});
    end
end

% Correlation L & R (32 corrcoef values (chan to chan, between chan pairs), 2 figures of topoplots)
% Compute corrleation
% cmin=-1;cmax=1;
Corr4x5=nan(2,4,5,32); % Corr coef for 2 syn_types x 4 states x 5 freq x 32 chan
f5_LR_chan=cat(3,delta_LR_chan,theta_LR_chan,alpha_LR_chan,beta_LR_chan,gamma_LR_chan);
for t=1:2 % type of session (synch/o)
    % canvas(0.65,0.5);
    for s=1:4 % 4 states
        for f=1:5 % 5 freq
    % Delta
    % Left suj vs Right suj
    % subplot(4,5,5*(s-1)+1)
    % delta_corr=diag(corr(delta_LR_chan(Inds4_LR(synind(t,:),s),:),delta_LR_chan(Inds4_LR(36+synind(t,:),s),:)));
    % topoplot(delta_corr,channels);
    % colorbar;
    % clim([cmin cmax]);
        fselected_LR_chan=[];
        fselected_LR_chan=f5_LR_chan(:,:,f);
        Corr4x5(t,s,f,:)=...
            diag(corr(...
            [fselected_LR_chan(Inds4_LR(synind(t,:),s),:)],...
            [fselected_LR_chan(Inds4_LR(36+synind(t,:),s),:)]...
            ));
        % the above corr is worng, it only include half of the data in the unidirectional case? might
        % need to correct them into 3 states as well
        end
    end
end
% Plot
cmin=-0.6;cmax=0.6;
for t=1:2;
    canvas(0.3,0.5);
    for s=1:4;
        for f=1:5;
subplot(4,5,5*(s-1)+f)
topoplot(Corr4x5(t,s,f,:),channels,'nosedir','+X');
clim([cmin cmax]);
        end
    end
    colormap(hnc)
    if true % annotation
    cb=colorbar;
    cb.AxisLocation = 'out';
    cb.Position = [0.92 0.15 0.01 0.75];
    delete(findall(gcf,'type','annotation'))
    h0=annotation('textbox',[0.17 0.95 0.05 0.03],'string','Delta','color',[0 0 0]);
    h1=annotation('textbox',[0.33 0.95 0.05 0.03],'string','Theta','color',[0 0 0]);
    h2=annotation('textbox',[0.5 0.95 0.05 0.03],'string','Alpha','color',[0 0 0]);
    h3=annotation('textbox',[0.66 0.95 0.05 0.03],'string','Beta','color',[0 0 0]);
    h4=annotation('textbox',[0.81 0.95 0.05 0.03],'string','Gamma','color',[0 0 0]);
    v0=annotation('textbox',[0.14 0.15 0.05 0.03],'string','Mutual','color',condicolors(4,:));
    v1=annotation('textbox',[0.14 0.37 0.05 0.03],'string','Following','color',condicolors(3,:));
    v2=annotation('textbox',[0.14 0.59 0.05 0.03],'string','Leading','color',condicolors(2,:));
    v3=annotation('textbox',[0.14 0.81 0.05 0.03],'string','Uncouple','color',condicolors(1,:));
    set(v0,'Rotation',90);set(v1,'Rotation',90);set(v2,'Rotation',90);set(v3,'Rotation',90);
    set([h0 h1 h2 h3 h4 v0 v1 v2 v3], 'fitboxtotext','on',...
    'edgecolor','none')
    set(gcf,'color','w'); % set background white for copying in ubuntu
    end
    % title annotation
    % delete(sg);
    if t == 1
        sg=annotation('textbox',[0.3 0.01 0.4 0.07],'string',...
        {'Corr of sum-EEG power (-500ms) in synch ^{* PLOT 13}',char(datetime('now'))},'Color',syn2colors(1,:));
    else % t == 2; 
        sg=annotation('textbox',[0.3 0.01 0.4 0.07],'string',...
        {'Corr of sum-EEG power (-500ms) in synco ^{* PLOT 13}',char(datetime('now'))},'Color',syn2colors(2,:));
    end
end

% Correlation L&R (32x32 corrcoef values, 2 figures of imagesc)
% Compute corrleation
% cmin=-1;cmax=1;
Corr4x5=[];
Corr4x5=nan(2,3,5,32,32); % Corr coef for 2 syn_types x 3 states x 5 freq x 32 chan x32 chan
f5_LR_chan=cat(3,delta_LR_chan,theta_LR_chan,alpha_LR_chan,beta_LR_chan,gamma_LR_chan);
for t=1:2 % type of session (synch/o)
    for s=1:3 % 3 states
        for f=1:5 % 5 freq
        fselected_LR_chan=[];
        fselected_LR_chan=f5_LR_chan(:,:,f);
        if s==1; % corr of L and R for uncoupled state
        Corr4x5(t,s,f,:,:)=...
            corr(fselected_LR_chan(Inds4_LR(synind(t,:),s),:),fselected_LR_chan(Inds4_LR(36+synind(t,:),s),:));
        elseif s==3; % corr of L and R for bidirectional state
        Corr4x5(t,s,f,:,:)=...
            corr(fselected_LR_chan(Inds4_LR(synind(t,:),4),:),fselected_LR_chan(Inds4_LR(36+synind(t,:),4),:));
        else s==2; % corr of leader and follwer for unidiretional state
        Corr4x5(t,s,f,:,:)=...
            corr(...
            [fselected_LR_chan(Inds4_LR(synind(t,:),2),:); fselected_LR_chan(Inds4_LR(36+synind(t,:),2),:)],...% leaders
            [fselected_LR_chan(Inds4_LR(36+synind(t,:),3),:); fselected_LR_chan(Inds4_LR(synind(t,:),3),:)]... % followers
            );    
        end
        end
    end
end
% Plot
cmin=-0.6;cmax=0.6;
for t=1:2;
    canvas(0.3,0.7);
    for s=1:3;
        for f=1:5;
        subplot(3,5,5*(s-1)+f)
        imagesc(squeeze(Corr4x5(t,s,f,:,:)));
        clim([cmin cmax]);
        if s==1 | s==3;
        imagesc([squeeze(Corr4x5(t,s,f,:,:))+(squeeze(Corr4x5(t,s,f,:,:)))']./2);
        clim([cmin cmax]);
        ylabel('chans (Participant A)');xlabel('chans (Participant B)');
        else s==2;
        imagesc(squeeze(Corr4x5(t,s,f,:,:)));
        clim([cmin cmax]);
        ylabel('chans (Leader)');xlabel('chans (Follower)');
        end
        end
    end
    colormap(hnc)
    if true % annotation
    cb=colorbar;
    cb.AxisLocation = 'out';
    cb.Position = [0.92 0.15 0.01 0.75];
    delete(findall(gcf,'type','annotation'))
    h0=annotation('textbox',[0.17 0.95 0.05 0.03],'string','Delta','color',[0 0 0]);
    h1=annotation('textbox',[0.33 0.95 0.05 0.03],'string','Theta','color',[0 0 0]);
    h2=annotation('textbox',[0.5 0.95 0.05 0.03],'string','Alpha','color',[0 0 0]);
    h3=annotation('textbox',[0.66 0.95 0.05 0.03],'string','Beta','color',[0 0 0]);
    h4=annotation('textbox',[0.81 0.95 0.05 0.03],'string','Gamma','color',[0 0 0]);
    v0=annotation('textbox',[0.1 0.15 0.05 0.03],'string','Bidirectional','color',condicolors(4,:));
    v1=annotation('textbox',[0.1 0.45 0.05 0.03],'string','Unidirectional','color',condicolors(3,:));
    % v2=annotation('textbox',[0.1 0.59 0.05 0.03],'string','Leading','color',condicolors(2,:));
    v3=annotation('textbox',[0.1 0.75 0.05 0.03],'string','Independent','color',condicolors(1,:));
    set(v0,'Rotation',90);set(v1,'Rotation',90);% set(v2,'Rotation',90);
    set(v3,'Rotation',90);
    % set([h0 h1 h2 h3 h4 v0 v1 v2 v3], 'fitboxtotext','on','edgecolor','none')
    set([h0 h1 h2 h3 h4 v0 v1 v3], 'fitboxtotext','on','edgecolor','none')
    set(gcf,'color','w'); % set background white for copying in ubuntu
    end
    % title annotation
    % delete(sg);
    if t == 1
        sg=annotation('textbox',[0.3 0.01 0.4 0.05],'string',...
        {'Corr of sum-EEG power (-500ms) in synch ^{* SECT 13}',char(datetime('now'))},'Color',syn2colors(1,:));
    else % t == 2; 
        sg=annotation('textbox',[0.3 0.01 0.4 0.05],'string',...
        {'Corr of sum-EEG power (-500ms) in synco ^{* SECT 13}',char(datetime('now'))},'Color',syn2colors(2,:));
    end
end

% Correlation L/R individually (32x32 individual corrcoef values, 4 figures)
% Compute corrleation
% cmin=-1;cmax=1;
Corr4x5=[];
Corr4x5=nan(2,2,4,5,32,32); % Corr coef for 2 syn_types x 2 subjects (L/R) X 4 states x 5 freq x 32 chan x32 chan
f5_LR_chan=cat(3,delta_LR_chan,theta_LR_chan,alpha_LR_chan,beta_LR_chan,gamma_LR_chan);
for t=1:2 % type of session (synch/o)
    for p=1:2 % two subjects
        for s=1:4 % 4 states
            for f=1:5 % 5 freq
            fselected_LR_chan=[];
            fselected_LR_chan=f5_LR_chan(:,:,f);
            if p==1; % Left subject
                Corr4x5(t,p,s,f,:,:)=...
                corr(fselected_LR_chan(Inds4_LR(synind(t,:),s),:),fselected_LR_chan(Inds4_LR(synind(t,:),s),:));
            else p==2; % Right subject
                Corr4x5(t,p,s,f,:,:)=...
                corr(fselected_LR_chan(Inds4_LR(36+synind(t,:),s),:),fselected_LR_chan(Inds4_LR(36+synind(t,:),s),:));
            end
            end
        end
    end
end
% Plot
cmin=-1;cmax=1;
for t=1:2;
    for p=1:2;
    canvas(0.3,0.7);
    for s=1:4;
        for f=1:5;
        subplot(4,5,5*(s-1)+f)
        imagesc(squeeze(Corr4x5(t,p,s,f,:,:)));
        clim([cmin cmax]);
        ylabel('chans');xlabel('chans');
        end
    end
    colormap(hnc)
    if true % annotation
    cb=colorbar;
    cb.AxisLocation = 'out';
    cb.Position = [0.92 0.15 0.01 0.75];
    delete(findall(gcf,'type','annotation'))
    h0=annotation('textbox',[0.17 0.95 0.05 0.03],'string','Delta','color',[0 0 0]);
    h1=annotation('textbox',[0.33 0.95 0.05 0.03],'string','Theta','color',[0 0 0]);
    h2=annotation('textbox',[0.5 0.95 0.05 0.03],'string','Alpha','color',[0 0 0]);
    h3=annotation('textbox',[0.66 0.95 0.05 0.03],'string','Beta','color',[0 0 0]);
    h4=annotation('textbox',[0.81 0.95 0.05 0.03],'string','Gamma','color',[0 0 0]);
    v0=annotation('textbox',[0.1 0.15 0.05 0.03],'string','Mutual','color',condicolors(4,:));
    v1=annotation('textbox',[0.1 0.37 0.05 0.03],'string','Following','color',condicolors(3,:));
    v2=annotation('textbox',[0.1 0.59 0.05 0.03],'string','Leading','color',condicolors(2,:));
    v3=annotation('textbox',[0.1 0.81 0.05 0.03],'string','Uncouple','color',condicolors(1,:));
    set(v0,'Rotation',90);set(v1,'Rotation',90);set(v2,'Rotation',90);set(v3,'Rotation',90);
    set([h0 h1 h2 h3 h4 v0 v1 v2 v3], 'fitboxtotext','on',...
    'edgecolor','none')
    set(gcf,'color','w'); % set background white for copying in ubuntu
    end
    % title annotation
    % delete(sg);
    if t == 1 & p==1;
        sg=annotation('textbox',[0.3 0.01 0.4 0.05],'string',...
        {'Corr of sum-EEG power (-500ms) in L synch ^{* SECT 13}',char(datetime('now'))},'Color',syn2colors(1,:));
    elseif t == 1 & p==2;
        sg=annotation('textbox',[0.3 0.01 0.4 0.05],'string',...
        {'Corr of sum-EEG power (-500ms) in R synch ^{* SECT 13}',char(datetime('now'))},'Color',syn2colors(1,:));
    elseif t == 2 & p==1;
        sg=annotation('textbox',[0.3 0.01 0.4 0.05],'string',...
        {'Corr of sum-EEG power (-500ms) in L synco ^{* SECT 13}',char(datetime('now'))},'Color',syn2colors(2,:));
    else t == 2 & p==2; 
        sg=annotation('textbox',[0.3 0.01 0.4 0.05],'string',...
        {'Corr of sum-EEG power (-500ms) in R synco ^{* SECT 13}',char(datetime('now'))},'Color',syn2colors(2,:));
    end
    end
end

% Correlation L&R combined (32x32 individual corrcoef values, 2 figures)
% Compute corrleation
% cmin=-1;cmax=1;
Corr4x5=[];
Corr4x5=nan(2,4,5,32,32); % Corr coef for 2 syn_types X 4 states x 5 freq x 32 chan x32 chan
f5_LR_chan=cat(3,delta_LR_chan,theta_LR_chan,alpha_LR_chan,beta_LR_chan,gamma_LR_chan);
for t=1:2 % type of session (synch/o)
    for s=1:4 % 4 states
        for f=1:5 % 5 freq
        fselected_LR_chan=[];
        fselected_LR_chan=f5_LR_chan(:,:,f);
        Corr4x5(t,s,f,:,:)=...
        corr(fselected_LR_chan(Inds4_LR([synind(t,:) 36+synind(t,:)],s),:),fselected_LR_chan(Inds4_LR([synind(t,:) 36+synind(t,:)],s),:));
        end
    end
end
% Plot
cmin=-1;cmax=1;
for t=1:2;
    canvas(0.3,0.7);
    for s=1:4;
        for f=1:5;
        subplot(4,5,5*(s-1)+f)
        imagesc(squeeze(Corr4x5(t,s,f,:,:)));
        clim([cmin cmax]);
        ylabel('chans');xlabel('chans');
        end
    end
    colormap(hnc)
    if true % annotation
    cb=colorbar;
    cb.AxisLocation = 'out';
    cb.Position = [0.92 0.15 0.01 0.75];
    delete(findall(gcf,'type','annotation'))
    h0=annotation('textbox',[0.17 0.95 0.05 0.03],'string','Delta','color',[0 0 0]);
    h1=annotation('textbox',[0.33 0.95 0.05 0.03],'string','Theta','color',[0 0 0]);
    h2=annotation('textbox',[0.5 0.95 0.05 0.03],'string','Alpha','color',[0 0 0]);
    h3=annotation('textbox',[0.66 0.95 0.05 0.03],'string','Beta','color',[0 0 0]);
    h4=annotation('textbox',[0.81 0.95 0.05 0.03],'string','Gamma','color',[0 0 0]);
    v0=annotation('textbox',[0.1 0.15 0.05 0.03],'string','Mutual','color',condicolors(4,:));
    v1=annotation('textbox',[0.1 0.37 0.05 0.03],'string','Following','color',condicolors(3,:));
    v2=annotation('textbox',[0.1 0.59 0.05 0.03],'string','Leading','color',condicolors(2,:));
    v3=annotation('textbox',[0.1 0.81 0.05 0.03],'string','Uncouple','color',condicolors(1,:));
    set(v0,'Rotation',90);set(v1,'Rotation',90);set(v2,'Rotation',90);set(v3,'Rotation',90);
    set([h0 h1 h2 h3 h4 v0 v1 v2 v3], 'fitboxtotext','on',...
    'edgecolor','none')
    set(gcf,'color','w'); % set background white for copying in ubuntu
    end
    % title annotation
    % delete(sg);
    if t == 1;
        sg=annotation('textbox',[0.3 0.01 0.4 0.05],'string',...
        {'Corr of sum-EEG power (-500ms) in synch ^{* SECT 13}',char(datetime('now'))},'Color',syn2colors(1,:));
    else t == 2;
        sg=annotation('textbox',[0.3 0.01 0.4 0.05],'string',...
        {'Corr of sum-EEG power (-500ms) in synco ^{* SECT 13}',char(datetime('now'))},'Color',syn2colors(2,:));
    end
end

% Organize H_all for corr
H_all; % (2xnumSesx12) for all sessions from SECT 10-1 (matched int)
H_all_L=squeeze(H_all(1,:,:));
H_all_R=squeeze(H_all(2,:,:));
% squeeze into 1 vector from the 96 blocks for each subject, for corr with pow in each chan
H_all_L=reshape(H_all_L',[],1);% 144x1 (each element from one block in time sequence) 
H_all_R=reshape(H_all_R',[],1);
% Combine L and R
H_all_LR=[H_all_L;H_all_R]; % to be used for corr and PLS
%% PLOT 13 corr of sum-EEG power (-500ms) and H-int all sessions
% Compute the correlation between sum-EEG pow and H-int
for c=1:32
    delta_L_H_L_corr(c)=corr(delta_LR_chan(1:144,c),H_all_L);
    theta_L_H_L_corr(c)=corr(theta_LR_chan(1:144,c),H_all_L);
    alpha_L_H_L_corr(c)=corr(alpha_LR_chan(1:144,c),H_all_L);
    beta_L_H_L_corr(c)=corr(beta_LR_chan(1:144,c),H_all_L);
    gamma_L_H_L_corr(c)=corr(gamma_LR_chan(1:144,c),H_all_L);
    delta_R_H_R_corr(c)=corr(delta_LR_chan(145:288,c),H_all_R);
    theta_R_H_R_corr(c)=corr(theta_LR_chan(145:288,c),H_all_R);
    alpha_R_H_R_corr(c)=corr(alpha_LR_chan(145:288,c),H_all_R);
    beta_R_H_R_corr(c)=corr(beta_LR_chan(145:288,c),H_all_R);
    gamma_R_H_R_corr(c)=corr(gamma_LR_chan(145:288,c),H_all_R);
end
% Topoplots (5x2)
canvas(0.2,0.6);
for topo5x2=1;
cmin=-0.5;cmax=0.5;
subplot(5,2,1);topoplot(delta_L_H_L_corr,channels,'nosedir','+X','style','map');title('delta-L & H-L');colorbar;colormap('jet');clim([cmin cmax]);
subplot(5,2,3);topoplot(theta_L_H_L_corr,channels,'nosedir','+X','style','map');title('theta-L & H-L');colorbar;colormap('jet');clim([cmin cmax]);
subplot(5,2,5);topoplot(alpha_L_H_L_corr,channels,'nosedir','+X','style','map');title('alpha-L & H-L');colorbar;colormap('jet');clim([cmin cmax]);
subplot(5,2,7);topoplot(beta_L_H_L_corr,channels,'nosedir','+X','style','map');title('beta-L & H-L');colorbar;colormap('jet');clim([cmin cmax]);
subplot(5,2,9);topoplot(gamma_L_H_L_corr,channels,'nosedir','+X','style','map');title('gamma-L & H-L');colorbar;colormap('jet');clim([cmin cmax]);
subplot(5,2,2);topoplot(delta_R_H_R_corr,channels,'nosedir','+X','style','map');title('delta-R & H-R');colorbar;colormap('jet');clim([cmin cmax]);
subplot(5,2,4);topoplot(theta_R_H_R_corr,channels,'nosedir','+X','style','map');title('theta-R & H-R');colorbar;colormap('jet');clim([cmin cmax]);
subplot(5,2,6);topoplot(alpha_R_H_R_corr,channels,'nosedir','+X','style','map');title('alpha-R & H-R');colorbar;colormap('jet');clim([cmin cmax]);
subplot(5,2,8);topoplot(beta_R_H_R_corr,channels,'nosedir','+X','style','map');title('beta-R & H-R');colorbar;colormap('jet');clim([cmin cmax]);
subplot(5,2,10);topoplot(gamma_R_H_R_corr,channels,'nosedir','+X','style','map');title('gamma-R & H-R');colorbar;colormap('jet');clim([cmin cmax]);
sgtitle({['corr of sum-EEG (-500ms) and H-int ^{* PLOT 13}'],char(datetime('now'))});
end

% Combine L and R in correlation;
for c=1:32
    delta_LR_H_LR_corr(c)=corr(delta_LR_chan(:,c),H_all_LR);
    theta_LR_H_LR_corr(c)=corr(theta_LR_chan(:,c),H_all_LR);
    alpha_LR_H_LR_corr(c)=corr(alpha_LR_chan(:,c),H_all_LR);
    beta_LR_H_LR_corr(c)=corr(beta_LR_chan(:,c),H_all_LR);
    gamma_LR_H_LR_corr(c)=corr(gamma_LR_chan(:,c),H_all_LR);
end
% Combine L and R in topoplots (1x5)
canvas(0.5,0.25);
for topo1x5=1;
cmin=-0.3;cmax=0.3;
subplot(1,5,1);topoplot(delta_LR_H_LR_corr,channels,'nosedir','+X','style','map');title('delta & H');colorbar;colormap('jet');clim([cmin cmax]);
subplot(1,5,2);topoplot(theta_LR_H_LR_corr,channels,'nosedir','+X','style','map');title('theta & H');colorbar;colormap('jet');clim([cmin cmax]);
subplot(1,5,3);topoplot(alpha_LR_H_LR_corr,channels,'nosedir','+X','style','map');title('alpha & H');colorbar;colormap('jet');clim([cmin cmax]);
subplot(1,5,4);topoplot(beta_LR_H_LR_corr,channels,'nosedir','+X','style','map');title('beta & H');colorbar;colormap('jet');clim([cmin cmax]);
subplot(1,5,5);topoplot(gamma_LR_H_LR_corr,channels,'nosedir','+X','style','map');title('gamma & H');colorbar;colormap('jet');clim([cmin cmax]);
sgtitle({['corr of sum-EEG (-500ms) and H-int ^{* PLOT 13}'],char(datetime('now'))});
colormap(hnc)
end

% Correlation in 4 states (combined synch/o)
delta_LR_H_LR_4corr=zeros(4,32);
theta_LR_H_LR_4corr=zeros(4,32);
alpha_LR_H_LR_4corr=zeros(4,32);
beta_LR_H_LR_4corr=zeros(4,32);
gamma_LR_H_LR_4corr=zeros(4,32);
for s=1:4
    for c=1:32
        delta_LR_H_LR_4corr(s,c)=corr(delta_LR_chan(Inds4_LR(:,s),c),H_all_LR(Inds4_LR(:,s)));
        theta_LR_H_LR_4corr(s,c)=corr(theta_LR_chan(Inds4_LR(:,s),c),H_all_LR(Inds4_LR(:,s)));
        alpha_LR_H_LR_4corr(s,c)=corr(alpha_LR_chan(Inds4_LR(:,s),c),H_all_LR(Inds4_LR(:,s)));
        beta_LR_H_LR_4corr(s,c)=corr(beta_LR_chan(Inds4_LR(:,s),c),H_all_LR(Inds4_LR(:,s)));
        gamma_LR_H_LR_4corr(s,c)=corr(gamma_LR_chan(Inds4_LR(:,s),c),H_all_LR(Inds4_LR(:,s)));
    end
end
% Combine L and R in 4 states(4x5)
canvas(0.3,0.5);
cmin=-0.7;cmax=0.7;
for s=1:4
    subplot(4,5,5*(s-1)+1);
    topoplot(delta_LR_H_LR_4corr(s,:),channels,'nosedir','+X','style','fill','numcontour',12);
    % title([states4names{s} ': delta & H'],'Color',condicolors(s,:));
    % colorbar;colormap('jet');
    clim([cmin cmax]);
    subplot(4,5,5*(s-1)+2);
    topoplot(theta_LR_H_LR_4corr(s,:),channels,'nosedir','+X','style','fill','numcontour',12);
    % title([states4names{s} ': theta & H'],'Color',condicolors(s,:));
    % colorbar;colormap('jet');
    clim([cmin cmax]);
    subplot(4,5,5*(s-1)+3);
    topoplot(alpha_LR_H_LR_4corr(s,:),channels,'nosedir','+X','style','fill','numcontour',12);
    % title([states4names{s} ': alpha & H'],'Color',condicolors(s,:));
    % colorbar;colormap('jet');
    clim([cmin cmax]);
    subplot(4,5,5*(s-1)+4);
    topoplot(beta_LR_H_LR_4corr(s,:),channels,'nosedir','+X','style','fill','numcontour',12);
    % title([states4names{s} ': beta & H'],'Color',condicolors(s,:));
    % colorbar;colormap('jet');
    clim([cmin cmax]);
    subplot(4,5,5*(s-1)+5);
    topoplot(gamma_LR_H_LR_4corr(s,:),channels,'nosedir','+X','style','fill','numcontour',12);
    % title([states4names{s} ': gamma & H'],'Color',condicolors(s,:));
    % colorbar;colormap('jet');
    clim([cmin cmax]);
end
% sgtitle('4states: Corr of sum-EEG (-500ms) and H-int ^{* PLOT 13}')
colormap(hotncold(12))
% https://www.mathworks.com/help/matlab/ref/matlab.graphics.illustration.colorbar-properties.html
cb=colorbar;
cb.AxisLocation = 'out';
cb.Position = [0.92 0.15 0.01 0.75];
% dim is [x position, y position of starting point, width, height] for
% i=1:5 h(i)=annotation('textbox',[(i-1)/5+15/100 0/5+19/20 2/100
% 1/100],'string','A','color',[0 0 0]); end for i=1:4
% v(i)=annotation('textbox',[0/5+2/20 (i-1)/5+1/20 2/100
% 1/100],'string','H','color',...
%     condicolors(i,:));
% end set([h v], 'fitboxtotext','on',...
%     'edgecolor','none')
% clf h(1)
delete(findall(gcf,'type','annotation'))
h0=annotation('textbox',[0.17 0.95 0.05 0.03],'string','Delta','color',[0 0 0])
h1=annotation('textbox',[0.33 0.95 0.05 0.03],'string','Theta','color',[0 0 0])
h2=annotation('textbox',[0.5 0.95 0.05 0.03],'string','Alpha','color',[0 0 0])
h3=annotation('textbox',[0.66 0.95 0.05 0.03],'string','Beta','color',[0 0 0])
h4=annotation('textbox',[0.81 0.95 0.05 0.03],'string','Gamma','color',[0 0 0])
v0=annotation('textbox',[0.14 0.15 0.05 0.03],'string','Mutual','color',condicolors(4,:))
v1=annotation('textbox',[0.14 0.37 0.05 0.03],'string','Following','color',condicolors(3,:))
v2=annotation('textbox',[0.14 0.59 0.05 0.03],'string','Leading','color',condicolors(2,:))
v3=annotation('textbox',[0.14 0.81 0.05 0.03],'string','Uncoupled','color',condicolors(1,:))
set(v0,'Rotation',90);set(v1,'Rotation',90);set(v2,'Rotation',90);set(v3,'Rotation',90);
delete(sg)
sg=annotation('textbox',[0.3 0.01 0.4 0.07],'string',...
    {['Correlation of sum-EEG (-500ms) and H-int ^{* PLOT 13}'],char(datetime('now'))});
set([h0 h1 h2 h3 h4 v0 v1 v2 v3], 'fitboxtotext','on',...
    'edgecolor','none')
set(gcf,'color','w'); % set background white for copying in ubuntu


% Correlation in 4 states (separate synch/o)
delta_LR_H_LR_4corr=zeros(2,4,32);% 2 sync types
theta_LR_H_LR_4corr=zeros(2,4,32);
alpha_LR_H_LR_4corr=zeros(2,4,32);
beta_LR_H_LR_4corr=zeros(2,4,32);
gamma_LR_H_LR_4corr=zeros(2,4,32);
for t=1:2
for s=1:4
    for c=1:32
        delta_LR_H_LR_4corr(t,s,c)=corr(delta_LR_chan(Inds4_LR([synind(t,:) 36+synind(t,:)],s),c),H_all_LR(Inds4_LR([synind(t,:) 36+synind(t,:)],s)));
        theta_LR_H_LR_4corr(t,s,c)=corr(theta_LR_chan(Inds4_LR([synind(t,:) 36+synind(t,:)],s),c),H_all_LR(Inds4_LR([synind(t,:) 36+synind(t,:)],s)));
        alpha_LR_H_LR_4corr(t,s,c)=corr(alpha_LR_chan(Inds4_LR([synind(t,:) 36+synind(t,:)],s),c),H_all_LR(Inds4_LR([synind(t,:) 36+synind(t,:)],s)));
        beta_LR_H_LR_4corr(t,s,c)=corr(beta_LR_chan(Inds4_LR([synind(t,:) 36+synind(t,:)],s),c),H_all_LR(Inds4_LR([synind(t,:) 36+synind(t,:)],s)));
        gamma_LR_H_LR_4corr(t,s,c)=corr(gamma_LR_chan(Inds4_LR([synind(t,:) 36+synind(t,:)],s),c),H_all_LR(Inds4_LR([synind(t,:) 36+synind(t,:)],s)));
    end
end
end
% Combine L and R in 4 states(4x5)
for t=1:2
canvas(0.3,0.5);
cmin=-0.7;cmax=0.7;
for s=1:4
    subplot(4,5,5*(s-1)+1);
    topoplot(delta_LR_H_LR_4corr(t,s,:),channels,'nosedir','+X','style','fill','numcontour',12);
    % title([states4names{s} ': delta & H'],'Color',condicolors(s,:));
    % colorbar;colormap('jet');
    clim([cmin cmax]);
    subplot(4,5,5*(s-1)+2);
    topoplot(theta_LR_H_LR_4corr(t,s,:),channels,'nosedir','+X','style','fill','numcontour',12);
    % title([states4names{s} ': theta & H'],'Color',condicolors(s,:));
    % colorbar;colormap('jet');
    clim([cmin cmax]);
    subplot(4,5,5*(s-1)+3);
    topoplot(alpha_LR_H_LR_4corr(t,s,:),channels,'nosedir','+X','style','fill','numcontour',12);
    % title([states4names{s} ': alpha & H'],'Color',condicolors(s,:));
    % colorbar;colormap('jet');
    clim([cmin cmax]);
    subplot(4,5,5*(s-1)+4);
    topoplot(beta_LR_H_LR_4corr(t,s,:),channels,'nosedir','+X','style','fill','numcontour',12);
    % title([states4names{s} ': beta & H'],'Color',condicolors(s,:));
    % colorbar;colormap('jet');
    clim([cmin cmax]);
    subplot(4,5,5*(s-1)+5);
    topoplot(gamma_LR_H_LR_4corr(t,s,:),channels,'nosedir','+X','style','fill','numcontour',12);
    % title([states4names{s} ': gamma & H'],'Color',condicolors(s,:));
    % colorbar;colormap('jet');
    clim([cmin cmax]);
end
% sgtitle('4states: Corr of sum-EEG (-500ms) and H-int ^{* PLOT 13}')
colormap(hotncold(12))
% https://www.mathworks.com/help/matlab/ref/matlab.graphics.illustration.colorbar-properties.html
cb=colorbar;
cb.AxisLocation = 'out';
cb.Position = [0.92 0.15 0.01 0.75];
% dim is [x position, y position of starting point, width, height] for
% i=1:5 h(i)=annotation('textbox',[(i-1)/5+15/100 0/5+19/20 2/100
% 1/100],'string','A','color',[0 0 0]); end for i=1:4
% v(i)=annotation('textbox',[0/5+2/20 (i-1)/5+1/20 2/100
% 1/100],'string','H','color',...
%     condicolors(i,:));
% end set([h v], 'fitboxtotext','on',...
%     'edgecolor','none')
% clf h(1)
delete(findall(gcf,'type','annotation'))
h0=annotation('textbox',[0.17 0.95 0.05 0.03],'string','Delta','color',[0 0 0])
h1=annotation('textbox',[0.33 0.95 0.05 0.03],'string','Theta','color',[0 0 0])
h2=annotation('textbox',[0.5 0.95 0.05 0.03],'string','Alpha','color',[0 0 0])
h3=annotation('textbox',[0.66 0.95 0.05 0.03],'string','Beta','color',[0 0 0])
h4=annotation('textbox',[0.81 0.95 0.05 0.03],'string','Gamma','color',[0 0 0])
v0=annotation('textbox',[0.14 0.15 0.05 0.03],'string','Mutual','color',condicolors(4,:))
v1=annotation('textbox',[0.14 0.37 0.05 0.03],'string','Following','color',condicolors(3,:))
v2=annotation('textbox',[0.14 0.59 0.05 0.03],'string','Leading','color',condicolors(2,:))
v3=annotation('textbox',[0.14 0.81 0.05 0.03],'string','Uncoupled','color',condicolors(1,:))
set(v0,'Rotation',90);set(v1,'Rotation',90);set(v2,'Rotation',90);set(v3,'Rotation',90);
% delete(sg)
if t==1;
sg=annotation('textbox',[0.3 0.01 0.4 0.07],'string',...
    {['Correlation of sum-EEG (-500ms) and H-int in synch ^{* PLOT 13}'],char(datetime('now'))},'Color',syn2colors(1,:));
else t==2;
sg=annotation('textbox',[0.3 0.01 0.4 0.07],'string',...
    {['Correlation of sum-EEG (-500ms) and H-int in synco^{* PLOT 13}'],char(datetime('now'))},'Color',syn2colors(2,:));
end
set([h0 h1 h2 h3 h4 v0 v1 v2 v3], 'fitboxtotext','on',...
    'edgecolor','none')
set(gcf,'color','w'); % set background white for copying in ubuntu
end

%% PLOT 13-0 corr of sum-EEG power (-500ms) and H-diff all sessions
% Correlation in 4 states  (combined synch/o)
delta_LR_H_LR_4corr=zeros(4,32);
theta_LR_H_LR_4corr=zeros(4,32);
alpha_LR_H_LR_4corr=zeros(4,32);
beta_LR_H_LR_4corr=zeros(4,32);
gamma_LR_H_LR_4corr=zeros(4,32);
for s=1:4
    for c=1:32
        delta_LR_H_LR_4corr(s,c)=corr(delta_LR_chan(Inds4_LR(:,s),c),H_all_diff_LR(Inds4_LR(:,s))); % from SECT 13 & PLOT 16-0
        theta_LR_H_LR_4corr(s,c)=corr(theta_LR_chan(Inds4_LR(:,s),c),H_all_diff_LR(Inds4_LR(:,s)));
        alpha_LR_H_LR_4corr(s,c)=corr(alpha_LR_chan(Inds4_LR(:,s),c),H_all_diff_LR(Inds4_LR(:,s)));
        beta_LR_H_LR_4corr(s,c)=corr(beta_LR_chan(Inds4_LR(:,s),c),H_all_diff_LR(Inds4_LR(:,s)));
        gamma_LR_H_LR_4corr(s,c)=corr(gamma_LR_chan(Inds4_LR(:,s),c),H_all_diff_LR(Inds4_LR(:,s)));
    end
end
% Combine L and R in 4 states(4x5)
canvas(0.3,0.5);
cmin=-0.7;cmax=0.7;
for s=1:4
    subplot(4,5,5*(s-1)+1);
    topoplot(delta_LR_H_LR_4corr(s,:),channels,'nosedir','+X','style','fill','numcontour',12);
    % title([states4names{s} ': delta & H'],'Color',condicolors(s,:));
    % colorbar;colormap('jet');
    clim([cmin cmax]);
    subplot(4,5,5*(s-1)+2);
    topoplot(theta_LR_H_LR_4corr(s,:),channels,'nosedir','+X','style','fill','numcontour',12);
    % title([states4names{s} ': theta & H'],'Color',condicolors(s,:));
    % colorbar;colormap('jet');
    clim([cmin cmax]);
    subplot(4,5,5*(s-1)+3);
    topoplot(alpha_LR_H_LR_4corr(s,:),channels,'nosedir','+X','style','fill','numcontour',12);
    % title([states4names{s} ': alpha & H'],'Color',condicolors(s,:));
    % colorbar;colormap('jet');
    clim([cmin cmax]);
    subplot(4,5,5*(s-1)+4);
    topoplot(beta_LR_H_LR_4corr(s,:),channels,'nosedir','+X','style','fill','numcontour',12);
    % title([states4names{s} ': beta & H'],'Color',condicolors(s,:));
    % colorbar;colormap('jet');
    clim([cmin cmax]);
    subplot(4,5,5*(s-1)+5);
    topoplot(gamma_LR_H_LR_4corr(s,:),channels,'nosedir','+X','style','fill','numcontour',12);
    % title([states4names{s} ': gamma & H'],'Color',condicolors(s,:));
    % colorbar;colormap('jet');
    clim([cmin cmax]);
end
% sgtitle('4states: Corr of sum-EEG (-500ms) and H-int ^{* PLOT 13}')
colormap(hotncold(12))
if true; %annotations
% https://www.mathworks.com/help/matlab/ref/matlab.graphics.illustration.colorbar-properties.html
cb=colorbar;
cb.AxisLocation = 'out';
cb.Position = [0.92 0.15 0.01 0.75];
% dim is [x position, y position of starting point, width, height] for
% i=1:5 h(i)=annotation('textbox',[(i-1)/5+15/100 0/5+19/20 2/100
% 1/100],'string','A','color',[0 0 0]); end for i=1:4
% v(i)=annotation('textbox',[0/5+2/20 (i-1)/5+1/20 2/100
% 1/100],'string','H','color',...
%     condicolors(i,:));
% end set([h v], 'fitboxtotext','on',...
%     'edgecolor','none')
% clf h(1)
delete(findall(gcf,'type','annotation'))
h0=annotation('textbox',[0.17 0.95 0.05 0.03],'string','Delta','color',[0 0 0])
h1=annotation('textbox',[0.33 0.95 0.05 0.03],'string','Theta','color',[0 0 0])
h2=annotation('textbox',[0.5 0.95 0.05 0.03],'string','Alpha','color',[0 0 0])
h3=annotation('textbox',[0.66 0.95 0.05 0.03],'string','Beta','color',[0 0 0])
h4=annotation('textbox',[0.81 0.95 0.05 0.03],'string','Gamma','color',[0 0 0])
v0=annotation('textbox',[0.14 0.15 0.05 0.03],'string','Mutual','color',condicolors(4,:))
v1=annotation('textbox',[0.14 0.37 0.05 0.03],'string','Following','color',condicolors(3,:))
v2=annotation('textbox',[0.14 0.59 0.05 0.03],'string','Leading','color',condicolors(2,:))
v3=annotation('textbox',[0.14 0.81 0.05 0.03],'string','Uncoupled','color',condicolors(1,:))
set(v0,'Rotation',90);set(v1,'Rotation',90);set(v2,'Rotation',90);set(v3,'Rotation',90);
delete(sg)
sg=annotation('textbox',[0.3 0.01 0.4 0.07],'string',...
    {['Correlation of sum-EEG (-500ms) and H-diff ^{* PLOT 13}'],char(datetime('now'))});
set([h0 h1 h2 h3 h4 v0 v1 v2 v3], 'fitboxtotext','on',...
    'edgecolor','none')
set(gcf,'color','w'); % set background white for copying in ubuntu
end

% Correlation in 4 states (separate synch/o)
delta_LR_H_LR_4corr=zeros(2,4,32);% 2 sync types
theta_LR_H_LR_4corr=zeros(2,4,32);
alpha_LR_H_LR_4corr=zeros(2,4,32);
beta_LR_H_LR_4corr=zeros(2,4,32);
gamma_LR_H_LR_4corr=zeros(2,4,32);
for t=1:2
for s=1:4
    for c=1:32
        delta_LR_H_LR_4corr(t,s,c)=corr(delta_LR_chan(Inds4_LR([synind(t,:) 36+synind(t,:)],s),c),H_all_diff_LR(Inds4_LR([synind(t,:) 36+synind(t,:)],s)));
        theta_LR_H_LR_4corr(t,s,c)=corr(theta_LR_chan(Inds4_LR([synind(t,:) 36+synind(t,:)],s),c),H_all_diff_LR(Inds4_LR([synind(t,:) 36+synind(t,:)],s)));
        alpha_LR_H_LR_4corr(t,s,c)=corr(alpha_LR_chan(Inds4_LR([synind(t,:) 36+synind(t,:)],s),c),H_all_diff_LR(Inds4_LR([synind(t,:) 36+synind(t,:)],s)));
        beta_LR_H_LR_4corr(t,s,c)=corr(beta_LR_chan(Inds4_LR([synind(t,:) 36+synind(t,:)],s),c),H_all_diff_LR(Inds4_LR([synind(t,:) 36+synind(t,:)],s)));
        gamma_LR_H_LR_4corr(t,s,c)=corr(gamma_LR_chan(Inds4_LR([synind(t,:) 36+synind(t,:)],s),c),H_all_diff_LR(Inds4_LR([synind(t,:) 36+synind(t,:)],s)));
    end
end
end
% Combine L and R in 4 states(4x5)
for t=1:2
canvas(0.3,0.5);
cmin=-0.7;cmax=0.7;
for s=1:4
    subplot(4,5,5*(s-1)+1);
    topoplot(delta_LR_H_LR_4corr(t,s,:),channels,'nosedir','+X','style','fill','numcontour',12);
    % title([states4names{s} ': delta & H'],'Color',condicolors(s,:));
    % colorbar;colormap('jet');
    clim([cmin cmax]);
    subplot(4,5,5*(s-1)+2);
    topoplot(theta_LR_H_LR_4corr(t,s,:),channels,'nosedir','+X','style','fill','numcontour',12);
    % title([states4names{s} ': theta & H'],'Color',condicolors(s,:));
    % colorbar;colormap('jet');
    clim([cmin cmax]);
    subplot(4,5,5*(s-1)+3);
    topoplot(alpha_LR_H_LR_4corr(t,s,:),channels,'nosedir','+X','style','fill','numcontour',12);
    % title([states4names{s} ': alpha & H'],'Color',condicolors(s,:));
    % colorbar;colormap('jet');
    clim([cmin cmax]);
    subplot(4,5,5*(s-1)+4);
    topoplot(beta_LR_H_LR_4corr(t,s,:),channels,'nosedir','+X','style','fill','numcontour',12);
    % title([states4names{s} ': beta & H'],'Color',condicolors(s,:));
    % colorbar;colormap('jet');
    clim([cmin cmax]);
    subplot(4,5,5*(s-1)+5);
    topoplot(gamma_LR_H_LR_4corr(t,s,:),channels,'nosedir','+X','style','fill','numcontour',12);
    % title([states4names{s} ': gamma & H'],'Color',condicolors(s,:));
    % colorbar;colormap('jet');
    clim([cmin cmax]);
end
% sgtitle('4states: Corr of sum-EEG (-500ms) and H-int ^{* PLOT 13}')
colormap(hotncold(12))
% https://www.mathworks.com/help/matlab/ref/matlab.graphics.illustration.colorbar-properties.html
cb=colorbar;
cb.AxisLocation = 'out';
cb.Position = [0.92 0.15 0.01 0.75];
% dim is [x position, y position of starting point, width, height] for
% i=1:5 h(i)=annotation('textbox',[(i-1)/5+15/100 0/5+19/20 2/100
% 1/100],'string','A','color',[0 0 0]); end for i=1:4
% v(i)=annotation('textbox',[0/5+2/20 (i-1)/5+1/20 2/100
% 1/100],'string','H','color',...
%     condicolors(i,:));
% end set([h v], 'fitboxtotext','on',...
%     'edgecolor','none')
% clf h(1)
delete(findall(gcf,'type','annotation'))
h0=annotation('textbox',[0.17 0.95 0.05 0.03],'string','Delta','color',[0 0 0])
h1=annotation('textbox',[0.33 0.95 0.05 0.03],'string','Theta','color',[0 0 0])
h2=annotation('textbox',[0.5 0.95 0.05 0.03],'string','Alpha','color',[0 0 0])
h3=annotation('textbox',[0.66 0.95 0.05 0.03],'string','Beta','color',[0 0 0])
h4=annotation('textbox',[0.81 0.95 0.05 0.03],'string','Gamma','color',[0 0 0])
v0=annotation('textbox',[0.14 0.15 0.05 0.03],'string','Mutual','color',condicolors(4,:))
v1=annotation('textbox',[0.14 0.37 0.05 0.03],'string','Following','color',condicolors(3,:))
v2=annotation('textbox',[0.14 0.59 0.05 0.03],'string','Leading','color',condicolors(2,:))
v3=annotation('textbox',[0.14 0.81 0.05 0.03],'string','Uncoupled','color',condicolors(1,:))
set(v0,'Rotation',90);set(v1,'Rotation',90);set(v2,'Rotation',90);set(v3,'Rotation',90);
% delete(sg)
if t==1;
sg=annotation('textbox',[0.3 0.01 0.4 0.07],'string',...
    {['Correlation of sum-EEG (-500ms) and H-diff in synch ^{* PLOT 13-0}'],char(datetime('now'))},'Color',syn2colors(1,:));
else t==2;
sg=annotation('textbox',[0.3 0.01 0.4 0.07],'string',...
    {['Correlation of sum-EEG (-500ms) and H-diff in synco^{* PLOT 13-0}'],char(datetime('now'))},'Color',syn2colors(2,:));
end
set([h0 h1 h2 h3 h4 v0 v1 v2 v3], 'fitboxtotext','on',...
    'edgecolor','none')
set(gcf,'color','w'); % set background white for copying in ubuntu
end

%% PLOT 13-1 corr of sum-EEG power (+ 500ms) and H-int all sessions
zEEG500_LL; zEEG500_RR; % for all sessions from SECT 12
delta_LL;theta_LL;alpha_LL;beta_LL;gamma_LL;
delta_RR;theta_RR;alpha_RR;beta_RR;gamma_RR;
% sum the band power in each channel in each block (8x12=96)
delta_LL_sum=[];theta_LL_sum=[];alpha_LL_sum=[];beta_LL_sum=[];gamma_LL_sum=[];
delta_RR_sum=[];theta_RR_sum=[];alpha_RR_sum=[];beta_RR_sum=[];gamma_RR_sum=[];
for s=1:numSes
    for b=1:12
        delta_LL_sum(s,b,:)=sum([delta_LL{s,b}]);
        theta_LL_sum(s,b,:)=sum([theta_LL{s,b}]);
        alpha_LL_sum(s,b,:)=sum([alpha_LL{s,b}]);
        beta_LL_sum(s,b,:)=sum([beta_LL{s,b}]);
        gamma_LL_sum(s,b,:)=sum([gamma_LL{s,b}]);
        delta_RR_sum(s,b,:)=sum([delta_RR{s,b}]);
        theta_RR_sum(s,b,:)=sum([theta_RR{s,b}]);
        alpha_RR_sum(s,b,:)=sum([alpha_RR{s,b}]);
        beta_RR_sum(s,b,:)=sum([beta_RR{s,b}]);
        gamma_RR_sum(s,b,:)=sum([gamma_RR{s,b}]);
    end
end
% squeeze into 32 vectors for corr with H
delta_LL_chan=[];theta_LL_chan=[];alpha_LL_chan=[];beta_LL_chan=[];gamma_LL_chan=[];
delta_RR_chan=[];theta_RR_chan=[];alpha_RR_chan=[];beta_RR_chan=[];gamma_RR_chan=[];
for c=1:32
    delta_LL_chan(:,c)=reshape(delta_LL_sum(:,:,c)',[],1); % 144x32 (each element from one block in time sequence) 
    theta_LL_chan(:,c)=reshape(theta_LL_sum(:,:,c)',[],1);
    alpha_LL_chan(:,c)=reshape(alpha_LL_sum(:,:,c)',[],1);
    beta_LL_chan(:,c)=reshape(beta_LL_sum(:,:,c)',[],1);
    gamma_LL_chan(:,c)=reshape(gamma_LL_sum(:,:,c)',[],1);
    delta_RR_chan(:,c)=reshape(delta_RR_sum(:,:,c)',[],1);
    theta_RR_chan(:,c)=reshape(theta_RR_sum(:,:,c)',[],1);
    alpha_RR_chan(:,c)=reshape(alpha_RR_sum(:,:,c)',[],1);
    beta_RR_chan(:,c)=reshape(beta_RR_sum(:,:,c)',[],1);
    gamma_RR_chan(:,c)=reshape(gamma_RR_sum(:,:,c)',[],1);
end
% Combine L and R
delta_LLRR_chan=[delta_LL_chan;delta_RR_chan];
theta_LLRR_chan=[theta_LL_chan;theta_RR_chan];
alpha_LLRR_chan=[alpha_LL_chan;alpha_RR_chan];
beta_LLRR_chan=[beta_LL_chan;beta_RR_chan];
gamma_LLRR_chan=[gamma_LL_chan;gamma_RR_chan];

% % Organize H_all for corr (same as in PLOT 13)
% H_all; % (2x 8x12) for all sessions from SECT 10-1 
% H_all_L=squeeze(H_all(1,:,:));
% H_all_R=squeeze(H_all(2,:,:));
% % squeeze into 1 vector from the 96 blocks for each subject, for corr with pow in each chan
% H_all_L=reshape(H_all_L',[],1);% 96x1 (each element from one block in time sequence) 
% H_all_R=reshape(H_all_R',[],1);
% % Combine L and R
% H_all_LR=[H_all_L;H_all_R];

% Compute the correlation between sum-EEG pow and H-int
for c=1:32
    delta_LL_H_L_corr(c)=corr(delta_LL_chan(:,c),H_all_L);
    theta_LL_H_L_corr(c)=corr(theta_LL_chan(:,c),H_all_L);
    alpha_LL_H_L_corr(c)=corr(alpha_LL_chan(:,c),H_all_L);
    beta_LL_H_L_corr(c)=corr(beta_LL_chan(:,c),H_all_L);
    gamma_LL_H_L_corr(c)=corr(gamma_LL_chan(:,c),H_all_L);
    delta_RR_H_R_corr(c)=corr(delta_RR_chan(:,c),H_all_R);
    theta_RR_H_R_corr(c)=corr(theta_RR_chan(:,c),H_all_R);
    alpha_RR_H_R_corr(c)=corr(alpha_RR_chan(:,c),H_all_R);
    beta_RR_H_R_corr(c)=corr(beta_RR_chan(:,c),H_all_R);
    gamma_RR_H_R_corr(c)=corr(gamma_RR_chan(:,c),H_all_R);
end
% Topoplots (5x2)
canvas(0.2,0.5);
cmin=-0.1;cmax=0.3;
subplot(5,2,1);topoplot(delta_LL_H_L_corr,channels,'nosedir','+X','style','map');title('delta-L & H-L');colorbar;colormap('jet');clim([cmin cmax]);
subplot(5,2,3);topoplot(theta_LL_H_L_corr,channels,'nosedir','+X','style','map');title('theta-L & H-L');colorbar;colormap('jet');clim([cmin cmax]);
subplot(5,2,5);topoplot(alpha_LL_H_L_corr,channels,'nosedir','+X','style','map');title('alpha-L & H-L');colorbar;colormap('jet');clim([cmin cmax]);
subplot(5,2,7);topoplot(beta_LL_H_L_corr,channels,'nosedir','+X','style','map');title('beta-L & H-L');colorbar;colormap('jet');clim([cmin cmax]);
subplot(5,2,9);topoplot(gamma_LL_H_L_corr,channels,'nosedir','+X','style','map');title('gamma-L & H-L');colorbar;colormap('jet');clim([cmin cmax]);
subplot(5,2,2);topoplot(delta_RR_H_R_corr,channels,'nosedir','+X','style','map');title('delta-R & H-R');colorbar;colormap('jet');clim([cmin cmax]);
subplot(5,2,4);topoplot(theta_RR_H_R_corr,channels,'nosedir','+X','style','map');title('theta-R & H-R');colorbar;colormap('jet');clim([cmin cmax]);
subplot(5,2,6);topoplot(alpha_RR_H_R_corr,channels,'nosedir','+X','style','map');title('alpha-R & H-R');colorbar;colormap('jet');clim([cmin cmax]);
subplot(5,2,8);topoplot(beta_RR_H_R_corr,channels,'nosedir','+X','style','map');title('beta-R & H-R');colorbar;colormap('jet');clim([cmin cmax]);
subplot(5,2,10);topoplot(gamma_RR_H_R_corr,channels,'nosedir','+X','style','map');title('gamma-R & H-R');colorbar;colormap('jet');clim([cmin cmax]);
sgtitle({['corr of sum-EEG (+500ms) and H-int'],char(datetime('now'))});

% Combine L and R in correlation;
for c=1:32
    delta_LLRR_H_LR_corr(c)=corr(delta_LLRR_chan(:,c),H_all_LR);
    theta_LLRR_H_LR_corr(c)=corr(theta_LLRR_chan(:,c),H_all_LR);
    alpha_LLRR_H_LR_corr(c)=corr(alpha_LLRR_chan(:,c),H_all_LR);
    beta_LLRR_H_LR_corr(c)=corr(beta_LLRR_chan(:,c),H_all_LR);
    gamma_LLRR_H_LR_corr(c)=corr(gamma_LLRR_chan(:,c),H_all_LR);
end
% Combine L and R in topoplots (1x5)
canvas(0.5,0.1);
cmin=-0.1;cmax=0.3;
subplot(1,5,1);topoplot(delta_LLRR_H_LR_corr,channels,'nosedir','+X','style','map');title('delta & H');colorbar;colormap('jet');clim([cmin cmax]);
subplot(1,5,2);topoplot(theta_LLRR_H_LR_corr,channels,'nosedir','+X','style','map');title('theta & H');colorbar;colormap('jet');clim([cmin cmax]);
subplot(1,5,3);topoplot(alpha_LLRR_H_LR_corr,channels,'nosedir','+X','style','map');title('alpha & H');colorbar;colormap('jet');clim([cmin cmax]);
subplot(1,5,4);topoplot(beta_LLRR_H_LR_corr,channels,'nosedir','+X','style','map');title('beta & H');colorbar;colormap('jet');clim([cmin cmax]);
subplot(1,5,5);topoplot(gamma_LLRR_H_LR_corr,channels,'nosedir','+X','style','map');title('gamma & H');colorbar;colormap('jet');clim([cmin cmax]);
sgtitle({['corr of sum-EEG (+500ms) and H-int ^{* PLOT 13-1}'],char(datetime('now'))});
colormap(hnc)

% Correlation in 4 states
delta_LLRR_H_LR_4corr=zeros(4,32);
theta_LLRR_H_LR_4corr=zeros(4,32);
alpha_LLRR_H_LR_4corr=zeros(4,32);
beta_LLRR_H_LR_4corr=zeros(4,32);
gamma_LLRR_H_LR_4corr=zeros(4,32);
for s=1:4
    for c=1:32
        delta_LLRR_H_LR_4corr(s,c)=corr(delta_LLRR_chan(Inds4_LR(:,s),c),H_all_LR(Inds4_LR(:,s)));
        theta_LLRR_H_LR_4corr(s,c)=corr(theta_LLRR_chan(Inds4_LR(:,s),c),H_all_LR(Inds4_LR(:,s)));
        alpha_LLRR_H_LR_4corr(s,c)=corr(alpha_LLRR_chan(Inds4_LR(:,s),c),H_all_LR(Inds4_LR(:,s)));
        beta_LLRR_H_LR_4corr(s,c)=corr(beta_LLRR_chan(Inds4_LR(:,s),c),H_all_LR(Inds4_LR(:,s)));
        gamma_LLRR_H_LR_4corr(s,c)=corr(gamma_LLRR_chan(Inds4_LR(:,s),c),H_all_LR(Inds4_LR(:,s)));
    end
end
% Combine L and R in 4 states(4x5)
canvas(0.5,0.6);
cmin=-0.2;cmax=0.4;
for s=1:4
    subplot(4,5,5*(s-1)+1);
    topoplot(delta_LLRR_H_LR_4corr(s,:),channels,'nosedir','+X','style','map');
    title([states4names{s} ': delta & H'],'Color',condicolors(s,:));
    colorbar;colormap('jet');clim([cmin cmax]);
    subplot(4,5,5*(s-1)+2);
    topoplot(theta_LLRR_H_LR_4corr(s,:),channels,'nosedir','+X','style','map');
    title([states4names{s} ': theta & H'],'Color',condicolors(s,:));
    colorbar;colormap('jet');clim([cmin cmax]);
    subplot(4,5,5*(s-1)+3);
    topoplot(alpha_LLRR_H_LR_4corr(s,:),channels,'nosedir','+X','style','map');
    title([states4names{s} ': alpha & H'],'Color',condicolors(s,:));
    colorbar;colormap('jet');clim([cmin cmax]);
    subplot(4,5,5*(s-1)+4);
    topoplot(beta_LLRR_H_LR_4corr(s,:),channels,'nosedir','+X','style','map');
    title([states4names{s} ': beta & H'],'Color',condicolors(s,:));
    colorbar;colormap('jet');clim([cmin cmax]);
    subplot(4,5,5*(s-1)+5);
    topoplot(gamma_LLRR_H_LR_4corr(s,:),channels,'nosedir','+X','style','map');
    title([states4names{s} ': gamma & H'],'Color',condicolors(s,:));
    colorbar;colormap('jet');clim([cmin cmax]);
end
sgtitle({['4states: Corr of sum-EEG (+500ms) and H-int'],char(datetime('now'))});
colormap(hnc)
%% PLOT 14 Compute H-EEG (- 500ms) all sessions & 4 states
% H derived from EEG power
delta_L;theta_L;alpha_L;beta_L;gamma_L;
delta_R;theta_R;alpha_R;beta_R;gamma_R;
H_delta_L=[];H_theat_L=[];H_alpha_L=[];H_beta_L=[];H_gamma_L=[];
H_delta_R=[];H_theat_R=[];H_alpha_R=[];H_beta_R=[];H_gamma_R=[];
for s=1:numSes
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
canvas(0.5,0.2);
cmin=0.3;cmax=0.7;
subplot(1,5,1);topoplot(mean(H_delta_LR_chan),channels,'nosedir','+X','style','map');title('H-delta');colorbar;colormap('jet');clim([cmin cmax]);
subplot(1,5,2);topoplot(mean(H_theta_LR_chan),channels,'nosedir','+X','style','map');title('H-theta');colorbar;colormap('jet');clim([cmin cmax]);
subplot(1,5,3);topoplot(mean(H_alpha_LR_chan),channels,'nosedir','+X','style','map');title('H-alpha');colorbar;colormap('jet');clim([cmin cmax]);
subplot(1,5,4);topoplot(mean(H_beta_LR_chan),channels,'nosedir','+X','style','map');title('H-beta');colorbar;colormap('jet');clim([cmin cmax]);
subplot(1,5,5);topoplot(mean(H_gamma_LR_chan),channels,'nosedir','+X','style','map');title('H-gamma');colorbar;colormap('jet');clim([cmin cmax]);
sgtitle({['mean H-EEG (- 500ms) ^{* PLOT 14}'],char(datetime('now'))});
colormap(hnc)
% only H-EEG of beta and gamma show complexity (H>0.5)

% check how many H-EEG > 0.5 according to 4 states (showing complexity)
% plot in 4 states
canvas(0.5,0.5);
cmin=0.3;cmax=0.7;
for c=1:4
subplot(4,5,(c-1)*5+1);
topoplot(mean(H_delta_LR_chan(Inds4_LR(:,c),:)),channels,'nosedir','+X','style','map');
title(['H-delta ' states4names{c}],'Color',condicolors(c,:));
colorbar;colormap('jet');clim([cmin cmax]);
subplot(4,5,(c-1)*5+2);
topoplot(mean(H_theta_LR_chan(Inds4_LR(:,c),:)),channels,'nosedir','+X','style','map');
title(['H-theta ' states4names{c}],'Color',condicolors(c,:));
colorbar;colormap('jet');clim([cmin cmax]);
subplot(4,5,(c-1)*5+3);
topoplot(mean(H_alpha_LR_chan(Inds4_LR(:,c),:)),channels,'nosedir','+X','style','map');
title(['H-alpha ' states4names{c}],'Color',condicolors(c,:));
colorbar;colormap('jet');clim([cmin cmax]);
subplot(4,5,(c-1)*5+4);
topoplot(mean(H_beta_LR_chan(Inds4_LR(:,c),:)),channels,'nosedir','+X','style','map');
title(['H-beta ' states4names{c}],'Color',condicolors(c,:));
colorbar;colormap('jet');clim([cmin cmax]);
subplot(4,5,(c-1)*5+5);
topoplot(mean(H_gamma_LR_chan(Inds4_LR(:,c),:)),channels,'nosedir','+X','style','map');
title(['H-gamma ' states4names{c}],'Color',condicolors(c,:));
colorbar;colormap('jet');clim([cmin cmax]);
end
sgtitle({['mean H-EEG (- 500ms) for each condition ^{* PLOT 14}'],char(datetime('now'))});
colormap(hnc)
% only H-EEG of beta and gamma show complexity (H>0.5)

%% PLOT 14-1 Compute H-EEG (+ 500ms) all sessions & 4 states
% H derived from EEG power
delta_LL;theta_LL;alpha_LL;beta_LL;gamma_LL;
delta_RR;theta_RR;alpha_RR;beta_RR;gamma_RR;
H_delta_LL=[];H_theat_LL=[];H_alpha_LL=[];H_beta_LL=[];H_gamma_LL=[];
H_delta_RR=[];H_theat_RR=[];H_alpha_RR=[];H_beta_RR=[];H_gamma_RR=[];
for s=1:numSes
    for b=1:12
        for c=1:32
            [~,H_delta_LL(s,b,c)]=DFA_main(delta_LL{s,b}(:,c));
            [~,H_theta_LL(s,b,c)]=DFA_main(theta_LL{s,b}(:,c));
            [~,H_alpha_LL(s,b,c)]=DFA_main(alpha_LL{s,b}(:,c));
            [~,H_beta_LL(s,b,c)]=DFA_main(beta_LL{s,b}(:,c));
            [~,H_gamma_LL(s,b,c)]=DFA_main(gamma_LL{s,b}(:,c));
            [~,H_delta_RR(s,b,c)]=DFA_main(delta_RR{s,b}(:,c));
            [~,H_theta_RR(s,b,c)]=DFA_main(theta_RR{s,b}(:,c));
            [~,H_alpha_RR(s,b,c)]=DFA_main(alpha_RR{s,b}(:,c));
            [~,H_beta_RR(s,b,c)]=DFA_main(beta_RR{s,b}(:,c));
            [~,H_gamma_RR(s,b,c)]=DFA_main(gamma_RR{s,b}(:,c)); 
        end
    end
end
% squeeze into 32 vectors for corr with H
H_delta_LL_chan=[];H_theat_LL_chan=[];H_alpha_LL_chan=[];H_beta_LL_chan=[];H_gamma_LL_chan=[];
H_delta_RR_chan=[];H_theat_RR_chan=[];H_alpha_RR_chan=[];H_beta_RR_chan=[];H_gamma_RR_chan=[];
for c=1:32
    H_delta_LL_chan(:,c) = reshape(H_delta_LL(:,:,c)',[],1);
    H_theta_LL_chan(:,c)  = reshape(H_theta_LL(:,:,c)',[],1);
    H_alpha_LL_chan(:,c)  = reshape(H_alpha_LL(:,:,c)',[],1);
    H_beta_LL_chan(:,c)   = reshape(H_beta_LL(:,:,c)',[],1);
    H_gamma_LL_chan(:,c)  = reshape(H_gamma_LL(:,:,c)',[],1);
    H_delta_RR_chan(:,c)  = reshape(H_delta_RR(:,:,c)',[],1);
    H_theta_RR_chan(:,c)  = reshape(H_theta_RR(:,:,c)',[],1);
    H_alpha_RR_chan(:,c)  = reshape(H_alpha_RR(:,:,c)',[],1);
    H_beta_RR_chan(:,c)   = reshape(H_beta_RR(:,:,c)',[],1);
    H_gamma_RR_chan(:,c)   = reshape(H_gamma_RR(:,:,c)',[],1);
end
% combine L and R
H_delta_LLRR_chan=[H_delta_LL_chan;H_delta_RR_chan];
H_theta_LLRR_chan=[H_theta_LL_chan;H_theta_RR_chan];
H_alpha_LLRR_chan=[H_alpha_LL_chan;H_alpha_RR_chan];
H_beta_LLRR_chan=[H_beta_LL_chan;H_beta_RR_chan];
H_gamma_LLRR_chan=[H_gamma_LL_chan;H_gamma_RR_chan];

% check how many H-EEG > 0.5 (showing complexity)
% mean(H_delta_LLRR_chan)
% mean(H_theta_LLRR_chan)
% mean(H_alpha_LLRR_chan)
% mean(H_beta_LLRR_chan)
% mean(H_gamma_LLRR_chan)
canvas(0.5,0.2);
cmin=0.3;cmax=0.7;
subplot(1,5,1);topoplot(mean(H_delta_LLRR_chan),channels,'nosedir','+X','style','map');title('H-delta');colorbar;colormap('jet');clim([cmin cmax]);
subplot(1,5,2);topoplot(mean(H_theta_LLRR_chan),channels,'nosedir','+X','style','map');title('H-theta');colorbar;colormap('jet');clim([cmin cmax]);
subplot(1,5,3);topoplot(mean(H_alpha_LLRR_chan),channels,'nosedir','+X','style','map');title('H-alpha');colorbar;colormap('jet');clim([cmin cmax]);
subplot(1,5,4);topoplot(mean(H_beta_LLRR_chan),channels,'nosedir','+X','style','map');title('H-beta');colorbar;colormap('jet');clim([cmin cmax]);
subplot(1,5,5);topoplot(mean(H_gamma_LLRR_chan),channels,'nosedir','+X','style','map');title('H-gamma');colorbar;colormap('jet');clim([cmin cmax]);
sgtitle({['mean H-EEG (+ 500ms) ^{* PLOT 14-1}'],char(datetime('now'))});
colormap(hnc)
% only H-EEG of beta and gamma show complexity (H>0.5)

% check how many H-EEG > 0.5 according to 4 states (showing complexity)
% plot in 4 states
canvas(0.5,0.5);
cmin=0.3;cmax=0.7;
for c=1:4
subplot(4,5,(c-1)*5+1);
topoplot(mean(H_delta_LLRR_chan(Inds4_LR(:,c),:)),channels,'nosedir','+X','style','map');
title(['H-delta ' states4names{c}],'Color',condicolors(c,:));
colorbar;colormap('jet');clim([cmin cmax]);
subplot(4,5,(c-1)*5+2);
topoplot(mean(H_theta_LLRR_chan(Inds4_LR(:,c),:)),channels,'nosedir','+X','style','map');
title(['H-theta ' states4names{c}],'Color',condicolors(c,:));
colorbar;colormap('jet');clim([cmin cmax]);
subplot(4,5,(c-1)*5+3);
topoplot(mean(H_alpha_LLRR_chan(Inds4_LR(:,c),:)),channels,'nosedir','+X','style','map');
title(['H-alpha ' states4names{c}],'Color',condicolors(c,:));
colorbar;colormap('jet');clim([cmin cmax]);
subplot(4,5,(c-1)*5+4);
topoplot(mean(H_beta_LLRR_chan(Inds4_LR(:,c),:)),channels,'nosedir','+X','style','map');
title(['H-beta ' states4names{c}],'Color',condicolors(c,:));
colorbar;colormap('jet');clim([cmin cmax]);
subplot(4,5,(c-1)*5+5);
topoplot(mean(H_gamma_LLRR_chan(Inds4_LR(:,c),:)),channels,'nosedir','+X','style','map');
title(['H-gamma ' states4names{c}],'Color',condicolors(c,:));
colorbar;colormap('jet');clim([cmin cmax]);
end
sgtitle({['mean H-EEG (+ 500ms) for each condition ^{* PLOT 14-1}'],char(datetime('now'))});
colormap(hnc)
% only H-EEG of beta and gamma show complexity (H>0.5)

%% PLOT 15 corr of H-EEG (- 500ms) and H-intall sessions & 4 states
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
canvas(0.5,0.2);
cmin=-0.7;cmax=0.7;
subplot(1,5,1);topoplot(H_delta_LR_H_LR_corr,channels,'nosedir','+X','style','map');
title('H-delta & H-int');colorbar;colormap('jet');clim([cmin cmax]);
subplot(1,5,2);topoplot(H_theta_LR_H_LR_corr,channels,'nosedir','+X','style','map');
title('H-theta & H-int');colorbar;colormap('jet');clim([cmin cmax]);
subplot(1,5,3);topoplot(H_alpha_LR_H_LR_corr,channels,'nosedir','+X','style','map');
title('H-alpha & H-int');colorbar;colormap('jet');clim([cmin cmax]);
subplot(1,5,4);topoplot(H_beta_LR_H_LR_corr,channels,'nosedir','+X','style','map');
title('H-beta & H-int');colorbar;colormap('jet');clim([cmin cmax]);
subplot(1,5,5);topoplot(H_gamma_LR_H_LR_corr,channels,'nosedir','+X','style','map');
title('H-gamma & H-int');colorbar;colormap('jet');clim([cmin cmax]);
sgtitle({['corr of H-EEG (- 500ms) and H-interval ^{* PLOT 15}'],char(datetime('now'))});
colormap(hnc)

% Compute correlation between H-pow and H-interval in 4 states
H_delta_LR_H_LR_corr4=zeros(4,32);
H_theta_LR_H_LR_corr4=zeros(4,32);
H_alpha_LR_H_LR_corr4=zeros(4,32);
H_beta_LR_H_LR_corr4=zeros(4,32);
H_gamma_LR_H_LR_corr4=zeros(4,32);
for s=1:4
    for c=1:32
        H_delta_LR_H_LR_corr4(s,c)=corr(H_delta_LR_chan(Inds4_LR(:,s),c),H_all_LR(Inds4_LR(:,s)));
        H_theta_LR_H_LR_corr4(s,c)=corr(H_theta_LR_chan(Inds4_LR(:,s),c),H_all_LR(Inds4_LR(:,s)));
        H_alpha_LR_H_LR_corr4(s,c)=corr(H_alpha_LR_chan(Inds4_LR(:,s),c),H_all_LR(Inds4_LR(:,s)));
        H_beta_LR_H_LR_corr4(s,c)=corr(H_beta_LR_chan(Inds4_LR(:,s),c),H_all_LR(Inds4_LR(:,s)));
        H_gamma_LR_H_LR_corr4(s,c)=corr(H_gamma_LR_chan(Inds4_LR(:,s),c),H_all_LR(Inds4_LR(:,s)));
    end
end
% Combine L and R in 4 states (4x5)
canvas(0.5,0.8);
cmin=-0.7;cmax=0.7;
for s=1:4
subplot(4,5,(s-1)*5+1);
topoplot(H_delta_LR_H_LR_corr4(s,:),channels,'nosedir','+X','style','map');
title(['H-delta ' states4names{s} ' & H-int'],'Color',condicolors(s,:));
colorbar;colormap('jet');clim([cmin cmax]);
subplot(4,5,(s-1)*5+2);
topoplot(H_theta_LR_H_LR_corr4(s,:),channels,'nosedir','+X','style','map');
title(['H-theta ' states4names{s} ' & H-int'],'Color',condicolors(s,:));
colorbar;colormap('jet');clim([cmin cmax]);
subplot(4,5,(s-1)*5+3);
topoplot(H_alpha_LR_H_LR_corr4(s,:),channels,'nosedir','+X','style','map');
title(['H-alpha ' states4names{s} ' & H-int'],'Color',condicolors(s,:));
colorbar;colormap('jet');clim([cmin cmax]);
subplot(4,5,(s-1)*5+4);
topoplot(H_beta_LR_H_LR_corr4(s,:),channels,'nosedir','+X','style','map');
title(['H-beta ' states4names{s} ' & H-int'],'Color',condicolors(s,:));
colorbar;colormap('jet');clim([cmin cmax]);
subplot(4,5,(s-1)*5+5);
topoplot(H_gamma_LR_H_LR_corr4(s,:),channels,'nosedir','+X','style','map');
title(['H-gamma ' states4names{s} ' & H-int'],'Color',condicolors(s,:));
colorbar;colormap('jet');clim([cmin cmax]);
end
sgtitle({['corr of H-EEG (- 500ms) and H-interval in 4 states ^{* PLOT 15}'],char(datetime('now'))});
colormap(hnc)
%% PLOT 15-1 corr of H-EEG (+ 500ms) and H-intall sessions & 4 states
% Do the correlation between H-pow and H-interval
% Combine L and R for correlation;
for c=1:32
    H_delta_LLRR_H_LR_corr(c)=corr(H_delta_LLRR_chan(:,c),H_all_LR);
    H_theta_LLRR_H_LR_corr(c)=corr(H_theta_LLRR_chan(:,c),H_all_LR);
    H_alpha_LLRR_H_LR_corr(c)=corr(H_alpha_LLRR_chan(:,c),H_all_LR);
    H_beta_LLRR_H_LR_corr(c)=corr(H_beta_LLRR_chan(:,c),H_all_LR);
    H_gamma_LLRR_H_LR_corr(c)=corr(H_gamma_LLRR_chan(:,c),H_all_LR);
end

% Combine L and R in topoplots (1x5)
canvas(0.5,0.2);
cmin=-0.7;cmax=0.7;
subplot(1,5,1);topoplot(H_delta_LLRR_H_LR_corr,channels,'nosedir','+X','style','map');
title('H-delta & H-int');colorbar;colormap('jet');clim([cmin cmax]);
subplot(1,5,2);topoplot(H_theta_LLRR_H_LR_corr,channels,'nosedir','+X','style','map');
title('H-theta & H-int');colorbar;colormap('jet');clim([cmin cmax]);
subplot(1,5,3);topoplot(H_alpha_LLRR_H_LR_corr,channels,'nosedir','+X','style','map');
title('H-alpha & H-int');colorbar;colormap('jet');clim([cmin cmax]);
subplot(1,5,4);topoplot(H_beta_LLRR_H_LR_corr,channels,'nosedir','+X','style','map');
title('H-beta & H-int');colorbar;colormap('jet');clim([cmin cmax]);
subplot(1,5,5);topoplot(H_gamma_LLRR_H_LR_corr,channels,'nosedir','+X','style','map');
title('H-gamma & H-int');colorbar;colormap('jet');clim([cmin cmax]);
sgtitle({['corr of H-EEG (+ 500ms) and H-interval  ^{* PLOT 15-1}'],char(datetime('now'))});
colormap(hnc)

% Compute correlation between H-pow and H-interval in 4 states
H_delta_LLRR_H_LR_corr4=zeros(4,32);
H_theta_LLRR_H_LR_corr4=zeros(4,32);
H_alpha_LLRR_H_LR_corr4=zeros(4,32);
H_beta_LLRR_H_LR_corr4=zeros(4,32);
H_gamma_LLRR_H_LR_corr4=zeros(4,32);
for s=1:4
    for c=1:32
        H_delta_LLRR_H_LR_corr4(s,c)=corr(H_delta_LLRR_chan(Inds4_LR(:,s),c),H_all_LR(Inds4_LR(:,s)));
        H_theta_LLRR_H_LR_corr4(s,c)=corr(H_theta_LLRR_chan(Inds4_LR(:,s),c),H_all_LR(Inds4_LR(:,s)));
        H_alpha_LLRR_H_LR_corr4(s,c)=corr(H_alpha_LLRR_chan(Inds4_LR(:,s),c),H_all_LR(Inds4_LR(:,s)));
        H_beta_LLRR_H_LR_corr4(s,c)=corr(H_beta_LLRR_chan(Inds4_LR(:,s),c),H_all_LR(Inds4_LR(:,s)));
        H_gamma_LLRR_H_LR_corr4(s,c)=corr(H_gamma_LLRR_chan(Inds4_LR(:,s),c),H_all_LR(Inds4_LR(:,s)));
    end
end
% Combine L and R in 4 states (4x5)
canvas(0.5,0.8);
cmin=-0.7;cmax=0.7;
for s=1:4
subplot(4,5,(s-1)*5+1);
topoplot(H_delta_LLRR_H_LR_corr4(s,:),channels,'nosedir','+X','style','map');
title(['H-delta ' states4names{s} ' & H-int'],'Color',condicolors(s,:));
colorbar;colormap('jet');clim([cmin cmax]);
subplot(4,5,(s-1)*5+2);
topoplot(H_theta_LLRR_H_LR_corr4(s,:),channels,'nosedir','+X','style','map');
title(['H-theta ' states4names{s} ' & H-int'],'Color',condicolors(s,:));
colorbar;colormap('jet');clim([cmin cmax]);
subplot(4,5,(s-1)*5+3);
topoplot(H_alpha_LLRR_H_LR_corr4(s,:),channels,'nosedir','+X','style','map');
title(['H-alpha ' states4names{s} ' & H-int'],'Color',condicolors(s,:));
colorbar;colormap('jet');clim([cmin cmax]);
subplot(4,5,(s-1)*5+4);
topoplot(H_beta_LLRR_H_LR_corr4(s,:),channels,'nosedir','+X','style','map');
title(['H-beta ' states4names{s} ' & H-int'],'Color',condicolors(s,:));
colorbar;colormap('jet');clim([cmin cmax]);
subplot(4,5,(s-1)*5+5);
topoplot(H_gamma_LLRR_H_LR_corr4(s,:),channels,'nosedir','+X','style','map');
title(['H-gamma ' states4names{s} ' & H-int'],'Color',condicolors(s,:));
colorbar;colormap('jet');clim([cmin cmax]);
end
sgtitle({['corr of H-EEG (+ 500ms) and H-interval in 4 states  ^{* PLOT 15-1}'],char(datetime('now'))});
colormap(hnc)

%% PLOT 16 PLS (sum-EEG -> H-int)
% % refer to previous code in sync.m
% open sync.m
% open syn_tryout.m
% % powforpls (28 matches x 3200 predictors) 50freqsx64chan=3200
% % cor(j).BP' (28 corr ceof x 1) 
% [R2,reg,ypred] = npls_pred(powforpls,cor(j).BP',1);
% weights = reshape(reg{1},50,64); % reg{1} (3200 x 1)
% R2 = R2;
% ypred = ypred;

% sum-EEG -> H-int
% Power
pow5forpls=[];pow5forpls2=[];pow5forpls3=[];
% % fix the scale in the data (% already fixed in PLOT 13)
% delta_LR_chan = delta_LR_chan./(ones(288,1)*std(delta_LR_chan)); %organized in PLOT 13
% theta_LR_chan = theta_LR_chan./(ones(288,1)*std(theta_LR_chan));
% alpha_LR_chan = alpha_LR_chan./(ones(288,1)*std(alpha_LR_chan));
% beta_LR_chan = beta_LR_chan./(ones(288,1)*std(beta_LR_chan));
% gamma_LR_chan = gamma_LR_chan./(ones(288,1)*std(gamma_LR_chan));
% arrage the band powers in "trials x chan x freq" (288 x 32 x 5)
pow5forpls=cat(3,delta_LR_chan,theta_LR_chan,alpha_LR_chan,beta_LR_chan,gamma_LR_chan);
% then switch dismenstion to "freq x chan x trials" (5 x 32 x 288)
pow5forpls2=permute(pow5forpls,[3 2 1]);
% then stack them in 2-d " trials x (freq x chan)" 
pow5forpls3=reshape(pow5forpls2,5*32,288)'; % 288x160 elements extracted colums-wise
% (colums: delta chan1, theta chan 1, alpha chan1, beta chan1, gamma chan1, delta chan2...)
% clear workspace
clear pow5forpls pow5forpls2

% H-int
H_all_LR; % form PLOT 13

addpath /home/zhibinz2/Documents/GitHub/matlab/ramesh/plsmodel
addpath(genpath('/home/zhibinz2/Documents/GitHub/matlab/external/')); 
addpath /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/PLS

% (ALL states: 5freq x 32chan = 160 predictors x 288 trials)
% -updated to use mynpls_pred function
reg=[];ypred_fit=[];X_MC=[];Y_MC=[];
clear plsmodel;
[reg,ypred_fit,X_MC,Y_MC] = mynpls_pred(pow5forpls3,H_all_LR);
ssOr=sum(sum((Y_MC-mean(Y_MC)).^2))
AIC=nan(1,10);R2s=nan(1,10);ssERs=nan(1,10);
for Fac=1:10;
    R2=[];ssEr=[];
    ssEr=sum(sum((Y_MC-ypred_fit(:,Fac)).^2))
    ssERs(Fac)=ssEr;
    R2=100*(1-ssEr/ssOr)
    R2s(Fac)=R2;
    % Cal AIC
    AIC(Fac)=log(ssEr)+2*(Fac);
end
canvas(0.3,0.25);
subplot(2,1,1);
plot(1:10,AIC,'.','MarkerSize',20);xlabel('nFac');ylabel('AIC');
xlim([0 11]);
subtitle('All states');
subplot(2,1,2);
plot(1:10,R2s,'r.','MarkerSize',20);xlabel('nFac');ylabel('R^2');
xlim([0 11]);
subtitle('All states');
sgtitle({['AIC & R^2 for PLS model in all statues: sum-EEG(-500ms) -> H-int^{ *PLOT16}'],char(datetime('now'))});
% Try cross validation (myxvalidation)
% select the Fac with minimal AIC
Fac=1
% do the #-fold cross validation
ssEr=[];R2=[];
X=pow5forpls3;Y=H_all_LR;
nSet=1;
[ssEr,R2]  = myxvalidation(X,Y,Fac,nSet);
% reshape regs{1:4}{Fac} 
clear plsmodel
% reshape from 160 x 1 back to 5 x 32 (freq x chan)
plsmodel.weights=reshape(reg{Fac},5,32); 
canvas(0.2,0.2)
cmin=-8e-3;cmax=8e-3;
imagesc(plsmodel.weights);colorbar; % caxis([-2 2]*10E-7); % by default, imagesc reverse the Y 
yticks([1:5]);yticklabels({'Delta','Theta','Alpha','Beta','Gamma'});
colormap('jet'); clim([cmin cmax]);
set(gca, 'YDir','normal');
xticks([1:32]);xticklabels([labels]);xtickangle(90);grid on;
title({['PLS model (R2= ' num2str(round(R2,1)) ') in all statues: sum-EEG(-500ms) -> H-int  ^{* PLOT 16}'],char(datetime('now'))});
set(gcf,'color','w'); % set background white for copying in ubuntu


% (Each of the 4 states: 5freq x 32 chan = 160 predictors x 48
% trials)-updated to use mynpls_pred function
cd /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/PLS
regs=[];AIC=nan(2,10);R2s=nan(4,10);ssERs=nan(4,10);
for c=1:4 % four states
    R2=[];reg=[];ypred=[];
    [reg,ypred_fit,X_MC,Y_MC] = mynpls_pred(pow5forpls3(Inds4_LR(:,c),:),H_all_LR(Inds4_LR(:,c)));
    regs{c}=reg;
    ssOr=sum(sum((Y_MC-mean(Y_MC)).^2))
    for Fac=1:10;
        R2=[];ssEr=[];
        ssEr=sum(sum((Y_MC-ypred_fit(:,Fac)).^2))
        ssERs(c,Fac)=ssEr;
        R2=100*(1-ssEr/ssOr)
        R2s(c,Fac)=R2;
        % Cal AIC
        AIC(c,Fac)=log(ssEr)+2*(Fac);
    end
end
figure;
for c=1:4
    subplot(2,4,c);
    plot(1:10,AIC(c,:),'.','MarkerSize',20);xlabel('nFac');ylabel('AIC');
    ylim([0 17]);xlim([0 11]);
    subtitle(states4names{c},'Color',condicolors(c,:));
    subplot(2,4,4+c);
    plot(1:10,R2s(c,:),'r.','MarkerSize',20);xlabel('nFac');ylabel('R^2');
    ylim([10 110]);xlim([0 11]);
    subtitle(states4names{c},'Color',condicolors(c,:));
end
sgtitle({['AIC & R^2 for PLS model in all 4 statues: sum-EEG(-500ms) -> H-int^{ *PLOT16}'],char(datetime('now'))});

% Try cross validation (myxvalidation)
% select the Fac with minimal AIC
Fac=1
ssErs=[];R2s=[];
% do the cross validation
for c=1:4
    ssEr=[];R2=[];
    X=pow5forpls3(Inds4_LR(:,c),:);Y=H_all_LR(Inds4_LR(:,c));
    [ssEr,R2]  = myxvalidation(X,Y,Fac,nSet);
    ssErs(c)=ssEr;
    R2s(c)=R2;
end
% reshape regs{1:4}{Fac} 
clear plsmodel
for c=1:4
    plsmodel(c).weights=reshape(regs{c}{Fac},5,32); 
end
% Plot the 4 states
canvas(0.4,0.4);
cmin=-0.01;cmax=0.01;
for c=1:4
    subplot(2,2,c);
    imagesc(plsmodel(c).weights);colorbar; 
    yticks([1:5]);yticklabels({'Delta','Theta','Alpha','Beta','Gamma'});
    set(gca, 'YDir','normal');
    xticks([1:32]);xticklabels([labels]);xtickangle(90);
    colormap('jet'); clim([cmin cmax]);
    % caxis([-0.1 0.1]);
    title([states4names{c} ': PLS model (R^2= ' num2str(round(R2s(c),1)) '  Fac= ' num2str(Fac) ') '], ...
        'Color',condicolors(c,:));
    grid on;
end
sgtitle({['PLS model: sum-EEG(-500ms) -> H-int ^{* PLOT 16}'],char(datetime('now'))});
set(gcf,'color','w'); % set background white for copying in ubuntu

% topoplot for uncouple and mutual (now 4 states)
addpath /home/zhibinz2/Documents/GitHub/matlab-archive/hnlcode/common/gen_code/color
hnc = hotncold(100);
band5names={'Delta','Theta','Alpha','Beta','Gamma'};
states2names={'Uncouple','Mutual'};
states4names; 
% figure;
canvas(0.3,0.8)
for states4=1;
c=[1 2 3 4];
for s=1:4
    for b=1:5
        subplot(4,5,(s-1)*5+b)
        topoplot(plsmodel(c(s)).weights(b,:),channels,'nosedir','+X','style','map');
        % colorbar; 
        % colormap(hnc); % colormap('jet');
        clim([cmin cmax]);
%         if b==3;title({['PLS model: sum-EEG(-500ms) -> H-int ^{* PLOT 16}'], ...
%                 [states4names{c(s)} '(R2= ' num2str(round(R2s(c(s)),1)) ...
%             '  Fac= ' num2str(Fac) ') ']},'Color',condicolors(c(s),:));end
        % subtitle(band5names(b));
    end
end
colormap(hnc)
% https://www.mathworks.com/help/matlab/ref/matlab.graphics.illustration.colorbar-properties.html
cb=colorbar;
cb.AxisLocation = 'out';
cb.Position = [0.92 0.15 0.01 0.75];
set(gcf,'color','w'); % set background white for copying in ubuntu
delete(findall(gcf,'type','annotation'))
h0=annotation('textbox',[0.17 0.95 0.05 0.03],'string','Delta','color',[0 0 0])
h1=annotation('textbox',[0.33 0.95 0.05 0.03],'string','Theta','color',[0 0 0])
h2=annotation('textbox',[0.5 0.95 0.05 0.03],'string','Alpha','color',[0 0 0])
h3=annotation('textbox',[0.66 0.95 0.05 0.03],'string','Beta','color',[0 0 0])
h4=annotation('textbox',[0.82 0.95 0.05 0.03],'string','Gamma','color',[0 0 0])
v0=annotation('textbox',[0.1 0.15 0.05 0.03],'string',...
    {[states4names{c(4)}], [' (R^2= ' num2str(round(R2s(c(4)),1)) ...
            '  Fac= ' num2str(Fac) ') ']}...
            ,'color',condicolors(4,:));
set(v0,'Rotation',90);
v1=annotation('textbox',[0.1 0.35 0.05 0.03],'string',...
    {[states4names{c(3)}], [' (R^2= ' num2str(round(R2s(c(3)),1)) ...
            '  Fac= ' num2str(Fac) ') ']} ...
            ,'color',condicolors(3,:));
set(v1,'Rotation',90);
v2=annotation('textbox',[0.1 0.56 0.05 0.03],'string',...
    {[states4names{c(2)}], [' (R^2= ' num2str(round(R2s(c(2)),1)) ...
            '  Fac= ' num2str(Fac) ') ']} ...
            ,'color',condicolors(2,:));
set(v2,'Rotation',90);
v3=annotation('textbox',[0.1 0.79 0.05 0.03],'string',...
    {[states4names{c(1)}], [' (R^2= ' num2str(round(R2s(c(1)),1)) ...
            '  Fac= ' num2str(Fac) ') ']} ...
            ,'color',condicolors(1,:));
set(v3,'Rotation',90);
end
sg=annotation('textbox',[0.3 0.01 0.4 0.05],'string',...
    {['PLS model: sum-EEG(-500ms) -> H-int ^{* PLOT 16}'],char(datetime('now'))});
set([h0 h1 h2 h3 h4 v0 v1 v2 v3], 'fitboxtotext','on',...
    'edgecolor','none')
% topoplot only for mutual ( 1 x 5 freq)
canvas(0.3,0.2)
for states1=1;
cmin=-0.01;cmax=0.01;
c=[4];
for s=1    
    for b=1:5
        subplot(1,5,(s-1)*5+b)
        topoplot(plsmodel(c(s)).weights(b,:),channels,'nosedir','+X','style','map');
        clim([cmin cmax]);
    end
end
colormap(hnc)
cb=colorbar;
cb.AxisLocation = 'out';
cb.Position = [0.92 0.15 0.01 0.75];
set(gcf,'color','w'); % set background white for copying in ubuntu
if true %annotations
delete(findall(gcf,'type','annotation'))
h0=annotation('textbox',[0.17 0.95 0.05 0.03],'string','Delta','color',[0 0 0])
h1=annotation('textbox',[0.33 0.95 0.05 0.03],'string','Theta','color',[0 0 0])
h2=annotation('textbox',[0.5 0.95 0.05 0.03],'string','Alpha','color',[0 0 0])
h3=annotation('textbox',[0.66 0.95 0.05 0.03],'string','Beta','color',[0 0 0])
h4=annotation('textbox',[0.82 0.95 0.05 0.03],'string','Gamma','color',[0 0 0])
v0=annotation('textbox',[0.1 0.4 0.05 0.03],'string',...
    {[states4names{c(1)}], [' (R^2= ' num2str(round(R2s(c(1)),1)) ...
            '  Fac= ' num2str(Fac) ') ']}...
            ,'color',condicolors(4,:));
set(v0,'Rotation',90);
set([h0 h1 h2 h3 h4 v0 v1 v2 v3], 'fitboxtotext','on',...
    'edgecolor','none')
end
end
sg=annotation('textbox',[0.3 0.01 0.3 0.15],'string',...
    {['PLS model: sum-EEG(-500ms) -> H-int ^{* PLOT 16}'],char(datetime('now'))});



% 3 states - combine uncouple state and leading state and call it "independent"
% cannot combine, subject can distinguish between uncouple and leading
states3names={'Independent','Following','Mutual'}; 
Inds3={[Inds4_LR(:,1);Inds4_LR(:,2)],Inds4_LR(:,3),Inds4_LR(:,4)};
canvas(0.6,0.2);
cmin=-8e-3;cmax=8e-3;
for c=1:3 % four states
    R2=[];reg=[];ypred=[];
    clear plsmodel;
    [R2,reg,ypred] = npls_pred(pow5forpls3(Inds3{c},:),H_all_LR(Inds3{c}),1);
    plsmodel(c).weights=reshape(reg{1},5,32);
    plsmodel(c).R2 = R2;
    plsmodel(c).ypred = ypred;
    subplot(1,3,c);
    imagesc(plsmodel(c).weights);colorbar; 
    yticks([1:5]);yticklabels({'Delta','Theta','Alpha','Beta','Gamma'});
    set(gca, 'YDir','normal');
    xticks([1:32]);xticklabels([labels]);xtickangle(90);
    colormap('jet'); clim([cmin cmax]);% caxis([-0.1 0.1]);
    title([states3names{c} ': PLS model (R2= ' num2str(round(R2,1)) ')'],'Color',condicolors(c,:));
    grid on;
end
sgtitle({['PLS model: sum-EEG(-500ms) -> H-int'],char(datetime('now'))});

%% PLOT 16-0 PLS (sum-EEG -> H-int-diff)
pow5forpls3;
H_all_LR;
% H diff
H_all_diff_LR=[(H_all_LR(1:144)-H_all_LR(145:288)); (H_all_LR(145:288)-H_all_LR(1:144))];


% (ALL states: 5freq x 32chan = 160 predictors x 288 trials)
% -updated to use mynpls_pred function
reg=[];ypred_fit=[];X_MC=[];Y_MC=[];
clear plsmodel;
[reg,ypred_fit,X_MC,Y_MC] = mynpls_pred(pow5forpls3,H_all_diff_LR);
ssOr=sum(sum((Y_MC-mean(Y_MC)).^2))
AIC=nan(1,10);R2s=nan(1,10);ssERs=nan(1,10);
for Fac=1:10;
    R2=[];ssEr=[];
    ssEr=sum(sum((Y_MC-ypred_fit(:,Fac)).^2))
    ssERs(Fac)=ssEr;
    R2=100*(1-ssEr/ssOr)
    R2s(Fac)=R2;
    % Cal AIC
    AIC(Fac)=log(ssEr)+2*(Fac);
end
canvas(0.3,0.25);
subplot(2,1,1);
plot(1:10,AIC,'.','MarkerSize',20);xlabel('nFac');ylabel('AIC');
xlim([0 11]);
subtitle('All states');
subplot(2,1,2);
plot(1:10,R2s,'r.','MarkerSize',20);xlabel('nFac');ylabel('R^2');
xlim([0 11]);
subtitle('All states');
sgtitle({['AIC & R^2 for PLS model in all statues: sum-EEG(-500ms) -> H-diff^{ *PLOT 16-0}'],char(datetime('now'))});
% Try cross validation (myxvalidation)
% select the Fac with minimal AIC
Fac=1
% do the #-fold cross validation and plot the results
if true % cross validation
ssEr=[];R2=[];
X=pow5forpls3;Y=H_all_diff_LR;
nSet=1;
[ssEr,R2]  = myxvalidation(X,Y,Fac,nSet);
% reshape regs{1:4}{Fac} 
clear plsmodel
% reshape from 160 x 1 back to 5 x 32 (freq x chan)
plsmodel.weights=reshape(reg{Fac},5,32); 
canvas(0.2,0.2)
cmin=-8e-3;cmax=8e-3;
imagesc(plsmodel.weights);colorbar; % caxis([-2 2]*10E-7); % by default, imagesc reverse the Y 
yticks([1:5]);yticklabels({'Delta','Theta','Alpha','Beta','Gamma'});
colormap('jet'); clim([cmin cmax]);
set(gca, 'YDir','normal');
xticks([1:32]);xticklabels([labels]);xtickangle(90);grid on;
title({['PLS model (R2= ' num2str(round(R2,1)) ') in all statues: sum-EEG(-500ms) -> H-diff  ^{* PLOT 16}'],char(datetime('now'))});
set(gcf,'color','w'); % set background white for copying in ubuntu
end

% PLS on 4 states 
cd /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/PLS
regs=[];AIC=nan(2,10);R2s=nan(4,10);ssERs=nan(4,10);
for c=1:4 % four states
    reg=[];ypred=[];ypred_fit=[];
    [reg,ypred_fit,X_MC,Y_MC] = mynpls_pred(pow5forpls3(Inds4_LR(:,c),:),H_all_diff_LR(Inds4_LR(:,c)));
    regs{c}=reg;
    ssOr=sum(sum((Y_MC-mean(Y_MC)).^2))
    for Fac=1:10;
        R2=[];ssEr=[];
        ssEr=sum(sum((Y_MC-ypred_fit(:,Fac)).^2))
        ssERs(c,Fac)=ssEr;
        R2=100*(1-ssEr/ssOr)
        R2s(c,Fac)=R2;
        % Cal AIC
        AIC(c,Fac)=log(ssEr)+2*(Fac);
    end
end
figure;
for c=1:4
    subplot(2,4,c);
    plot(1:10,AIC(c,:),'.','MarkerSize',20);xlabel('nFac');ylabel('AIC');
    ylim([0 17]);xlim([0 11]);
    subtitle(states4names{c},'Color',condicolors(c,:));
    subplot(2,4,4+c);
    plot(1:10,R2s(c,:),'r.','MarkerSize',20);xlabel('nFac');ylabel('R^2');
    ylim([10 110]);xlim([0 11]);
    subtitle(states4names{c},'Color',condicolors(c,:));
end
sgtitle({['AIC & R^2 for PLS model in all 4 statues: H-EEG(-500ms) -> H-diff ^{ *PLOT16-0}'],char(datetime('now'))});
% Try cross validation (myxvalidation)
% select the Fac with minimal AIC
Fac=1;
ssErs=[];R2s=[];
% do the cross validation
for c=1:4
    ssEr=[];R2=[];
    X=pow5forpls3(Inds4_LR(:,c),:);Y=H_all_diff_LR(Inds4_LR(:,c));
    [ssEr,R2]  = myxvalidation(X,Y,Fac,nSet);
    ssErs(c)=ssEr;
    R2s(c)=R2;
end
% reshape regs{1:4}{Fac} 
clear plsmodel
for c=1:4
    plsmodel(c).weights=reshape(regs{c}{Fac},5,32); 
end
% Plot the 4 states
canvas(0.4,0.4);
cmin=-0.1;cmax=0.1;
for c=1:4
    subplot(2,2,c);
    imagesc(plsmodel(c).weights);colorbar; 
    yticks([1:5]);yticklabels({'Delta','Theta','Alpha','Beta','Gamma'});
    set(gca, 'YDir','normal');
    xticks([1:32]);xticklabels([labels]);xtickangle(90);
    colormap('jet'); clim([cmin cmax]);
    % caxis([-0.1 0.1]);
    title([states4names{c} ': PLS model (R^2= ' num2str(round(R2s(c),1)) '  Fac= ' num2str(Fac) ') '], ...
        'Color',condicolors(c,:));
    grid on;
end
sgtitle({['PLS model: H-EEG(-500ms) -> H-diff ^{* PLOT 16-0}'],char(datetime('now'))});
set(gcf,'color','w'); % set background white for copying in ubuntu


% PLS on 4 states separating synch and synco
% Inds_synch4_LR; Inds_synco4_LR; % from SECT 11-2
syn4_LR=cat(3,Inds_synch4_LR,Inds_synco4_LR)
regs=[];AIC=nan(2,4,10);R2s=nan(2,4,10);ssERs=nan(2,4,10);
for t=1:2 % two syn types
    select_Inds4_LR=syn4_LR(:,:,t);
for c=1:4 % four states
    reg=[];ypred=[];ypred_fit=[];ssOr=[];
    [reg,ypred_fit,X_MC,Y_MC] = mynpls_pred(pow5forpls3(select_Inds4_LR(:,c),:),H_all_diff_LR(select_Inds4_LR(:,c)));
    regs{t,c}=reg;
    ssOr=sum(sum((Y_MC-mean(Y_MC)).^2))
    for Fac=1:10;
        R2=[];ssEr=[];
        ssEr=sum(sum((Y_MC-ypred_fit(:,Fac)).^2))
        ssERs(t,c,Fac)=ssEr;
        R2=100*(1-ssEr/ssOr)
        R2s(t,c,Fac)=R2;
        % Cal AIC
        AIC(t,c,Fac)=log(ssEr)+2*(Fac);
    end
end
end
for t=1:2
    figure;
for c=1:4
    subplot(2,4,c);
    plot(1:10,squeeze(AIC(t,c,:)),'.','MarkerSize',20);xlabel('nFac');ylabel('AIC');
    ylim([0 17]);xlim([0 11]);
    subtitle(states4names{c},'Color',condicolors(c,:));
    subplot(2,4,4+c);
    plot(1:10,squeeze(R2s(t,c,:)),'r.','MarkerSize',20);xlabel('nFac');ylabel('R^2');
    ylim([10 110]);xlim([0 11]);
    subtitle(states4names{c},'Color',condicolors(c,:));
end
if t==1;
   sgtitle({['AIC & R^2 for PLS model in all 4 statues for synch sessions: H-EEG(-500ms) -> H-diff ^{ *PLOT16-0}'],...
       char(datetime('now'))},'Color',syn2colors(t,:));
else t==2;
   sgtitle({['AIC & R^2 for PLS model in all 4 statues for synco sessions: H-EEG(-500ms) -> H-diff ^{ *PLOT16-0}'],...
       char(datetime('now'))},'Color',syn2colors(t,:));
end
end
% Try cross validation (myxvalidation)
% select the Fac with minimal AIC
Fac=1;
ssErs=[];R2s=[];
% do the cross validation
for t=1:2
    select_Inds4_LR=syn4_LR(:,:,t);
for c=1:4
    ssEr=[];R2=[];
    X=pow5forpls3(select_Inds4_LR(:,c),:);Y=H_all_diff_LR(select_Inds4_LR(:,c));
    [ssEr,R2]  = myxvalidation(X,Y,Fac,nSet);
    ssErs(t,c)=ssEr;
    R2s(t,c)=R2;
end
end
% reshape regs{1:4}{Fac} 
clear plsmodel
for t=1:2
    % select_Inds4_LR=syn4_LR(:,:,t);
for c=1:4
    plsmodel(t,c).weights=reshape(regs{t,c}{Fac},5,32); 
end
end
% Plot the 4 states
cmin=-0.02;cmax=0.02;
for t=1:2
    % select_Inds4_LR=syn4_LR(:,:,t);
    canvas(0.4,0.4);
for c=1:4
    subplot(2,2,c);
    imagesc(plsmodel(t,c).weights);colorbar; 
    yticks([1:5]);yticklabels({'Delta','Theta','Alpha','Beta','Gamma'});
    set(gca, 'YDir','normal');
    xticks([1:32]);xticklabels([labels]);xtickangle(90);
    colormap('jet'); clim([cmin cmax]);
    % caxis([-0.1 0.1]);
    title([states4names{c} ': PLS model (R^2= ' num2str(round(R2s(t,c),1)) '  Fac= ' num2str(Fac) ') '], ...
        'Color',condicolors(c,:));
    grid on;
end
if t==1;
sgtitle({['PLS model for synch: H-EEG(-500ms) -> H-diff ^{* PLOT 16-0}'],...
    char(datetime('now'))},'Color',syn2colors(t,:));
else t==2;
sgtitle({['PLS model for synco: H-EEG(-500ms) -> H-diff ^{* PLOT 16-0}'],...
    char(datetime('now'))},'Color',syn2colors(t,:));
end
set(gcf,'color','w'); % set background white for copying in ubuntu
end
%% PLOT 16-1 PLS (H-EEG -> H-int)
% H-EEG -> H-int
H5forpls=[];H5forpls2=[];H5forpls3=[];
% arrage the band powers in "trials x chan x freq" (288 x 32 x 5)
H5forpls=cat(3,H_delta_LR_chan,H_theta_LR_chan,H_alpha_LR_chan,H_beta_LR_chan,H_gamma_LR_chan);
% then switch dismenstion to "freq x chan x trials" (5 x 32 x 288)
H5forpls2=permute(H5forpls,[3 2 1]);
% then stack them in 2-d " trials x (freq x chan)" 
H5forpls3=reshape(H5forpls2,5*32,288)'; % 288x160 elements extracted colums-wise
% (colums: delta chan1, theta chan 1, alpha chan1, beta chan1, gamma chan1, delta chan2...)

% H-int
H_all_LR; % form PLOT 13

addpath /home/zhibinz2/Documents/GitHub/matlab/ramesh/plsmodel
addpath(genpath('/home/zhibinz2/Documents/GitHub/matlab/external/')); 

% (ALL states: 5freq x 32chan = 160 predictors x 192 trials)
% -updated to use mynpls_pred function
reg=[];ypred_fit=[];X_MC=[];Y_MC=[];
clear plsmodel;
[reg,ypred_fit,X_MC,Y_MC] = mynpls_pred(H5forpls3,H_all_LR);
ssOr=sum(sum((Y_MC-mean(Y_MC)).^2))
AIC=nan(1,10);R2s=nan(1,10);ssERs=nan(1,10);
for Fac=1:10;
    R2=[];ssEr=[];
    ssEr=sum(sum((Y_MC-ypred_fit(:,Fac)).^2))
    ssERs(Fac)=ssEr;
    R2=100*(1-ssEr/ssOr)
    R2s(Fac)=R2;
    % Cal AIC
    AIC(Fac)=log(ssEr)+2*(Fac);
end
figure;
subplot(2,1,1);
plot(1:10,AIC,'.','MarkerSize',20);xlabel('nFac');ylabel('AIC');
xlim([0 11]);
subtitle('All states');
subplot(2,1,2);
plot(1:10,R2s,'r.','MarkerSize',20);xlabel('nFac');ylabel('R^2');
xlim([0 11]);
subtitle('All states');
sgtitle({['AIC & R^2 for PLS model in all statues: H-EEG(-500ms) -> H-int^{ *PLOT16-1}'],char(datetime('now'))});
% Try cross validation (myxvalidation)
% select the Fac with minimal AIC
Fac=1
% do the 8-fold cross validation
ssEr=[];R2=[];
X=H5forpls3;Y=H_all_LR;
[ssEr,R2]  = myxvalidation(X,Y,Fac,nSet);
% reshape regs{1:4}{Fac} 
clear plsmodel
% reshape from 160 x 1 back to 5 x 32 (freq x chan)
plsmodel.weights=reshape(reg{Fac},5,32); 
canvas(0.2,0.2)
cmin=-0.15;cmax=0.15;
imagesc(plsmodel.weights);colorbar; % caxis([-2 2]*10E-7); % by default, imagesc reverse the Y 
yticks([1:5]);yticklabels({'Delta','Theta','Alpha','Beta','Gamma'});
colormap('jet'); clim([cmin cmax]);
set(gca, 'YDir','normal');
xticks([1:32]);xticklabels([labels]);xtickangle(90);grid on;
title({['PLS model (R2= ' num2str(round(R2,1)) ') in all statues: H-EEG(-500ms) -> H-int  ^{* PLOT 16-1}'],char(datetime('now'))});
set(gcf,'color','w'); % set background white for copying in ubuntu

% (Each of the 4 states: 5freq x 32 chan = 160 predictors x 48 trials)
% -updated to use mynpls_pred function
cd /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/PLS
regs=[];AIC=nan(2,10);R2s=nan(4,10);ssERs=nan(4,10);
for c=1:4 % four states
    reg=[];ypred=[];ypred_fit=[];
    [reg,ypred_fit,X_MC,Y_MC] = mynpls_pred(H5forpls3(Inds4_LR(:,c),:),H_all_LR(Inds4_LR(:,c)));
    regs{c}=reg;
    ssOr=sum(sum((Y_MC-mean(Y_MC)).^2))
    for Fac=1:10;
        R2=[];ssEr=[];
        ssEr=sum(sum((Y_MC-ypred_fit(:,Fac)).^2))
        ssERs(c,Fac)=ssEr;
        R2=100*(1-ssEr/ssOr)
        R2s(c,Fac)=R2;
        % Cal AIC
        AIC(c,Fac)=log(ssEr)+2*(Fac);
    end
end
figure;
for c=1:4
    subplot(2,4,c);
    plot(1:10,AIC(c,:),'.','MarkerSize',20);xlabel('nFac');ylabel('AIC');
    ylim([0 17]);xlim([0 11]);
    subtitle(states4names{c},'Color',condicolors(c,:));
    subplot(2,4,4+c);
    plot(1:10,R2s(c,:),'r.','MarkerSize',20);xlabel('nFac');ylabel('R^2');
    ylim([10 110]);xlim([0 11]);
    subtitle(states4names{c},'Color',condicolors(c,:));
end
sgtitle({['AIC & R^2 for PLS model in all 4 statues: H-EEG(-500ms) -> H-int ^{ *PLOT16-1}'],char(datetime('now'))});
% Try cross validation (myxvalidation)
% select the Fac with minimal AIC
Fac=1;
ssErs=[];R2s=[];
% do the cross validation
for c=1:4
    ssEr=[];R2=[];
    X=H5forpls3(Inds4_LR(:,c),:);Y=H_all_LR(Inds4_LR(:,c));
    [ssEr,R2]  = myxvalidation(X,Y,Fac,nSet);
    ssErs(c)=ssEr;
    R2s(c)=R2;
end
% reshape regs{1:4}{Fac} 
clear plsmodel
for c=1:4
    plsmodel(c).weights=reshape(regs{c}{Fac},5,32); 
end
% Plot the 4 states
canvas(0.4,0.4);
cmin=-0.1;cmax=0.1;
for c=1:4
    subplot(2,2,c);
    imagesc(plsmodel(c).weights);colorbar; 
    yticks([1:5]);yticklabels({'Delta','Theta','Alpha','Beta','Gamma'});
    set(gca, 'YDir','normal');
    xticks([1:32]);xticklabels([labels]);xtickangle(90);
    colormap('jet'); clim([cmin cmax]);
    % caxis([-0.1 0.1]);
    title([states4names{c} ': PLS model (R^2= ' num2str(round(R2s(c),1)) '  Fac= ' num2str(Fac) ') '], ...
        'Color',condicolors(c,:));
    grid on;
end
sgtitle({['PLS model: H-EEG(-500ms) -> H-int ^{* PLOT 16-1}'],char(datetime('now'))});
set(gcf,'color','w'); % set background white for copying in ubuntu
% topoplot for uncouple and mutual
addpath /home/zhibinz2/Documents/GitHub/matlab-archive/hnlcode/common/gen_code/color
hnc = hotncold(100);
band5names={'Delta','Theta','Alpha','Beta','Gamma'};
states2names={'Uncouple','Mutual'};
states4names;
figure; % canvas(0.9,0.5)
c=[1 4];
for s=1:2
    for b=1:5
        subplot(2,5,(s-1)*5+b)
        topoplot(plsmodel(c(s)).weights(b,:),channels,'nosedir','+X','style','map');colorbar; 
        colormap(hnc); % colormap('jet');
        clim([cmin cmax]);
        if b==3;title({['PLS model: H-EEG(-500ms) -> H-int ^{* PLOT 16-1}'], ...
                [states4names{c(s)} '(R2= ' num2str(round(R2s(c(s)),1)) ...
            '  Fac= ' num2str(Fac) ') '],char(datetime('now'))},'Color',condicolors(c(s),:));end
        subtitle(band5names(b));
    end
end
set(gcf,'color','w'); % set background white for copying in ubuntu


% 3 states - combine uncouple state and leading state and call it "independent"
% cannot combine, subject can distinguish between uncouple and leading
states3names={'independent','following','mutual'};
Inds3={[Inds4_LR(:,1);Inds4_LR(:,2)],Inds4_LR(:,3),Inds4_LR(:,4)};
canvas(0.6,0.2);
cmin=-0.15;cmax=0.15;
for c=1:3 % three states
    R2=[];reg=[];ypred=[];
    [R2,reg,ypred] = npls_pred(H5forpls3(Inds3{c},:),H_all_LR(Inds3{c}),1);
    plsmodel(c).weights=reshape(reg{1},5,32);
    plsmodel(c).R2 = R2;
    plsmodel(c).ypred = ypred;
    subplot(1,3,c);
    imagesc(plsmodel(c).weights);colorbar; 
    yticks([1:5]);yticklabels({'Delta','Theta','Alpha','Beta','Gamma'});
    set(gca, 'YDir','normal');
    xticks([1:32]);xticklabels([labels]);xtickangle(90);
    colormap('jet');clim([cmin cmax]);% caxis([-0.1 0.1]);
    title([states3names{c} ': PLS model (R2= ' num2str(round(R2,1)) ')'],'Color',condicolors(c,:));
    grid on;
end
sgtitle({['PLS model: H-EEG(-500ms) -> H-int'],char(datetime('now'))});

%% PLOT 17 Lasso (sum-EEG(-500ms)-> H-int)
% https://www.mathworks.com/help/stats/lassoglm.html
% openExample('stats/RemoveRedundantPredictorsUsingLassoMethodExample')
% open lassoglm.m

% sum-EEG -> H-int
% (ALL states: 5freq x 32chan = 160 predictors x 192 trials)
% B  = lassoglm(pow5forpls3,H_all_LR);
clear B FitInfo
[B,FitInfo]  = lassoglm(pow5forpls3,H_all_LR,'normal','CV',5);% try normal or gamma
% Find the coefficient vector for the 75th Lambda value in B.
% B(:,75)
imagesc(B);colorbar;
% Examine the cross-validation plot to see the effect of the Lambda regularization parameter.
lassoPlot(B,FitInfo,'plottype','CV'); 
legend('show') % Show legend
% locate the point with minimum cross-validation error plus one standard deviation
idxLambda1SE = FitInfo.Index1SE;
min1coefs = find(B(:,idxLambda1SE));
Lambda_select=idxLambda1SE;
if isempty(min1coefs)
   % locate the Lambda with minimum cross-validation error
    idxLambdaMinDeviance = FitInfo.IndexMinDeviance;
    mincoefs = find(B(:,idxLambdaMinDeviance))
    Lambda_select=idxLambdaMinDeviance;
end
% reshape from 160 x 1 back to 5 x 32 (freq x chan)
Lambda_coef = reshape(B(:,Lambda_select),5,32); 
% Lambda_coef = reshape(B(:,max(idxLambda1SE)),5,32);
cmin=-8e-3;cmax=8e-3;
canvas(0.15,0.15)
imagesc(Lambda_coef);colorbar;
yticks([1:5]);yticklabels({'Delta','Theta','Alpha','Beta','Gamma'});
set(gca, 'YDir','normal');
colormap('jet'); clim([cmin cmax]);
xticks([1:32]);xticklabels([labels]);xtickangle(90);
title({['Lassoglm in all 4 statues: sum-EEG(-500ms) -> H-int'],char(datetime('now'))});
subtitle(['coefficients for the ' num2str(Lambda_select) 'th Lambda value' ...
    ' (Deviance=' num2str(round(FitInfo.Deviance(Lambda_select),1)) ')']);
grid on;

% (Each of the 4 states: 5freq x 32 chan = 160 predictors x 48 trials)
canvas(0.4,0.4);
cmin=-8e-2;cmax=8e-2;
for c=1:4 % three states
    clear B FitInfo idxLambdaMinDeviance idxLambda1SE Lambda_coef
    [B,FitInfo]  = lassoglm(pow5forpls3(Inds4_LR(:,c),:),H_all_LR(Inds4_LR(:,c)),'normal','CV',length(H_all_LR(Inds4_LR(:,c))));% normal or gamma
    % locate the point with minimum cross-validation error plus one standard deviation
    idxLambda1SE = FitInfo.Index1SE;
    min1coefs = find(B(:,idxLambda1SE));
    Lambda_select=idxLambda1SE;
    if isempty(min1coefs)
       % locate the Lambda with minimum cross-validation error
        idxLambdaMinDeviance = FitInfo.IndexMinDeviance;
        mincoefs = find(B(:,idxLambdaMinDeviance));
        Lambda_select=idxLambdaMinDeviance;
    end
    % reshape from 160 x 1 back to 5 x 32 (freq x chan)
    % Lambda_coef = reshape(B(:,max(idxLambdaMinDeviance)),5,32); 
    Lambda_coef = reshape(B(:,Lambda_select),5,32);
    subplot(2,2,c);
    imagesc(Lambda_coef);colorbar;
    yticks([1:5]);yticklabels({'Delta','Theta','Alpha','Beta','Gamma'});
    set(gca, 'YDir','normal');
    colormap('jet'); clim([cmin cmax]);
    xticks([1:32]);xticklabels([labels]);xtickangle(90);
    title([states4names{c} ': Lassoglm'],'Color',condicolors(c,:));
    subtitle(['coefficients of ' num2str(idxLambda1SE) 'th Lambda' ...
         '(Deviance=' num2str(round(FitInfo.Deviance(Lambda_select),1)) ')']);
    grid on;
end
sgtitle({['Lassoglm: sum-EEG(-500ms) -> H-int'],char(datetime('now'))});
set(gcf,'color','w'); % set background white for copying in ubuntu

% 3 states - combine uncouple state and leading state and call it "independent"
states3names={'Independent','Following','Mutual'};
Inds3={[Inds4_LR(:,1);Inds4_LR(:,2)],Inds4_LR(:,3),Inds4_LR(:,4)};
canvas(0.6,0.2);
cmin=-8e-2;cmax=8e-2;
for c=1:3 % three states
    clear B FitInfo idxLambdaMinDeviance idxLambda1SE Lambda_coef
    [B,FitInfo]  = lassoglm(pow5forpls3(Inds3{c},:),H_all_LR(Inds3{c}),'normal','CV',5);%normal or gamma
    % locate the point with minimum cross-validation error plus one standard deviation
    idxLambda1SE = FitInfo.Index1SE;
    min1coefs = find(B(:,idxLambda1SE));
    Lambda_select=idxLambda1SE;
    if isempty(min1coefs)
       % locate the Lambda with minimum cross-validation error
        idxLambdaMinDeviance = FitInfo.IndexMinDeviance;
        mincoefs = find(B(:,idxLambdaMinDeviance));
        Lambda_select=idxLambdaMinDeviance;
    end
    % reshape from 160 x 1 back to 5 x 32 (freq x chan)
    % Lambda_coef = reshape(B(:,max(idxLambdaMinDeviance)),5,32); 
    Lambda_coef = reshape(B(:,Lambda_select),5,32);
    subplot(1,3,c);
    imagesc(Lambda_coef);colorbar;
    yticks([1:5]);yticklabels({'Delta','Theta','Alpha','Beta','Gamma'});
    set(gca, 'YDir','normal');
    colormap('jet'); clim([cmin cmax]);
    xticks([1:32]);xticklabels([labels]);xtickangle(90);
    title([states3names{c} ': Lassoglm'],'Color',condicolors(c,:));
    subtitle(['coefficients of ' num2str(max(idxLambda1SE)) 'th Lambda' ...
         '(Deviance=' num2str(round(FitInfo.Deviance(Lambda_select),1)) ')']);
    grid on;
end
sgtitle({['Lassoglm: sum-EEG(-500ms) -> H-int'],char(datetime('now'))});

%% PLOT 17-1 Lasso (H-EEG(-500ms) -> H-int)
% H-EEG -> H-int
% (ALL states: 5freq x 32chan = 160 predictors x 192 trials)
clear B FitInfo
[B,FitInfo]  = lassoglm(H5forpls3,H_all_LR,'normal','CV',5);% try normal or gamma
% Find the coefficient vector for the 75th Lambda value in B.
imagesc(B);colorbar;
% Examine the cross-validation plot to see the effect of the Lambda regularization parameter.
lassoPlot(B,FitInfo,'plottype','CV'); 
legend('show') % Show legend
% locate the point with minimum cross-validation error plus one standard deviation
idxLambda1SE = FitInfo.Index1SE;
min1coefs = find(B(:,idxLambda1SE));
Lambda_select=idxLambda1SE;
if isempty(min1coefs)
   % locate the Lambda with minimum cross-validation error
    idxLambdaMinDeviance = FitInfo.IndexMinDeviance;
    mincoefs = find(B(:,idxLambdaMinDeviance));
    Lambda_select=idxLambdaMinDeviance;
end
% reshape from 160 x 1 back to 5 x 32 (freq x chan)
Lambda_coef = reshape(B(:,Lambda_select),5,32); 
cmin=-8e-3;cmax=8e-3;
canvas(0.2,0.2)
imagesc(Lambda_coef);colorbar;
yticks([1:5]);yticklabels({'Delta','Theta','Alpha','Beta','Gamma'});
set(gca, 'YDir','normal');
colormap('jet'); clim([cmin cmax]);
xticks([1:32]);xticklabels([labels]);xtickangle(90);
title({['Lassoglm results in all 4 statues: H-EEG(-500ms) -> H-int'],char(datetime('now'))});
subtitle(['coefficients for the ' num2str(Lambda_select) 'th Lambda value' ...
    ' (Deviance=' num2str(round(FitInfo.Deviance(Lambda_select),1)) ')']);
% subtitle(['the coefficients for the ' num2str(max(idxLambda1SE)) 'th Lambda value']);
grid on;

% (Each of the 4 states: 5freq x 32 chan = 160 predictors x 48 trials)
canvas(0.4,0.4);
cmin=-8e-2;cmax=8e-2;
for c=1:4 % three states
    clear B FitInfo idxLambdaMinDeviance idxLambda1SE Lambda_coef
    [B,FitInfo]  = lassoglm(H5forpls3(Inds4_LR(:,c),:),H_all_LR(Inds4_LR(:,c)),'normal','CV',5);%normal or gamma
    % locate the point with minimum cross-validation error plus one standard deviation
    idxLambda1SE = FitInfo.Index1SE;
    min1coefs = find(B(:,idxLambda1SE));
    Lambda_select=idxLambda1SE;
    if isempty(min1coefs)
       % locate the Lambda with minimum cross-validation error
        idxLambdaMinDeviance = FitInfo.IndexMinDeviance;
        mincoefs = find(B(:,idxLambdaMinDeviance));
        Lambda_select=idxLambdaMinDeviance;
    end
    % reshape from 160 x 1 back to 5 x 32 (freq x chan)
    % Lambda_coef = reshape(B(:,max(idxLambdaMinDeviance)),5,32); 
    Lambda_coef = reshape(B(:,Lambda_select),5,32);
    subplot(2,2,c);
    imagesc(Lambda_coef);colorbar;
    yticks([1:5]);yticklabels({'Delta','Theta','Alpha','Beta','Gamma'});
    set(gca, 'YDir','normal');
    colormap('jet'); clim([cmin cmax]);
    xticks([1:32]);xticklabels([labels]);xtickangle(90);
    title([states4names{c} ': Lassoglm'],'Color',condicolors(c,:));
    subtitle(['coefficients of ' num2str(idxLambda1SE) 'th Lambda' ...
         '(Deviance=' num2str(round(FitInfo.Deviance(Lambda_select),1)) ')']);
    grid on;
end
sgtitle({['Lassoglm: H-EEG(-500ms) -> H-int'],char(datetime('now'))});

% 3 states - combine uncouple state and leading state and call it "independent"
states3names={'independent','following','mutual'};
Inds3={[Inds4_LR(:,1);Inds4_LR(:,2)],Inds4_LR(:,3),Inds4_LR(:,4)};
canvas(0.6,0.2);
cmin=-8e-2;cmax=8e-2;
for c=1:3 % four states
    clear B FitInfo idxLambdaMinDeviance idxLambda1SE Lambda_coef
    [B,FitInfo]  = lassoglm(H5forpls3(Inds3{c},:),H_all_LR(Inds3{c}),'normal','CV',5);%normal or gamma
    % locate the point with minimum cross-validation error plus one standard deviation
    idxLambda1SE = FitInfo.Index1SE;
    min1coefs = find(B(:,idxLambda1SE));
    Lambda_select=idxLambda1SE;
    if isempty(min1coefs)
       % locate the Lambda with minimum cross-validation error
        idxLambdaMinDeviance = FitInfo.IndexMinDeviance;
        mincoefs = find(B(:,idxLambdaMinDeviance));
        Lambda_select=idxLambdaMinDeviance;
    end
    % reshape from 160 x 1 back to 5 x 32 (freq x chan)
    Lambda_coef = reshape(B(:,Lambda_select),5,32);
    subplot(1,3,c);
    imagesc(Lambda_coef);colorbar;
    yticks([1:5]);yticklabels({'Delta','Theta','Alpha','Beta','Gamma'});
    set(gca, 'YDir','normal');
    colormap('jet'); clim([cmin cmax]);
    xticks([1:32]);xticklabels([labels]);xtickangle(90);
    title([states3names{c} ': Lassoglm'],'Color',condicolors(c,:));
    subtitle(['coefficients of ' num2str(Lambda_select) 'th Lambda']);
    grid on;
end
sgtitle({['Lassoglm: H-EEG(-500ms) -> H-int'],char(datetime('now'))});

%% SECT 18 Xcorr Peak and Lag(original order)
clear XcorrPeakLag XcorrPeak
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);
sessions={'synch','synco','synch','synco','synch','synco','synch','synco','synch','synco',...
    'synch','synco'};
XcorrPeakLag=zeros(numSes,12);XcorrPeak=zeros(numSes,12);
cd /ssd/zhibin/1overf/
tic
for r=1:numSes;
    clear intervals conditions sortorders
    runid = num2str(seeds(r,:));
    path = [runid '_2P/Cleaned_data/'];
    load([path  'clean_' runid '.mat']);
    for j = 1:12
        clear intL_good_dmean intR_good_dmean
        % remove the mean
        intL_good_dmean=intervals{j}(:,1)-mean(intervals{j}(:,1));
        intR_good_dmean=intervals{j}(:,2)-mean(intervals{j}(:,2));
        % Xcorr based on int_dmean_drm (before d removal)********
        r12=[];lags12=[];
        [r12,lags12]=xcorr(intL_good_dmean,intR_good_dmean,10,'normalized');
        XcorrPeakLag(r,j)=lags12(find(r12==max(r12)));
        XcorrPeak(r,j)=max(r12);
    end
end
toc

% Organize XcorrPeak for corr
XcorrPeak; % (numSesx12) for all sessions 
% squeeze into 1 vector from the 96 blocks for both subject, for corr with pow in each chan
XcorrPeak_all=reshape(XcorrPeak',[],1); % 96x1 (each element from one block in time sequence) 
% then republicate the the vector for the Left and right
XcorrPeak_all_LR=[XcorrPeak_all;XcorrPeak_all];

% Organize XcorrPeakLag for corr
XcorrPeakLag; % (numSesx12) for all sessions 
% squeeze into 1 vector from the 96 blocks for both subject, for corr with pow in each chan
XcorrPeakLag_all=reshape(XcorrPeakLag',[],1); % 96x1 (each element from one block in time sequence) 
% then republicate the the vector for the Left and right
XcorrPeakLag_all_LR=[XcorrPeakLag_all;XcorrPeakLag_all];

% Organize 
Xcorr10Lag; % 12 numses x 12 blocks x 21 lags % from SECT 9 (matched int. sorted order)
% Just for examing them
% We need original order for PLS 
% % Organized in the same order as H_all_LR as in PLOT13 
% Xcorr_L1=Xcorr10Lag(:,:,10);
% Xcorr_R1=Xcorr10Lag(:,:,12);
% % squeeze into 1 vector from the 144 blocks for each subject, for corr with pow in each chan
% Xcorr_L1=reshape(Xcorr_L1',[],1);% 144x1 (each element from one block in time sequence) 
% Xcorr_R1=reshape(Xcorr_R1',[],1);
% Xcorr_L1R1=[Xcorr_L1; Xcorr_R1]; 
% % We need original order for PLS

%% PLOT 18-1 PLS (sum-EEG -> XcorrPeak / XcorrPeakLag)
% sum-EEG Power from PLOT 16
% pow5forpls3; % 192 blocks x160 elements extracted colums-wise
% (colums: delta chan1, theta chan 1, alpha chan1, beta chan1, gamma chan1, delta chan2...)
% Re-combine L and R sum-EEG power for 64 predictors (32 chan L + 32 chan R) in PLS 
delta_LR_chan64=[delta_L_chan delta_R_chan]; % (96 x 64) from PLOT 13
theta_LR_chan64=[theta_L_chan theta_R_chan];
alpha_LR_chan64=[alpha_L_chan alpha_R_chan];
beta_LR_chan64=[beta_L_chan beta_R_chan];
gamma_LR_chan64=[gamma_L_chan gamma_R_chan];
% fix the scale in the data
delta_LR_chan64 = delta_LR_chan64./(ones(96,1)*std(delta_LR_chan64));
theta_LR_chan64 = theta_LR_chan64./(ones(96,1)*std(theta_LR_chan64));
alpha_LR_chan64 = alpha_LR_chan64./(ones(96,1)*std(alpha_LR_chan64));
beta_LR_chan64 = beta_LR_chan64./(ones(96,1)*std(beta_LR_chan64));
gamma_LR_chan64 = gamma_LR_chan64./(ones(96,1)*std(gamma_LR_chan64));
% combine 5 bands in a 3d matrix  "trials x chan x  freq" (144 x 64 x 5)
Band5_LR_chan64=cat(3,delta_LR_chan64,theta_LR_chan64,...
    alpha_LR_chan64,beta_LR_chan64,gamma_LR_chan64);
% switch dismenstion to "freq x chan x trials" (5 x 64 x 144)
Band5_LR_chan64_2=permute(Band5_LR_chan64,[3 2 1]);
% reshape into 96 trials x ( 5 freq x 64 chan)
% stack them in 2-d "trials x (freq x chan)"  (144 x 320)
Band5_LR_chan64_3=reshape(Band5_LR_chan64_2,5*64,144)';% 144x320 elements extracted colums-wise
% (colums: delta chan1, theta chan 1, alpha chan1, beta chan1, gamma chan1, delta chan2...)

% XcorrPeak_all_LR; % (288 x 1) from SECT18
XcorrPeak_all; % (144 x 1) from SECT18
% XcorrPeakLag_all_LR; % (288 x 1) from SECT18
XcorrPeakLag_all; % (144 x 1) from SECT18
% cobine 2 depdent measures in 2d matrix
XcorrPeakmat2=[XcorrPeak_all XcorrPeakLag_all];
XcorrPeakmat2names={'Xcorr Peak','Xcorr Peak Lag'};


% Create labels for L and R 
for i=1:32
    labelsL{i}=['L-' labels{i}];
    labelsR{i}=['R-' labels{i}];
end


% select a dependent measure
depen_select=1;
depen_select=2;
% select the scale
cmin=-12e-3;cmax=12e-3;
cmin=-0.07;cmax=0.07;


% (ALL states: 5freq x 32chan = 160 predictors x 96 trials)
% similar code as in PLOT 16
clear plsmodel;
canvas(0.5,0.3);
R2=[];reg=[];ypred=[];
[R2,reg,ypred] = npls_pred(Band5_LR_chan64_3,XcorrPeakmat2(:,depen_select),1);
% reshape from 320 x 1 back to 5 x 64 (freq x chan)
plsmodel(depen_select).weights = reshape(reg{1},5,64); 
plsmodel(depen_select).R2 = R2;
plsmodel(depen_select).ypred = ypred;
imagesc(plsmodel(depen_select).weights);colorbar; % caxis([-2 2]*10E-7); % by default, imagesc reverse the Y 
yticks([1:5]);yticklabels({'Delta','Theta','Alpha','Beta','Gamma'});
set(gca, 'YDir','normal');
colormap('jet'); 
% clim([cmin cmax]);
xticks([1:64]);xticklabels([labelsL labelsR]);xtickangle(90);grid on;
title({['PLS model (R2= ' num2str(round(R2,1)) ') in all 4 statues: sum-EEG(-500ms) -> '...
XcorrPeakmat2names{depen_select} ' ^{* PLOT 18-1}'],char(datetime('now'))});

% (Each of the 4 states: 5freq x 32 chan = 160 predictors x 48 trials)
% similar code as in PLOT 16
canvas(0.4,0.4);
plsmodel=[];
cmin=-3e-7;cmax=3e-7;
% cmin=-0.07;cmax=0.07;
Fac=1;
for c=1:4 % four states
    R2=[];reg=[];ypred=[];
    [R2,reg,ypred] = npls_pred(Band5_LR_chan64_3(Inds4_LR(1:24,c),:),XcorrPeakmat2(Inds4_LR(1:24,c),depen_select),Fac);
    plsmodel(c).weights=reshape(reg{Fac},5,64); 
    plsmodel(c).R2 = R2;
    plsmodel(c).ypred = ypred;
    subplot(2,2,c);
    imagesc(plsmodel(c).weights);colorbar; 
    yticks([1:5]);yticklabels({'Delta','Theta','Alpha','Beta','Gamma'});
    set(gca, 'YDir','normal');
    xticks([1:64]);xticklabels([labelsL labelsR]);xtickangle(90);grid on; % from PLOT 18-1 
    colormap('jet'); 
    % clim([cmin cmax]);% caxis([-0.1 0.1]);
    title([states4names{c} ': PLS model (R2= ' num2str(round(R2,1)) '  Fac= ' num2str(Fac) ') '], ...
        'Color',condicolors(c,:));
    grid on;
end
sgtitle({['PLS model: sum-EEG(-500ms) ->' XcorrPeakmat2names{depen_select} ' ^{* PLOT 18-1}'],char(datetime('now'))});
% topoplot for uncouple and mutual
addpath /home/zhibinz2/Documents/GitHub/matlab-archive/hnlcode/common/gen_code/color
hnc = hotncold(100);
band5names={'Delta','Theta','Alpha','Beta','Gamma'};
states2names={'Uncouple','Mutual'};
states4names;
figure; % canvas(0.9,0.5)
cmin=-8e-3;cmax=8e-3;
c=[1 4];
for s=1:2
    for b=1:5
        subplot(2,5,(s-1)*5+b)
        topoplot(plsmodel(c(s)).weights(b,:),channels,'nosedir','+X','style','map');colorbar; 
        colormap(hnc); % colormap('jet'); %  clim([cmin cmax]);
        if b==3;title({['PLS model: sum-EEG(-500ms) -> ' XcorrPeakmat2names{depen_select} ' ^{* PLOT 18-1}'], ...
                [states4names{c(s)} '(R2= ' num2str(round(plsmodel(c(s)).R2,1)) ...
            '  Fac= ' num2str(Fac) ') ']},'Color',condicolors(c(s),:));end
        subtitle(band5names(b));
    end
end
set(gcf,'color','w'); % set background white for copying in ubuntu

%% SECT 19 ...(unfinished): try out DFA on EEG
ans(:,15)
[D,Alpha1,n,F_n,FitValues]=DFA_main(ans(:,15));
%% SECT 20 ...(unfinished): continuous EEG 5 bands
% tryout in one session
s=1; b=1;
zEEG_L; zEEG_R; %Refter to SECT 12
EEG=zEEG_L{b};

tic
[deltaEEG{s,b}, thetaEEG{s,b}, alphaEEG{s,b}, betaEEG{s,b}, gammaEEG{s,b}]...
                =filter5band(EEG,sr);
toc
% processing time for filtering one block of EEG into 5 bands take 30s
% for 192 blocks, it is going to take 96 min

% do it in all sessions
clear
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816];
sessions={'synch','synco','synch','synco','synch','synco','synch','synco'};
cd /ssd/zhibin/1overf/

tic 
for s=1:8 % each session
    runid=num2str(seeds(s,:));
    clear dataL dataR
    load(['/ssd/zhibin/1overf/' runid '_2P/Cleaned_data/clean_' runid '.mat' ])
    zEEG_L=[];zEEG_R=[];
    deltaEEG=cell(2,12); thetaEEG=cell(2,12); alphaEEG=cell(2,12); betaEEG=cell(2,12); gammaEEG=cell(2,12); 
    for b=1:12 % each block  
        zEEG_L=zscore(dataL{b}(:,1:32),[],1);
        zEEG_R=zscore(dataR{b}(:,1:32),[],1);
        [deltaEEG{1,b}, thetaEEG{1,b}, alphaEEG{1,b}, betaEEG{1,b}, gammaEEG{1,b}]...
                =filter5band(zEEG_L,sr);
        [deltaEEG{2,b}, thetaEEG{2,b}, alphaEEG{2,b}, betaEEG{2,b}, gammaEEG{2,b}]...
                =filter5band(zEEG_L,sr);
    end
    cd(['/ssd/zhibin/1overf/' runid '_2P/Cleaned_data/']);
    save([runid '_5band.mat'],'deltaEEG','thetaEEG','alphaEEG','betaEEG','gammaEEG');
    cd /ssd/zhibin/1overf/
end
toc
% take about 2 hours to filter all sessions

% load one session
clear
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816];
sessions={'synch','synco','synch','synco','synch','synco','synch','synco'};
cd /ssd/zhibin/1overf/
% select one session
s=1;
runid=num2str(seeds(s,:));
load(['/ssd/zhibin/1overf/' runid '_2P/Cleaned_data/' runid '_5band.mat'])

figure;
subplot(2,1,1);
EEG=deltaEEG{1,11};
plot(EEG);
subplot(2,1,2);

% try Re reference
addpath /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing
% This function pass filtered_data and get back average referenced data.
[reRef_data] = reRef2(EEG);
plot(reRef_data);


% compute envelope using hilbert

% savitzy golay filter (smooth window 2000/f. f is the mid freq)

% average trial power in each band, corr with H-int

% downsample to 100 Hz

% H-EEG in each band, corr with H-int

%% SECT 21 GC: arfit and MVGC
cd /ssd/zhibin/1overf/all_session20220713_1005
load('int_dmean_drm.mat')
int_dmean_drm; % from SECT9 (matched int) (sorted order)
% 2 L/R X 12 blocks x 12 sessions

% orgainized the synchronization trials into 4 conditions
condi4Ind={[1:3],[4:6],[7:9],[10:12]};
condi4names={'Uncouple','L-lead','R-lead','Mutual'};
syn2Ind={[1:2:11],[2:2:12]};syn2names={'Synch','Synco'};

% test_data=int_dmean_drm(:,1:3,1:2:11); % uncouple synch
% test_data=int_dmean_drm(:,4:6,1:2:11); % L-lead synch
% test_data=int_dmean_drm(:,7:9,1:2:11); % R-lead synch
% test_data=int_dmean_drm(:,10:12,1:2:11); % R-lead synch

% arfit
% organize into time x 2 L/R matrix then apply arfit and plot GC
pmin=1; % minimal order
pmax=10; % maximal order
figure;
for syn=1:2
    for condi=1:4
        test_data=int_dmean_drm(:,condi4Ind{condi},syn2Ind{syn});
        % Method1a: concatenate into one big trial
        Int12LR=[];
        for b=1:size(test_data,2)
            for s=1:size(test_data,3)
                    Int12LR=[Int12LR; [test_data{1,b,s} test_data{2,b,s}]];
            end
        end
%         % Method1b: randomized the concatenation
%         Int12LR=[];
%         for b=randperm(size(test_data,2))
%             for s=randperm(size(test_data,3))
%                     Int12LR=[Int12LR; [test_data{1,b,s} test_data{2,b,s}]];
%             end
%         end
        if syn==1; subplot(2,4,condi);end
        if syn==2; subplot(2,4,4+condi);end
        [w, A, C, sbc, fpe, th]=arfit(Int12LR, pmin, pmax);
        % Granger Causality
        % GC of L->R
        GC_L2R=log(C(2,2)/C(1,2));
        % GC of R->L
        GC_R2L=log(C(1,1)/C(1,2));
        % plot
        bar([GC_L2R GC_R2L]);xticks([1 2]);xticklabels({'L->R','R->L'});
        title([condi4names{condi}]);subtitle(syn2names{syn});
        ylim([0 6]);
    end
end
sgtitle({['GC: arfit (concatenate data)'],char(datetime('now'))});
% sgtitle({['GC: arfit (permuted concatenate data)',char(datetime('now'))});

% MVGC
% organize into time x 2 L/R matrix then apply MVGC and plot
Fs=cell(2,4);
for syn=1:2
    for condi=1:4
        test_data=int_dmean_drm(:,condi4Ind{condi},syn2Ind{syn});
        % Method1a: concatenate into one big trial
        Int12LR=[];
        for b=1:size(test_data,2) % block
            for s=1:size(test_data,3) % session
                    Int12LR=[Int12LR; [test_data{1,b,s} test_data{2,b,s}]];
            end
        end
%         % Method1b: randomized the concatenation
%         Int12LR=[];
%         for b=randperm(size(test_data,2))
%             for s=randperm(size(test_data,3))
%                     Int12LR=[Int12LR; [test_data{1,b,s} test_data{2,b,s}]];
%             end
%         end
%         if syn==1; subplot(2,4,condi);end
%         if syn==2; subplot(2,4,4+condi);end
        X=[];
        X=permute(Int12LR,[2,1]);
        % Parameters
        ntrials   = size(X,3);     % number of trials
        nobs      = size(X,2);  % number of observations per trial
        regmode   = 'LWR';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
        icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default) 
        morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
        momax     = 10;     % maximum model order for model order estimation
        % recomended acmaxlags minimum = 105
        acmaxlags = 10;   % maximum autocovariance lags (empty for automatic calculation)
        tstat     = '';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
        alpha     = 0.05;   % significance level for significance test
        mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')
        % for plotting only
        fs        = 1.3;    % sample rate (Hz)
        fres      = [];     % frequency resolution (empty for automatic calculation)
        nvars = size(X,1); % number of variables
        % model order estimation
        [AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X,momax,icregmode);
        % Select model order.
        if     strcmpi(morder,'actual')
            morder = amo;
            fprintf('\nusing actual model order = %d\n',morder);
        elseif strcmpi(morder,'AIC')
            morder = moAIC;
            fprintf('\nusing AIC best model order = %d\n',morder);
        elseif strcmpi(morder,'BIC')
            morder = moBIC;
            fprintf('\nusing BIC best model order = %d\n',morder);
        else
            fprintf('\nusing specified model order = %d\n',morder);
        end
        % VAR model estimation 
        [A,SIG] = tsdata_to_var(X,morder,regmode);
        % Check for failed regression
        assert(~isbad(A),'VAR estimation failed');
        % Autocovariance calculation
        [G,info] = var_to_autocov(A,SIG,acmaxlags);
        var_acinfo(info,true); % report results (and bail out on error)
        % Granger causality calculation: time domain 
        F = autocov_to_pwcgc(G);
        % save F
        % F(2,1): L->R; F(1,2): R->L;
        Fs{syn,condi}=[F(2,1) F(1,2)];
    end
end

% 3 subplots 
if true % 3 subplots
ymax=0.06;
canvas(0.4, 0.4);
tiledlayout(1,3);
% Tile 1
nexttile
bar(mean([Fs{1,1};Fs{2,1}],'all'));xticks([1]);xticklabels({'L<->R'});
title('uncoupled');ylabel('GC');
ylim([0 ymax]);% xlim([0.25 1.75]);
% Tile 2
nexttile
bar([mean([Fs{1,2}(1);Fs{2,2}(1);Fs{1,3}(2);Fs{2,3}(2)]) ...
    mean([Fs{1,2}(2);Fs{2,2}(2);Fs{1,3}(1);Fs{2,3}(1)])]);xticks([1 2]);
xticklabels({'Leader->Follower','Follower->Leader'});
title('unidirectional');ylabel('GC');
ylim([0 ymax]);% xlim([0.25 2.75]);
% Tile 1
nexttile
bar(mean([Fs{1,4};Fs{2,4}],'all'));xticks([1]);xticklabels({'L<->R'});
title('bidirectional');ylabel('GC');
ylim([0 ymax]);% xlim([0.25 1.75]);
% sgtitle('GC: var-to-autocov (concatenate data)')
% sgtitle('GC: var-to-autocov (permuted concatenate data)')
sg=annotation('textbox',[0.3 0.01 0.4 0.05],'string',...
    {['MVGC (concatenate data) H-int ^{* SECT 21}'],char(datetime('now'))});
end

% 3 subplots combined as one 
for collapse=1;
ymax=0.06;
figure;
bar([mean([Fs{1,1};Fs{2,1}],'all')...
    mean([Fs{1,2}(1);Fs{2,2}(1);Fs{1,3}(2);Fs{2,3}(2)]) ...
    mean([Fs{1,2}(2);Fs{2,2}(2);Fs{1,3}(1);Fs{2,3}(1)]) ...
    mean([Fs{1,4};Fs{2,4}],'all')]);
xticks([1 2 3 4]);
xticklabels({'uncouple (L<->R)','Leader->Follower','Follower->Leader', 'mutual (L<->R)'});
ylim([0 ymax]);ylabel('GC');
sg=annotation('textbox',[0.3 0.01 0.4 0.05],'string',...
    'MVGC (concatenate data) H-int ^{* SECT 21}')
end

% organize into time x 2 L/R matrix then apply MVGC and collect p values (1x3 suplots)
% and plot
Fs=nan(1,4);Ps=nan(1,4);
for condi=1
    test_data=[];
    test_data=int_dmean_drm(:,condi4Ind{condi},:);
    L_mat=[];R_mat=[];
    L_mat=cell2mat(reshape(squeeze(test_data(1,:,:)),36,1));
    R_mat=cell2mat(reshape(squeeze(test_data(2,:,:)),36,1));
    X=[];
    X=[L_mat R_mat; R_mat L_mat];
    X=permute(X,[2,1]);
    [F,pval,sig] = myGCautocov(X);
    Fs(1)=F(2,1);
    Ps(1)=pval(2,1);
end
for condi=[2 3]
    test_data=[];
    test_data=int_dmean_drm(:,condi4Ind{2},:);
    L_mat=[];R_mat=[];
    L_mat=cell2mat(reshape(squeeze(test_data(1,:,:)),36,1));
    R_mat=cell2mat(reshape(squeeze(test_data(2,:,:)),36,1));
    X1=[];
    X1=[L_mat R_mat];
    test_data=[];
    test_data=int_dmean_drm(:,condi4Ind{3},:);
    L_mat=[];R_mat=[];
    L_mat=cell2mat(reshape(squeeze(test_data(1,:,:)),36,1));
    R_mat=cell2mat(reshape(squeeze(test_data(2,:,:)),36,1));
    X2=[];
    X2=[R_mat L_mat];
    X=[];
    X=[X1; X2];
    X=permute(X,[2,1]);
    [F,pval,sig] = myGCautocov(X);
    Fs(2)=F(2,1);Fs(3)=F(1,2);
    Ps(2)=pval(2,1);Ps(3)=pval(1,2);
end
for condi=4
    test_data=[];
    test_data=int_dmean_drm(:,condi4Ind{condi},:);
    L_mat=[];R_mat=[];
    L_mat=cell2mat(reshape(squeeze(test_data(1,:,:)),36,1));
    R_mat=cell2mat(reshape(squeeze(test_data(2,:,:)),36,1));
    X=[];
    X=[L_mat R_mat; R_mat L_mat];
    X=permute(X,[2,1]);
    [F,pval,sig] = myGCautocov(X);
    Fs(4)=F(2,1);
    Ps(4)=pval(2,1);
end
% 3 subplots 
for subplot3=1
ymax=0.04;
canvas(0.5, 0.45);
tiledlayout(1,3);
% Tile 1
nexttile
bar(Fs(1));xticks([1]);xticklabels({'L<->R'});
title('Uncoupled');ylabel('GC');
ylim([0 ymax]);% xlim([0.25 1.75]);
text(1-0.2, Fs(1)+0.005, sprintf('p=%.2f',Ps(1)),'Color',[0 0 0]);
% Tile 2
nexttile
bar(Fs(2:3));xticks([1 2]);
xticklabels({'Leader->Follower','Follower->Leader'});
title('Unidirectional');ylabel('GC');
ylim([0 ymax]);% xlim([0.25 2.75]);
text(1-0.2, Fs(2)+0.005, sprintf('p=%.2f',Ps(2)),'Color',[1 0 0]);
text(2-0.2, Fs(3)+0.005, sprintf('p=%.2f',Ps(3)),'Color',[0 0 0]);
% Tile 1
nexttile
bar(Fs(4));xticks([1]);xticklabels({'L<->R'});
title('Bidirectional');ylabel('GC');
ylim([0 ymax]);% xlim([0.25 1.75]);
text(1-0.2, Fs(4)+0.005, sprintf('p=%.2f',Ps(4)),'Color',[1 0 0]);
set(gcf,'color','w'); % set background white for copying in ubuntu
% sgtitle('GC: var-to-autocov (concatenate data)')
% sgtitle('GC: var-to-autocov (permuted concatenate data)')
delete(sg)
sg=annotation('textbox',[0.3 0.01 0.4 0.07],'string',...
    {['MVGC (concatenate data) H-int ^{* SECT 21}'],char(datetime('now'))});
end


% organize into time x 2 L/R matrix then apply MVGC and collect p values (2 suplots separate synch/o)
Fs=nan(2,4);Ps=nan(2,4); % 2 syn type x 4 conditions
if true % calculation of GC and P values
for condi=1
    for t=1:2
    test_data=[];
    test_data=int_dmean_drm(:,condi4Ind{condi},syn2Ind{t});
    L_mat=[];R_mat=[];
    L_mat=cell2mat(reshape(squeeze(test_data(1,:,:)),18,1));
    R_mat=cell2mat(reshape(squeeze(test_data(2,:,:)),18,1));
    X=[];
    X=[L_mat R_mat; R_mat L_mat];
    X=permute(X,[2,1]);
    [F,pval,sig] = myGCautocov(X);
    Fs(t,1)=F(2,1);
    Ps(t,1)=pval(2,1);
    end
end

for condi=[2 3]
    for t=1:2
    test_data=[];
    test_data=int_dmean_drm(:,condi4Ind{2},syn2Ind{t});
    L_mat=[];R_mat=[];
    L_mat=cell2mat(reshape(squeeze(test_data(1,:,:)),18,1));
    R_mat=cell2mat(reshape(squeeze(test_data(2,:,:)),18,1));
    X1=[];
    X1=[L_mat R_mat];

    test_data=[];
    test_data=int_dmean_drm(:,condi4Ind{3},syn2Ind{t});
    L_mat=[];R_mat=[];
    L_mat=cell2mat(reshape(squeeze(test_data(1,:,:)),18,1));
    R_mat=cell2mat(reshape(squeeze(test_data(2,:,:)),18,1));
    X2=[];
    X2=[R_mat L_mat];
    
    X=[];
    X=[X1; X2];
    X=permute(X,[2,1]);
    [F,pval,sig] = myGCautocov(X);
    Fs(t,2)=F(2,1);Fs(t,3)=F(1,2);
    Ps(t,2)=pval(2,1);Ps(t,3)=pval(1,2);
    end
end

for condi=4
    for t=1:2
    test_data=[];
    test_data=int_dmean_drm(:,condi4Ind{condi},syn2Ind{t});
    L_mat=[];R_mat=[];
    L_mat=cell2mat(reshape(squeeze(test_data(1,:,:)),18,1));
    R_mat=cell2mat(reshape(squeeze(test_data(2,:,:)),18,1));
    X=[];
    X=[L_mat R_mat; R_mat L_mat];
    X=permute(X,[2,1]);
    [F,pval,sig] = myGCautocov(X);
    Fs(t,4)=F(2,1);
    Ps(t,4)=pval(2,1);
    end
end
end

% 1 barplots 
if true % barplot 4groups x 2types
canvas(0.23, 0.4);
model_series = Fs';
clear b
b = bar(model_series, 'FaceColor','flat');
b(1).FaceColor=darkgreen;b(2).FaceColor=pink;
xticks(1:4);xticklabels({'Uncoupled','Leader -> Follower','Follower -> Leader','Mutual'});
xl = get(gca,'XTickLabel');  
set(gca,'XTickLabel',xl,'fontsize',20,'FontWeight','bold')
% set(gca,'XTickLabel',{'Uncoupled','Leading','Following','Mutual'},'fontsize',12,'FontWeight','bold');
% label p values
% delete(findall(gcf,'type','text'))
for t=1:2
    for p=1:4
        if Ps(t,p) <0.01; signifColor=red; else; signifColor=black; end
        if t==1
            % text(p-0.27, Fs(t,p)+0.005, sprintf('p=%.2f',Ps(t,p)),'Color',signifColor,'FontSize',16);
        else
            text(p+0.01, Fs(t,p)+0.005, sprintf('p=%.2f',Ps(t,p)),'Color',signifColor,'FontSize',16);
        end
    end
end
legend({'Synch','Synco'},'location','northwest','FontSize', 30);
ylabel('GC','FontSize',17); 
yl = get(gca,'YTickLabel');  
set(gca,'YTickLabel',yl,'fontsize',20,'FontWeight','bold');
% delete(sg)
% sg=annotation('textbox',[0.2 0.01 0.7 0.07],'string',...
%     {['MVGC (concatenate data) H-int ^{* SECT 21}'],char(datetime('now'))});
set(gcf,'color','w'); % set background white for copying in ubuntu
end


% 1 barplots (removed p value)
if true % barplot 4groups x 2types
canvas(0.23, 0.4);
model_series = Fs';
clear b
b = bar(model_series, 'FaceColor','flat');
b(1).FaceColor=darkgreen;b(2).FaceColor=pink;
xticks(1:4);xticklabels({'Uncoupled','Leader -> Follower','Follower -> Leader','Mutual'});
xl = get(gca,'XTickLabel');  
set(gca,'XTickLabel',xl,'fontsize',20,'FontWeight','bold')
% set(gca,'XTickLabel',{'Uncoupled','Leading','Following','Mutual'},'fontsize',12,'FontWeight','bold');
% label p values
% delete(findall(gcf,'type','text'))
% for t=1:2
%     for p=1:4
%         if Ps(t,p) <0.01; signifColor=red; else; signifColor=black; end
%         if t==1
%             text(p-0.27, Fs(t,p)+0.005, sprintf('*','Color',signifColor,'FontSize',16))
%             % text(p-0.27, Fs(t,p)+0.005, sprintf('p=%.2f',Ps(t,p)),'Color',signifColor,'FontSize',16);
%         else
%             text(p+0.01, Fs(t,p)+0.005, sprintf('*','Color',signifColor,'FontSize',16))
%             % text(p+0.01, Fs(t,p)+0.005, sprintf('p=%.2f',Ps(t,p)),'Color',signifColor,'FontSize',16);
%         end
%     end
% end
legend({'Synch','Synco'},'location','northwest','FontSize', 30);
ylabel('GC','FontSize',17); 
yl = get(gca,'YTickLabel');  
set(gca,'YTickLabel',yl,'fontsize',20,'FontWeight','bold');
% delete(sg)
% sg=annotation('textbox',[0.2 0.01 0.7 0.07],'string',...
%     {['MVGC (concatenate data) H-int ^{* SECT 21}'],char(datetime('now'))});
set(gcf,'color','w'); % set background white for copying in ubuntu
end

%% PLOT 21 EEG (-500ms) -> GC (MVGC): -Corr -PLS
% Organize GC for PLS
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);
sessions={'synch','synco','synch','synco','synch','synco','synch','synco','synch','synco',...
    'synch','synco'};
% organize conditions in all sessions into a matrix of 8 x 12
condition_all_mat=NaN(numSes,12); % (Refer to SECT 10-1)
% Compute GC_all for all trials for PLS
Fs=nan(2,numSes,12);
for s=1:numSes;
    clear intervals conditions
    runid=num2str(seeds(s,:));
    load(['/ssd/zhibin/1overf/' runid '_2P/Cleaned_data/clean_' runid '.mat' ]);
    condition_all_mat(s,:)=conditions;
    intL_good_dmean=[];intR_good_dmean=[];
    for b=1:12
        % remove the mean (no need for DFA)
        intL_good_dmean=intervals{b}(:,1)-mean(intervals{b}(:,1));
        intR_good_dmean=intervals{b}(:,2)-mean(intervals{b}(:,2));
        X=[];F=[];
        X=[intL_good_dmean intR_good_dmean];
        X=permute(X,[2,1]);
        [F,pval,sig] = myGCautocov(X);
        Fs(1,s,b)=F(2,1);
        X=[];F=[];
        X=[intR_good_dmean intL_good_dmean];
        X=permute(X,[2,1]);
        [F,pval,sig] = myGCautocov(X);
        Fs(2,s,b)=F(2,1);
    end
end
% Examine Fs in sorted order
[condition_all_sorted,ind2]=sort(condition_all_mat,2)
canvas(0.8,0.8);
for s=1:12
    for b=1:12
        subplot(12,12,(s-1)*12+b);
        clear barplot
        barplot=bar([Fs(1,s,ind2(s,b)) Fs(2,s,ind2(s,b))]);
        set(gca, 'XTickLabel', {'L -> R' 'R -> L'});
        barplot.FaceColor='flat';barplot.CData(1,:) = [1 0 0];barplot.CData(2,:) = [0 0 1];
        xlim([0.25 2.75]); ylim([0 1]);
        ylabel('GC');
        title(condi4names(condition_all_sorted(s,b)),'color',condicolors(condition_all_sorted(s,b),:))
    end
end
% Examine only uncouple Fs
canvas(0.2,0.8);
for s=1:12
    for b=1:3
        subplot(12,3,(s-1)*3+b);
        clear barplot
        barplot=bar([Fs(1,s,ind2(s,b)) Fs(2,s,ind2(s,b))]);
        set(gca, 'XTickLabel', {'L -> R' 'R -> L'});
        barplot.FaceColor='flat';barplot.CData(1,:) = [1 0 0];barplot.CData(2,:) = [0 0 1];
        xlim([0.25 2.75]); ylim([0 0.125]);
        ylabel('GC');
        title(condi4names(condition_all_sorted(s,b)),'color',condicolors(condition_all_sorted(s,b),:))
    end
end

% Organize Fs for PLS
% same way as we did with H in PLOT 13
Fs; % (2xnumSesx12) for all sessions (matched int) (original order)
Fs_L=squeeze(Fs(1,:,:));
Fs_R=squeeze(Fs(2,:,:));
% squeeze into 1 vector from the 96 blocks for each subject, for corr with pow in each chan
Fs_L=reshape(Fs_L',[],1);% 144x1 (each element from one block in time sequence) 
Fs_R=reshape(Fs_R',[],1);
% Combine L and R
Fs_LR=[Fs_L;Fs_R]; % to be used for corr and PLS

% run SECT 12 and PLOT 13 to organize EEG first (-500ms)
% Corr
% Compute the correlation between sum-EEG pow and BPint_xcorrSeries_LR
% for c=1:32
%     delta_L_F_L_corr(c)=corr(delta_L_chan(:,c),Fs_L);
%     theta_L_F_L_corr(c)=corr(theta_L_chan(:,c),Fs_L);
%     alpha_L_F_L_corr(c)=corr(alpha_L_chan(:,c),Fs_L);
%     beta_L_F_L_corr(c)=corr(beta_L_chan(:,c),Fs_L);
%     gamma_L_F_L_corr(c)=corr(gamma_L_chan(:,c),Fs_L);
%     delta_R_F_R_corr(c)=corr(delta_R_chan(:,c),Fs_R);
%     theta_R_F_R_corr(c)=corr(theta_R_chan(:,c),Fs_R);
%     alpha_R_F_R_corr(c)=corr(alpha_R_chan(:,c),Fs_R);
%     beta_R_F_R_corr(c)=corr(beta_R_chan(:,c),Fs_R);
%     gamma_R_F_R_corr(c)=corr(gamma_R_chan(:,c),Fs_R);
% end
% indices for 4 states from PLOT 10-1
Inds4_LR;
% Compute Correlation in 4 states
% % fix the scale in the data (as in PLOT 16)
% delta_LR_chan = delta_LR_chan./(ones(288,1)*std(delta_LR_chan)); %organized in SECT 13
% theta_LR_chan = theta_LR_chan./(ones(288,1)*std(theta_LR_chan));
% alpha_LR_chan = alpha_LR_chan./(ones(288,1)*std(alpha_LR_chan));
% beta_LR_chan = beta_LR_chan./(ones(288,1)*std(beta_LR_chan));
% gamma_LR_chan = gamma_LR_chan./(ones(288,1)*std(gamma_LR_chan));
delta_LR_F_LR_4corr=zeros(4,32);
theta_LR_F_LR_4corr=zeros(4,32);
alpha_LR_F_LR_4corr=zeros(4,32);
beta_LR_F_LR_4corr=zeros(4,32);
gamma_LR_F_LR_4corr=zeros(4,32);
for s=1:4
    for c=1:32
        delta_LR_F_LR_4corr(s,c)=corr(delta_LR_chan(Inds4_LR(:,s),c),Fs_LR(Inds4_LR(:,s)));
        theta_LR_F_LR_4corr(s,c)=corr(theta_LR_chan(Inds4_LR(:,s),c),Fs_LR(Inds4_LR(:,s)));
        alpha_LR_F_LR_4corr(s,c)=corr(alpha_LR_chan(Inds4_LR(:,s),c),Fs_LR(Inds4_LR(:,s)));
        beta_LR_F_LR_4corr(s,c)=corr(beta_LR_chan(Inds4_LR(:,s),c),Fs_LR(Inds4_LR(:,s)));
        gamma_LR_F_LR_4corr(s,c)=corr(gamma_LR_chan(Inds4_LR(:,s),c),Fs_LR(Inds4_LR(:,s)));
    end
end
% Combine L and R in 4 states(4x5)
canvas(0.3,0.5);
cmin=-0.4;cmax=0.4;
for s=1:4
    subplot(4,5,5*(s-1)+1);
    topoplot(delta_LR_F_LR_4corr(s,:),channels,'nosedir','+X','style','map');
    clim([cmin cmax]);
    subplot(4,5,5*(s-1)+2);
    topoplot(theta_LR_F_LR_4corr(s,:),channels,'nosedir','+X','style','map');
    clim([cmin cmax]);
    subplot(4,5,5*(s-1)+3);
    topoplot(alpha_LR_F_LR_4corr(s,:),channels,'nosedir','+X','style','map');
    clim([cmin cmax]);
    subplot(4,5,5*(s-1)+4);
    topoplot(beta_LR_F_LR_4corr(s,:),channels,'nosedir','+X','style','map');
    clim([cmin cmax]);
    subplot(4,5,5*(s-1)+5);
    topoplot(gamma_LR_F_LR_4corr(s,:),channels,'nosedir','+X','style','map');
    clim([cmin cmax]);
end
colormap(hnc)
% https://www.mathworks.com/help/matlab/ref/matlab.graphics.illustration.colorbar-properties.html
cb=colorbar;
cb.AxisLocation = 'out';
cb.Position = [0.92 0.15 0.01 0.75];
% annotations
delete(findall(gcf,'type','annotation'))
h0=annotation('textbox',[0.17 0.95 0.05 0.03],'string','Delta','color',[0 0 0])
h1=annotation('textbox',[0.33 0.95 0.05 0.03],'string','Theta','color',[0 0 0])
h2=annotation('textbox',[0.5 0.95 0.05 0.03],'string','Alpha','color',[0 0 0])
h3=annotation('textbox',[0.66 0.95 0.05 0.03],'string','Beta','color',[0 0 0])
h4=annotation('textbox',[0.81 0.95 0.05 0.03],'string','Gamma','color',[0 0 0])
v0=annotation('textbox',[0.14 0.15 0.05 0.03],'string','Mutual','color',condicolors(4,:))
v1=annotation('textbox',[0.14 0.37 0.05 0.03],'string','Following','color',condicolors(3,:))
v2=annotation('textbox',[0.14 0.59 0.05 0.03],'string','Leading','color',condicolors(2,:))
v3=annotation('textbox',[0.14 0.81 0.05 0.03],'string','Uncouple','color',condicolors(1,:))
set(v0,'Rotation',90);set(v1,'Rotation',90);set(v2,'Rotation',90);set(v3,'Rotation',90);
delete(sg)
sg=annotation('textbox',[0.3 0.01 0.4 0.07],'string',...
    {['Correlation of sum-EEG (-500ms) and GC ^{* PLOT 21}'],char(datetime('now'))});
set([h0 h1 h2 h3 h4 v0 v1 v2 v3], 'fitboxtotext','on',...
    'edgecolor','none')
set(gcf,'color','w'); % set background white for copying in ubuntu

% PLS
addpath /home/zhibinz2/Documents/GitHub/matlab/ramesh/plsmodel
addpath(genpath('/home/zhibinz2/Documents/GitHub/matlab/external/')); 
addpath /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/PLS
% get pow5forpls3 from PLOT 16
pow5forpls3;
% (ALL states: 5freq x 32chan = 160 predictors x 288 trials)
% -updated to use mynpls_pred function
reg=[];ypred_fit=[];X_MC=[];Y_MC=[];
clear plsmodel;
[reg,ypred_fit,X_MC,Y_MC] = mynpls_pred(pow5forpls3,Fs_LR);
ssOr=sum(sum((Y_MC-mean(Y_MC)).^2))
AIC=nan(1,10);R2s=nan(1,10);ssERs=nan(1,10);
for Fac=1:10;
    R2=[];ssEr=[];
    ssEr=sum(sum((Y_MC-ypred_fit(:,Fac)).^2))
    ssERs(Fac)=ssEr;
    R2=100*(1-ssEr/ssOr)
    R2s(Fac)=R2;
    % Cal AIC
    AIC(Fac)=log(ssEr)+2*(Fac);
end
figure;
subplot(2,1,1);
plot(1:10,AIC,'.','MarkerSize',20);xlabel('nFac');ylabel('AIC');
xlim([0 11]);
subtitle('All states');
subplot(2,1,2);
plot(1:10,R2s,'r.','MarkerSize',20);xlabel('nFac');ylabel('R^2');
xlim([0 11]);
subtitle('All states');
sgtitle({['AIC & R^2 for PLS model in all statues: sum-EEG(-500ms) -> GC ^{ *PLOT 21}'],char(datetime('now'))});
set(gcf,'color','w'); % set background white for copying in ubuntu
% Try cross validation (myxvalidation)
% select the Fac with minimal AIC
Fac=1
% do the #-fold cross validation
ssEr=[];R2=[];
X=pow5forpls3;Y=Fs_LR;
[ssEr,R2]  = myxvalidation(X,Y,Fac,nSet);
% reshape regs{1:4}{Fac} 
clear plsmodel
% reshape from 160 x 1 back to 5 x 32 (freq x chan)
plsmodel.weights=reshape(reg{Fac},5,32); 
canvas(0.2,0.2)
cmin=-0.01;cmax=0.01;
imagesc(plsmodel.weights);colorbar; % caxis([-2 2]*10E-7); % by default, imagesc reverse the Y 
yticks([1:5]);yticklabels({'Delta','Theta','Alpha','Beta','Gamma'});
colormap('jet'); clim([cmin cmax]);
set(gca, 'YDir','normal');
xticks([1:32]);xticklabels([labels]);xtickangle(90);grid on;
title({['PLS model (R^2= ' num2str(round(R2,1)) ') in all statues: ' ...
    'sum-EEG(-500ms) -> GC  ^{* PLOT 21}'],char(datetime('now'))});
set(gcf,'color','w'); % set background white for copying in ubuntu
% (Each of the 4 states: 5freq x 32 chan = 160 predictors x 72
% trials)-updated to use mynpls_pred function
cd /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/PLS
regs=[];AIC=nan(2,10);R2s=nan(4,10);ssERs=nan(4,10);
for c=1:4 % four states
    R2=[];reg=[];ypred=[];
    [reg,ypred_fit,X_MC,Y_MC] = mynpls_pred(pow5forpls3(Inds4_LR(:,c),:),Fs_LR(Inds4_LR(:,c)));
    regs{c}=reg;
    ssOr=sum(sum((Y_MC-mean(Y_MC)).^2))
    for Fac=1:10;
        R2=[];ssEr=[];
        ssEr=sum(sum((Y_MC-ypred_fit(:,Fac)).^2))
        ssERs(c,Fac)=ssEr;
        R2=100*(1-ssEr/ssOr)
        R2s(c,Fac)=R2;
        % Cal AIC
        AIC(c,Fac)=log(ssEr)+2*(Fac);
    end
end
figure;
for c=1:4
    subplot(2,4,c);
    plot(1:10,AIC(c,:),'.','MarkerSize',20);xlabel('nFac');ylabel('AIC');
    ylim([0 17]);xlim([0 11]);
    subtitle(states4names{c},'Color',condicolors(c,:));
    subplot(2,4,4+c);
    plot(1:10,R2s(c,:),'r.','MarkerSize',20);xlabel('nFac');ylabel('R^2');
    ylim([10 110]);xlim([0 11]);
    subtitle(states4names{c},'Color',condicolors(c,:));
end
sgtitle(['AIC & R^2 for PLS model in all 4 statues: sum-EEG(-500ms) -> GC ^{ *PLOT 21}']);
set(gcf,'color','w'); % set background white for copying in ubuntu
% Try cross validation (myxvalidation)
% select the Fac with minimal AIC
Fac=1
ssErs=[];R2s=[];
% do the cross validation
for c=1:4
    ssEr=[];R2=[];
    X=pow5forpls3(Inds4_LR(:,c),:);Y=BPint_xcorrSeries_LR(Inds4_LR(:,c));
    [ssEr,R2]  = myxvalidation(X,Y,Fac,nSet);
    ssErs(c)=ssEr;
    R2s(c)=R2;
end
% reshape regs{1:4}{Fac} 
clear plsmodel
for c=1:4
    plsmodel(c).weights=reshape(regs{c}{Fac},5,32); 
end
% Plot the 4 states
canvas(0.4,0.4);
cmin=-0.005;cmax=0.005;
for c=1:4
    subplot(2,2,c);
    imagesc(plsmodel(c).weights);colorbar; 
    yticks([1:5]);yticklabels({'Delta','Theta','Alpha','Beta','Gamma'});
    set(gca, 'YDir','normal');
    xticks([1:32]);xticklabels([labels]);xtickangle(90);
    colormap('jet'); clim([cmin cmax]);
    % caxis([-0.1 0.1]);
    title([states4names{c} ': PLS model (R^2= ' num2str(round(R2s(c),1)) '  Fac= ' num2str(Fac) ') '], ...
        'Color',condicolors(c,:));
    grid on;
end
sgtitle({['PLS model: sum-EEG(-500ms) -> GC ^{* PLOT 21}'],char(datetime('now'))});
set(gcf,'color','w'); % set background white for copying in ubuntu
colormap(hnc)
% topoplot for the 4 states
addpath /home/zhibinz2/Documents/GitHub/matlab-archive/hnlcode/common/gen_code/color
hnc = hotncold(100);
band5names={'Delta','Theta','Alpha','Beta','Gamma'};
states2names={'Uncouple','Mutual'};
states4names;
% figure;
canvas(0.3,0.8)
c=[1 2 3 4];
for s=1:4
    for b=1:5
        subplot(4,5,(s-1)*5+b)
        topoplot(plsmodel(c(s)).weights(b,:),channels,'nosedir','+X','style','map');
        clim([cmin cmax]);
    end
end
colormap(hnc)
% https://www.mathworks.com/help/matlab/ref/matlab.graphics.illustration.colorbar-properties.html
cb=colorbar;
cb.AxisLocation = 'out';
cb.Position = [0.92 0.15 0.01 0.75];
set(gcf,'color','w'); % set background white for copying in ubuntu
delete(findall(gcf,'type','annotation'))
h0=annotation('textbox',[0.17 0.95 0.05 0.03],'string','Delta','color',[0 0 0])
h1=annotation('textbox',[0.33 0.95 0.05 0.03],'string','Theta','color',[0 0 0])
h2=annotation('textbox',[0.5 0.95 0.05 0.03],'string','Alpha','color',[0 0 0])
h3=annotation('textbox',[0.66 0.95 0.05 0.03],'string','Beta','color',[0 0 0])
h4=annotation('textbox',[0.82 0.95 0.05 0.03],'string','Gamma','color',[0 0 0])
v0=annotation('textbox',[0.1 0.15 0.05 0.03],'string',...
    {[states4names{c(4)}], [' (R^2= ' num2str(round(R2s(c(4)),1)) ...
            '  Fac= ' num2str(Fac) ') ']}...
            ,'color',condicolors(4,:));
set(v0,'Rotation',90);
v1=annotation('textbox',[0.1 0.35 0.05 0.03],'string',...
    {[states4names{c(3)}], [' (R^2= ' num2str(round(R2s(c(3)),1)) ...
            '  Fac= ' num2str(Fac) ') ']} ...
            ,'color',condicolors(3,:));
set(v1,'Rotation',90);
v2=annotation('textbox',[0.1 0.56 0.05 0.03],'string',...
    {[states4names{c(2)}], [' (R^2= ' num2str(round(R2s(c(2)),1)) ...
            '  Fac= ' num2str(Fac) ') ']} ...
            ,'color',condicolors(2,:));
set(v2,'Rotation',90);
v3=annotation('textbox',[0.1 0.79 0.05 0.03],'string',...
    {[states4names{c(1)}], [' (R^2= ' num2str(round(R2s(c(1)),1)) ...
            '  Fac= ' num2str(Fac) ') ']} ...
            ,'color',condicolors(1,:));
set(v3,'Rotation',90);
sg=annotation('textbox',[0.3 0.01 0.4 0.05],'string',...
    {['PLS model: sum-EEG(-500ms) -> GC ^{* PLOT 21}'],char(datetime('now'))});
set([h0 h1 h2 h3 h4 v0 v1 v2 v3], 'fitboxtotext','on',...
    'edgecolor','none')

%% PLOT 22 Intervals histogram (p test with all intervals) (published)
% histogram of intertap intervals 
% organize intervals from all sessions
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);
cd /ssd/zhibin/1overf/
intervals_L={}; % 12 session x 12 trials
intervals_R={}; % 12 session x 12 trials
tic
for s=1:numSes
    runid=num2str(seeds(s,:));
    clear intervals
    load(['/ssd/zhibin/1overf/' runid '_2P/Cleaned_data/clean_' runid '.mat'],'intervals');
    for b=1:12
        intervals_L{s,b}=intervals{b}(:,1); % L subject 
        intervals_R{s,b}=intervals{b}(:,2); % R subject 
    end
end
toc % 76 sec

intervals_L;intervals_R;
% reorganzie cell sequence for intervals_all_L and intervals_all_R;
intervals_all_L=reshape(intervals_L',[],1);
intervals_all_R=reshape(intervals_R',[],1);
intervals_all_LR_all=[intervals_all_L; intervals_all_R]; % combine cell array
% organize intervals for 4 states in 2 syn types
syn2names={'Synch','Synco'};
intervals_2t_4ss={}; % 2 x 4 cells
tic
for t=1:2
    for ss=1:4
    cat_stateIntervals=[];
    cat_stateIntervals=cat(1,intervals_all_LR_all{Inds4_LR([synind(t,:) 36+synind(t,:)],ss)});
    intervals_2t_4ss{t,ss}=cat_stateIntervals;
    end
end
toc %

% transpose to (4 sessions x 2 groups)
intervals_2t_4ss=intervals_2t_4ss';
% compute mean and std (4 sessions x 2 groups)
intervals_2t_4ss_mean=[];intervals_2t_4ss_std=[];
intervals_2t_4ss_ste=[]; % standard error
for t=1:2
    for ss=1:4
        intervals_2t_4ss_mean(ss,t)=mean(intervals_2t_4ss{ss,t}./2); % and convert to ms
        intervals_2t_4ss_std(ss,t)=std(intervals_2t_4ss{ss,t}./2);
        intervals_2t_4ss_ste(ss,t)=std(intervals_2t_4ss{ss,t}./2)/sqrt(length(intervals_2t_4ss{ss,t}));
    end
end

% barplot with errorbar (published)
canvas(0.23, 0.4);
model_series = intervals_2t_4ss_mean;
model_error = intervals_2t_4ss_std;
b = bar(model_series, 'grouped');
b(1).FaceColor=darkgreen;b(2).FaceColor=pink;
hold on;
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(model_series);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars 
errorbar(x',model_series,model_error,'k','linestyle','none');
hold off
xticks(1:4);xticklabels({'Independent','Leader','Follower','Bidirectional'});
xl = get(gca,'XTickLabel');  
set(gca,'XTickLabel',xl,'fontsize',17,'FontWeight','bold')
xlim([0.5 4.5]);ylim([400 1000]);
ylabel('Mean tapping interval (ms)','FontSize',15); 
yl = get(gca,'YTickLabel');  
set(gca,'YTickLabel',yl,'fontsize',17,'FontWeight','bold')
set(gcf,'color','w'); % set background white for copying in ubuntu
% delete(findall(gcf,'type','annotation'))
% sg=annotation('textbox',[0.05 0.01 0.5 0.09],'string',...
%     {['mean tapping intervals(matched int) ^{ *PLOT 22}' char(datetime('now'))]})
% sg.Rotation=90
lg=legend({'Synch','Synco'},'location','north');lg.FontSize=17;

% prepend different colors for xticks
% https://www.mathworks.com/matlabcentral/answers/28537-set-the-yticklabel-to-different-colors
%{
        plot([1:1000]/2,stateERP(:,selectChan),'color',condicolors(ss,:));
        hold on;

        title(['Chan' labels{selectChan}]); xlabel('time (ms)');
        xlim([0 500]);ylim([-0.2 0.25])


lg=legend(states4names,'location','eastoutside');
lg.Position = [0.9475 0.4 0.01 0.25];
set(gcf,'color','w'); % set background white for copying in ubuntu
sgtitle('ERP of 4 states: Synch (top row); Synco (bottom row)');
%}

% statistical test (ANOVA) with unequal size
% prepare a matrictwo vectors for unbalance ANOVA, one for values, one for group labels
invertal_2t_4ss_anova=[];
invertal_2t_4ss_group=[];
for t=1:2
    for ss=1:4
        invertal_2t_4ss_anova=[invertal_2t_4ss_anova;intervals_2t_4ss{ss,t}./2]; % and convert to ms
        invertal_2t_4ss_group=[invertal_2t_4ss_group; ones(size(intervals_2t_4ss{ss,t})).*(t*10+ss)];
    end
end
[p,tbl,stats] = anova1(invertal_2t_4ss_anova,invertal_2t_4ss_group);
% statistical test (Two sample t test)
p={};
for t=1:2
    for ss=1:4
        % [h,p,ci,stats] = ttest2([intervals_2t_4ss{4,1}./2],[intervals_2t_4ss{1,1}./2]);
        [~,p{t,ss}(1),~,~] = ttest2([intervals_2t_4ss{ss,t}./2],[intervals_2t_4ss{1,1}./2]);
        [~,p{t,ss}(2),~,~] = ttest2([intervals_2t_4ss{ss,t}./2],[intervals_2t_4ss{2,1}./2]);
        [~,p{t,ss}(3),~,~] = ttest2([intervals_2t_4ss{ss,t}./2],[intervals_2t_4ss{3,1}./2]);
        [~,p{t,ss}(4),~,~] = ttest2([intervals_2t_4ss{ss,t}./2],[intervals_2t_4ss{4,1}./2]);
        [~,p{t,ss}(5),~,~] = ttest2([intervals_2t_4ss{ss,t}./2],[intervals_2t_4ss{1,2}./2]);
        [~,p{t,ss}(6),~,~] = ttest2([intervals_2t_4ss{ss,t}./2],[intervals_2t_4ss{2,2}./2]);
        [~,p{t,ss}(7),~,~] = ttest2([intervals_2t_4ss{ss,t}./2],[intervals_2t_4ss{3,2}./2]);
        [~,p{t,ss}(8),~,~] = ttest2([intervals_2t_4ss{ss,t}./2],[intervals_2t_4ss{4,2}./2]);
    end
end
figure;
for t=1:2
    for ss=1:4
        subplot(2,4,(t-1)*4+ss)
        bar(p{t,ss})
        ylim([0,0.5])
    end
end

%{ 
[p,tbl,stats] = anova1([...
    [H_L_uncouple_synch H_R_uncouple_synch]' ...
    [H_L_Llead_synch H_R_Rlead_synch]' ... 
    [H_R_Llead_synch H_L_Rlead_synch]' ...
    [H_L_mutual_synch H_R_mutual_synch]' ...
    [H_L_uncouple_synco H_R_uncouple_synco]' ...
    [H_L_Llead_synco H_R_Rlead_synco]' ... 
    [H_R_Llead_synco H_L_Rlead_synco]' ...
    [H_L_mutual_synco H_R_mutual_synco]'...
    ]); % one-way ANOVA (work)
% Bonferroni correction
[results,means,~,gnames] = multcompare(stats,"CriticalValueType","bonferroni"); % multi-comparison
tbl = array2table([results],"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])
tb2 = array2table([means],"VariableNames", ["Mean","Standard Error"])
p_bonf=results(:,6);
% statistical test (Two sample t test)
[h,p(1),ci,stats] = ttest2([H_L_uncouple_synch H_R_uncouple_synch]',[H_L_uncouple_synco H_R_uncouple_synco]') % Uncouple
[h,p(2),ci,stats] = ttest2([H_L_Llead_synch H_R_Rlead_synch]',[H_L_Llead_synco H_R_Rlead_synco]') % Leading
[h,p(3),ci,stats] = ttest2([H_R_Llead_synch H_L_Rlead_synch]',[H_R_Llead_synco H_L_Rlead_synco]','Alpha',0.01) % Following
[h,p(4),ci,stats] = ttest2([H_L_mutual_synch H_R_mutual_synch]',[H_L_mutual_synco H_R_mutual_synco]') % Mutual
sprintf('P values for each group: %0.3f\n',p)
[h,p,ci,stats] = ttest2([H_L_mutual_synch H_R_mutual_synch H_L_Llead_synch H_R_Rlead_synch ...
    H_R_Llead_synch H_L_Rlead_synch H_L_mutual_synch H_R_mutual_synch]',...
    [H_L_mutual_synco H_R_mutual_synco H_L_Llead_synco H_R_Rlead_synco ...
    H_R_Llead_synco H_L_Rlead_synco H_L_mutual_synco H_R_mutual_synco]') % Synch vs Synco
%}

%% PLOT 22-1 Intervals histogram (p test with mean interval in each block)
intervals_all_LR_all; % combined cell array from PLOT 22
for i=1:288
    intervals_all_LR_all_m(i)=mean(intervals_all_LR_all{i});
end
% organize intervals for 4 states in 2 syn types
syn2names={'Synch','Synco'};
% compute mean and std (4 sessions x 2 groups)
intervals_2t_4ss_mean=[];intervals_2t_4ss_std=[];
intervals_2t_4ss_ste=[]; % standard error
for t=1:2
    for ss=1:4
    intervals_2t_4ss_mean(ss,t)=mean(intervals_all_LR_all_m(Inds4_LR([synind(t,:) 36+synind(t,:)],ss)))./2;
    intervals_2t_4ss_std(ss,t)=std(intervals_all_LR_all_m(Inds4_LR([synind(t,:) 36+synind(t,:)],ss)))./2;
    intervals_2t_4ss_ste(ss,t)=(std(intervals_all_LR_all_m(Inds4_LR([synind(t,:) 36+synind(t,:)],ss)))./2)/...
        sqrt(length(intervals_all_LR_all_m(Inds4_LR([synind(t,:) 36+synind(t,:)],ss))));
    end
end
% barplot with errorbar
figure; % canvas(0.23, 0.4);
model_series = intervals_2t_4ss_mean;
model_error = intervals_2t_4ss_std;
b = bar(model_series, 'grouped');
b(1).FaceColor=darkgreen;b(2).FaceColor=pink;
hold on;
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(model_series);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars 
errorbar(x',model_series,model_error,'k','linestyle','none');
hold off
xticks(1:4);xticklabels({'Uncoupled','Leading','Following','Mutual'});
xl = get(gca,'XTickLabel');  
set(gca,'XTickLabel',xl,'fontsize',17,'FontWeight','bold')
xlim([0.5 4.5]);ylim([400 1000]);
ylabel('Mean tapping interval (ms)','FontSize',15); 
yl = get(gca,'YTickLabel');  
set(gca,'YTickLabel',yl,'fontsize',17,'FontWeight','bold')
set(gcf,'color','w'); % set background white for copying in ubuntu
% delete(findall(gcf,'type','annotation'))
% sg=annotation('textbox',[0.05 0.01 0.5 0.09],'string',...
%     {['mean tapping intervals(matched int) ^{ *PLOT 22}' char(datetime('now'))]})
% sg.Rotation=90
lg=legend({'Synch','Synco'},'location','north');lg.FontSize=17;


% statistical test (ANOVA) with unequal size
% prepare a matrictwo vectors for unbalance ANOVA, one for values, one for group labels
invertal_2t_4ss_anova=[];
invertal_2t_4ss_group=[];
for t=1:2
    for ss=1:4
        invertal_2t_4ss_anova=[invertal_2t_4ss_anova; (intervals_all_LR_all_m(Inds4_LR([synind(t,:) 36+synind(t,:)],ss))./2)']; % and convert to ms
        invertal_2t_4ss_group=[invertal_2t_4ss_group; (ones(size(intervals_all_LR_all_m(Inds4_LR([synind(t,:) 36+synind(t,:)],ss)))).*(t*10+ss))'];
    end
end
[p,tbl,stats] = anova1(invertal_2t_4ss_anova,invertal_2t_4ss_group);
[p,tbl,stats] = anova2(invertal_2t_4ss_anova,invertal_2t_4ss_group);


% statistical test (Two sample t test)
p={};
for t=1:2
    for ss=1:4
        % [h,p,ci,stats] = ttest2([intervals_2t_4ss{4,1}./2],[intervals_2t_4ss{1,1}./2]);
        [~,p{t,ss}(1),~,~] = ttest2([intervals_2t_4ss{ss,t}./2],[intervals_2t_4ss{1,1}./2]);
        [~,p{t,ss}(2),~,~] = ttest2([intervals_2t_4ss{ss,t}./2],[intervals_2t_4ss{2,1}./2]);
        [~,p{t,ss}(3),~,~] = ttest2([intervals_2t_4ss{ss,t}./2],[intervals_2t_4ss{3,1}./2]);
        [~,p{t,ss}(4),~,~] = ttest2([intervals_2t_4ss{ss,t}./2],[intervals_2t_4ss{4,1}./2]);
        [~,p{t,ss}(5),~,~] = ttest2([intervals_2t_4ss{ss,t}./2],[intervals_2t_4ss{1,2}./2]);
        [~,p{t,ss}(6),~,~] = ttest2([intervals_2t_4ss{ss,t}./2],[intervals_2t_4ss{2,2}./2]);
        [~,p{t,ss}(7),~,~] = ttest2([intervals_2t_4ss{ss,t}./2],[intervals_2t_4ss{3,2}./2]);
        [~,p{t,ss}(8),~,~] = ttest2([intervals_2t_4ss{ss,t}./2],[intervals_2t_4ss{4,2}./2]);
    end
end
%{ 
% Bonferroni correction
[results,means,~,gnames] = multcompare(stats,"CriticalValueType","bonferroni"); % multi-comparison
tbl = array2table([results],"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])
tb2 = array2table([means],"VariableNames", ["Mean","Standard Error"])
p_bonf=results(:,6);
%}
%% Working from home
cd D:\360MoveData\Users\alienware\Documents\GitHub\Motor_cordination\1_over_f\data_analysis\2P_testing

% addpath
addpath D:\360MoveData\Users\alienware\Documents\GitHub\Motor_cordination\1_over_f\data_analysis\2P_testing

% load workspace
cd G:\1overf\all_session20220713_1005
load('20220713_1005workspace.mat')
load('FB_all_samples_all.mat')

% git commit
cd D:\360MoveData\Users\alienware\Documents\GitHub\Motor_cordination
run git_control.m   

%% saves all variables from the current workspace
tic
save([num2str(seed) 'workspace.mat']);
toc
% takes about 11 min for one session

cd /ssd/zhibin/1overf/all_session20220713_1005

% save
tic
save(['20220713_1005workspace.mat']);
toc

% load
tic 
load(['20220713_1005workspace.mat']);
toc
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
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);
sessions={'synch','synco','synch','synco','synch','synco','synch','synco','synch','synco',...
    'synch','synco'};
cd /ssd/zhibin/1overf/
for s=1:numSes
    mkdir(['/ssd/zhibin/1overf/' num2str(seeds(s,:)) '_2P/Cleaned_data/Plots'])
    cd(['/ssd/zhibin/1overf/' num2str(seeds(s,:)) '_2P/Cleaned_data'])
    movefile cleanup*.png Plots
    cd /ssd/zhibin/1overf/
end

%%
