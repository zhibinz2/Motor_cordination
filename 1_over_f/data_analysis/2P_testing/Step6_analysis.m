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

%     figureName=['Corr_DFA-trial-' num2str(i)];
%     saveas(gcf,figureName,'jpg');

end
close all;

%% PLOT 2: Auto&Xcorr and DFA and Autocov Granger Causality and PSA for each trial
cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/ARFIMA/Sim_plots
cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/DFA

cd /ssd/zhibin/1overf/20220713_2P/Segmented_data
cd /ssd/zhibin/1overf/20220721_2P/Segmented_data

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
maxfreq=50; win=10; % 10 Second with df of 0.1Hz

% Calculation on BP
y1=[];y2=[];y12=[];XcorrR=[];
nTrials=12;
alpha = 0.05; 
for i=1:nTrials
    % extract BP01
    y1(1,i).BP01=(BP(1).BP{i})'; % in boolean
    y2(1,i).BP01=(BP(2).BP{i})'; % in boolean 
   
    % extract BP intervals
    y1(1,i).BPint=Calinterval(y1(1,i).BP01)./sr; % in second
    y2(1,i).BPint=Calinterval(y2(i).BP01)./sr; % in second
    
    MinLength=min([length(y1(i).BPint) length(y2(i).BPint)]);

    % remove the mean
    y1(1,i).BPint = y1(i).BPint-mean(y1(i).BPint);
    y2(1,i).BPint = y2(i).BPint-mean(y2(i).BPint);
    
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

% Calculate EEG Power of the whole trial
for i=1:nTrials
    y1(1,i).EEG=zscore(EEG(1).EEG{i},[],1);
    y2(1,i).EEG=zscore(EEG(2).EEG{i},[],1);
    
    % create a function that reorganize EEG into time x chan x chunks
    % feed it to allspectra (for 2 s chunks) or spectra3 (for 10s chunks)
    [~,freqs, ~, y1(i).pow, ~] = spectra3(y1(i).EEG,sr,maxfreq,win);
    [~,freqs, ~, y2(i).pow, ~] = spectra3(y2(i).EEG,sr,maxfreq,win);
    % later examine pow eppow corr cprod in imagesc and create myallspectra
end
    
%% PLOT 2: plotting
cd /ssd/zhibin/1overf/20220713_2P/Segmented_data/Plots/Corr_DFA_GC_PSA
cd /ssd/zhibin/1overf/20220721_2P/Segmented_data/Plots/Corr_DFA_GC_PSA

% check condiSeq
condiSeq=allPerm;

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
