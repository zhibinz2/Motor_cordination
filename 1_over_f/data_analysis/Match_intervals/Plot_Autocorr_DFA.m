%% color
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
% showcolor=mediumblue;
% imagesc(cat(3,showcolor(1),showcolor(2),showcolor(3)));
%% Add my repo
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/DFA
%% sort order and plot results
seed=20220713;
seed=20220721;
seed=20220804;
seed=20220808;
seed=20220810;
seed=20220811;
seed=20220815;
seed=20220816;


[x,sortorder]=sort(conditions);
canvas(0.5,0.3);
for j = 1:12
    subplot(4,3,j)
    plot(intervals{sortorder(j)})
    title(conditions(sortorder(j)))
end
sgtitle([num2str(seed) ' ' sessionTypes{session} ' matched intervals'])
eval(['print ' outpath 'cleanup_all_12.png -dpng'])

%% Autocorr before and after
canvas(1,0.5);
for j=1:12
    % Before aligment
    intL{sortorder(j)};
    intR{sortorder(j)};
    % remove the mean
    intL_dmean=intL{sortorder(j)}-mean(intL{sortorder(j)});
    intR_dmean=intR{sortorder(j)}-mean(intR{sortorder(j)});
    
    subplot(4,12,j)
    acf=[];lags=[];bounds=[];
    [acf,lags,bounds] = autocorr(intL_dmean,20);
    stem(lags,acf,'r'); xlabel('Lag');ylabel('\rho(k)');
    ylim([-0.3 1]);
    hold on;
    h = line(lags,[bounds'].*ones(size(lags)));
    set(h,'color',[1 0 1]);
    title({[conditionNames{condiSeq(sortorder(j))}],['Before matching']},'Color',condicolors(condiSeq(sortorder(j)),:));
    subtitle(['intL (n=' num2str(length(intL_dmean)) '): Autocorr'],'Color',condicolors(condiSeq(sortorder(j)),:));
    
    subplot(4,12,24+j)
   
    stem(lags,acf,'b'); xlabel('Lag');ylabel('\rho(k)');
    ylim([-0.3 1]);% axis([0 MinLength/2 -0.3 1]) 
    hold on;
    h = line(lags,[bounds'].*ones(size(lags)));
    set(h,'color',[1 0 1]);
    title('Before matching','Color',condicolors(condiSeq(sortorder(j)),:));
    subtitle(['intR (n=' num2str(length(intR_dmean)) '): Autocorr'],'Color',condicolors(condiSeq(sortorder(j)),:));
    
    % After aligment
    intervals{sortorder(j)}(:,1);
    intervals{sortorder(j)}(:,2);
    % remove the mean
    intL_good_dmean=intervals{sortorder(j)}(:,1)-mean(intervals{sortorder(j)}(:,1));
    intR_good_dmean=intervals{sortorder(j)}(:,2)-mean(intervals{sortorder(j)}(:,2));
    
    % Autocorr
    subplot(4,12,12+j)
    acf=[];lags=[];bounds=[];
    [acf,lags,bounds] = autocorr(intL_good_dmean,20);
    stem(lags,acf,'r'); xlabel('Lag');ylabel('\rho(k)');
    ylim([-0.3 1]);
    hold on;
    h = line(lags,[bounds'].*ones(size(lags)));
    set(h,'color',[1 0 1]);
    title('After matching','Color',condicolors(condiSeq(sortorder(j)),:));
    subtitle(['intL (n=' num2str(length(intL_good_dmean)) '): Autocorr'],'Color',condicolors(condiSeq(sortorder(j)),:));
    
    
    subplot(4,12,36+j)
    acf=[];lags=[];bounds=[];
    [acf,lags,bounds] = autocorr(intR_good_dmean,20);
    stem(lags,acf,'b'); xlabel('Lag');ylabel('\rho(k)');
    ylim([-0.3 1]);% axis([0 MinLength/2 -0.3 1]) 
    hold on;
    h = line(lags,[bounds'].*ones(size(lags)));
    set(h,'color',[1 0 1]);
    title('After matching','Color',condicolors(condiSeq(sortorder(j)),:));
    subtitle(['intR (n=' num2str(length(intR_good_dmean)) '): Autocorr'],'Color',condicolors(condiSeq(sortorder(j)),:));
end
sgtitle([num2str(seed) ' ' sessionTypes{session}]);
eval(['print ' outpath 'Autocorr_Compare_all_12.png -dpng'])

%% DFA before and after
canvas(1,0.3);
for j=1:12
    % Before aligment
    intL{sortorder(j)};
    intR{sortorder(j)};
    % remove the mean
    intL_dmean=intL{sortorder(j)}-mean(intL{sortorder(j)});
    intR_dmean=intR{sortorder(j)}-mean(intR{sortorder(j)});
    % After aligment
    intervals{sortorder(j)}(:,1);
    intervals{sortorder(j)}(:,2);
    % remove the mean
    intL_good_dmean=intervals{sortorder(j)}(:,1)-mean(intervals{sortorder(j)}(:,1));
    intR_good_dmean=intervals{sortorder(j)}(:,2)-mean(intervals{sortorder(j)}(:,2));

    subplot(2,12,j)
    y1=[];y2=[];
    n=[];F_n=[];FitValues=[];
    [~,Alpha,n,F_n,FitValues]=DFA_main(intL_dmean);
    y1.Alpha=Alpha;
    y1.n=n;
    y1.F_n=F_n;
    y1.FitValues=FitValues;
    n=[];F_n=[];FitValues=[];
    [~,Alpha,n,F_n,FitValues]=DFA_main(intL_good_dmean);
    y2.Alpha=Alpha;
    y2.n=n;
    y2.F_n=F_n;
    y2.FitValues=FitValues;  

    plot(log10(y1.n),log10(y1.F_n),'rx');
    hold on; plot(log10(y1.n),y1.FitValues,'r');
    plot(log10(y2.n),log10(y2.F_n),'kx');
    hold on; plot(log10(y2.n),y2.FitValues,'k');
    % axis([0.6 2 -1.8 0]);
    % legend({'L',['H=' num2str(y1.Alpha)],'R',['H=' num2str(y2.Alpha)],},'Location','southeast','NumColumns',2);
    xlabel('Scale [log10(n)]') % win_lengths 
    ylabel('RMS [log10(F-n)]') % RMS values

    title({[conditionNames{condiSeq(sortorder(j))}],['intL (n=' num2str(length(intL_dmean)) '->' num2str(length(intL_good_dmean)) ')']},'Color',condicolors(condiSeq(sortorder(j)),:));
    subtitle({['H-before=' num2str(round(y1.Alpha,2))],['H-after=' num2str(round(y2.Alpha,2)) ' (black)']},'Color',condicolors(condiSeq(sortorder(j)),:))
    ylim([1.5 3.4]);


    subplot(2,12,12+j)
    y1=[];y2=[];
    n=[];F_n=[];FitValues=[];
    [~,Alpha,n,F_n,FitValues]=DFA_main(intR_dmean);
    y1.Alpha=Alpha;
    y1.n=n;
    y1.F_n=F_n;
    y1.FitValues=FitValues;
    n=[];F_n=[];FitValues=[];
    [~,Alpha,n,F_n,FitValues]=DFA_main(intR_good_dmean);
    y2.Alpha=Alpha;
    y2.n=n;
    y2.F_n=F_n;
    y2.FitValues=FitValues;   

    plot(log10(y1.n),log10(y1.F_n),'bx');
    hold on; plot(log10(y1.n),y1.FitValues,'b');
    plot(log10(y2.n),log10(y2.F_n),'kx');
    hold on; plot(log10(y2.n),y2.FitValues,'k');
    % axis([0.6 2 -1.8 0]);
    % legend({'L',['H=' num2str(y1.Alpha)],'R',['H=' num2str(y2.Alpha)],},'Location','southeast','NumColumns',2);
    xlabel('Scale [log10(n)]') % win_lengths 
    ylabel('RMS [log10(F-n)]') % RMS values
    title(['intR (n=' num2str(length(intR_dmean)) '->' num2str(length(intR_good_dmean)) ')'],'Color',condicolors(condiSeq(sortorder(j)),:));
    subtitle({['H-before=' num2str(round(y1.Alpha,2))], ['H-after=' num2str(round(y2.Alpha,2)) ' (black)']},'Color',condicolors(condiSeq(sortorder(j)),:));
    ylim([1.5 3.4]);

end
sgtitle([num2str(seed) ' ' sessionTypes{session} '  DFA']);
eval(['print ' outpath 'DFA_Compare_all_12.png -dpng'])

%% Go to the directory to check
cd([outpath])
close all;
