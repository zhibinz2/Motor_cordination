clear
%% addpath
addpath /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing
addpath /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/DFA

% % d estimate and removal
% % http://www.lucafaes.net/LMSE-MSE_VARFI.html
% addpath /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/MSE-VARFI
% % Granger Causality
% addpath(genpath('/home/zhibinz2/Documents/GitHub/MVGC1'));
% cd /home/zhibinz2/Documents/GitHub/MVGC1
% run startup
% addpath /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/VAR_Granger


%% Color Scheme
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

% % test color
% showcolor=black;
% imagesc(cat(3,showcolor(1),showcolor(2),showcolor(3)));

% colormap
cd /home/zhibinz2/Documents/GitHub/matlab-archive/hnlcode/common/gen_code/color
addpath /home/zhibinz2/Documents/GitHub/matlab-archive/hnlcode/common/gen_code/color
hnc = hotncold(100);
hottermap=hotter(100);
% 
% % test colormap
% mycolormap=hottermap;
% mymat=repmat(-0.99:0.02:0.99,100,1); % transitional color
% % mymat=(rand(100)-0.5)*2; % random from -1 to 1
% imagesc(mymat);colorbar;
% colormap(mycolormap);

%% Xcorr10Lag
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);

% (Before d removal)********
% d estimation from H (DFA method) based on the good matched intervals
Xcorr10Lag=nan(numSes,12,21);
cd /ssd/zhibin/1overf/
for r=1:numSes;
    clear intervals conditions sortorders
    runid = num2str(seeds(r,:));
    path = [runid '_2P/Cleaned_data/'];
    load([path  'clean_' runid '.mat'],'conditions','intervals');
    % sort order and plot results
    [x,sortorder]=sort(conditions);
    for j = 1:12
        clear intL_good_dmean intR_good_dmean
        % remove the mean
        intL_good_dmean=intervals{sortorder(j)}(:,1)-mean(intervals{sortorder(j)}(:,1));
        intR_good_dmean=intervals{sortorder(j)}(:,2)-mean(intervals{sortorder(j)}(:,2));
        % Xcorr based on int_dmean_drm (before d removal)********
        r12=[];lags12=[];
        [r12,lags12]=xcorr(intL_good_dmean,intR_good_dmean,10,'normalized');
        Xcorr10Lag(r,j,1:21)=r12;
    end
end

%%
% % indicies in the synch and synco time sequence
% synchind=[1:3 7:9 13:15 19:21 25:27 31:33]; % 3 trials x 6 sessions
% syncoind=[4:6 10:12 16:18 22:24 28:30 34:36]; % 3 trials x 6 sessions
% synind=[synchind;syncoind];
% 
% 
% % orgainized the synchronization trials into 4 conditions
% condi4Ind={[1:3],[4:6],[7:9],[10:12]};
% condi4names={'Uncouple','L-lead','R-lead','Mutual'};
syn2Ind={[1:2:11],[2:2:12]};syn2names={'Synch','Synco'};

%%
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
