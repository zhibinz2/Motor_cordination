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

%% condition_all
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);
condition_all=[];
for s=1:numSes
    clear conditions
    runid=num2str(seeds(s,:));
    load(['/ssd/zhibin/1overf/' runid '_2P/Cleaned_data/clean_' runid '.mat'],'conditions');
    condition_all(s,:)=conditions; % 12 session x 12 trials
end
% reshape into a vector in time sequence
condition_all=reshape(condition_all',[],1); % 144 x 1 

%%

% PLOT 22 Intervals histogram (p test with all intervals) (published)
% histogram of intertap intervals 
% organize intervals from all sessions

cd /ssd/zhibin/1overf/
intervals_L={}; % 12 session x 12 trials
intervals_R={}; % 12 session x 12 trials
for s=1:numSes
    runid=num2str(seeds(s,:));
    clear intervals
    load(['/ssd/zhibin/1overf/' runid '_2P/Cleaned_data/clean_' runid '.mat'],'intervals');
    for b=1:12
        intervals_L{s,b}=intervals{b}(:,1); % L subject 
        intervals_R{s,b}=intervals{b}(:,2); % R subject 
    end
end

%%
% indicies in the synch and synco time sequence
synchind=[1:3 7:9 13:15 19:21 25:27 31:33]; % 3 trials x 6 sessions
syncoind=[4:6 10:12 16:18 22:24 28:30 34:36]; % 3 trials x 6 sessions
synind=[synchind;syncoind];


% indices for 4 states
states4names={'Independent','Leading','Following','Bidrectional'};
% find the indices for each condition in L R conbined sequence (4 states)
% and organize for PLS
uncoupleInd_LR=[find(condition_all==1);12*numSes+find(condition_all==1)];
leadingInd_LR=[find(condition_all==2);12*numSes+find(condition_all==3)];
followingInd_LR=[find(condition_all==3);12*numSes+find(condition_all==2)];
mutualInd_LR=[find(condition_all==4);12*numSes+find(condition_all==4)];
Inds4_LR=[uncoupleInd_LR leadingInd_LR followingInd_LR mutualInd_LR];

%% Intervals
intervals_L;intervals_R;
% reorganzie cell sequence for intervals_all_L and intervals_all_R;
intervals_all_L=reshape(intervals_L',[],1);
intervals_all_R=reshape(intervals_R',[],1);
intervals_all_LR_all=[intervals_all_L; intervals_all_R]; % combine cell array
% organize intervals for 4 states in 2 syn types
syn2names={'Synch','Synco'};
intervals_2t_4ss={}; % 2 x 4 cells
for t=1:2
    for ss=1:4
        cat_stateIntervals=[];
        cat_stateIntervals=cat(1,intervals_all_LR_all{Inds4_LR([synind(t,:) 36+synind(t,:)],ss)});
        intervals_2t_4ss{t,ss}=cat_stateIntervals;
    end
end


% transpose to (4 sessions x 2 groups)
intervals_2t_4ss=intervals_2t_4ss';
% compute mean and std (4 sessions x 2 groups)
intervals_2t_4ss_mean=[];
intervals_2t_4ss_std=[];
intervals_2t_4ss_ste=[]; % standard error
intervals_2t_4ss_25=[]; % percentile
intervals_2t_4ss_75=[]; % percentile
for t=1:2
    for ss=1:4
        intervals_2t_4ss_mean(ss,t)=mean(intervals_2t_4ss{ss,t}./2); % and convert to ms
        intervals_2t_4ss_std(ss,t)=std(intervals_2t_4ss{ss,t}./2);
        intervals_2t_4ss_ste(ss,t)=std(intervals_2t_4ss{ss,t}./2)/sqrt(length(intervals_2t_4ss{ss,t}));
        intervals_2t_4ss_25(ss,t)=prctile(intervals_2t_4ss{ss,t}./2,25);
        intervals_2t_4ss_75(ss,t)=prctile(intervals_2t_4ss{ss,t}./2,75);
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
errorbar(x',model_series,model_error,'k','linestyle','none','LineWidth',2);
% Mark the 25-75 percentile on eror bars
for t=1:2
    for ss=1:4
        plot(x(t,ss),intervals_2t_4ss_75(ss,t),'b_','LineWidth',2)
        plot(x(t,ss),intervals_2t_4ss_25(ss,t),'b_','LineWidth',2)
    end
end
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



%% Boxchart https://www.mathworks.com/help/matlab/ref/boxchart.html
syn2names={'Synch','Synco'};
st4names={'Independent','Leader','Follower','Bidirectional'};
% organize the arrays into a table
intervals=[];syn={};st={};
for t=1:2
    for ss=1:4
        intervals_temp=intervals_2t_4ss{ss,t}./2; % and convert to ms
        len_temp=length(intervals_temp);
        intervals=[intervals; intervals_temp];
        syn=[syn;repelem(syn2names(t),len_temp)'];
        st=[st;repelem(st4names(ss),len_temp)'];
    end
end

% put the above 3 array into 3 columns of a table
tdl=table(intervals);
tdl.syn=categorical(syn);
tdl.st=categorical(st,st4names);
boxchart(tdl.st,tdl.intervals,'GroupByColor',tdl.syn);
ylabel('intervals (ms)');
legend

%% boxplot https://www.mathworks.com/matlabcentral/answers/566745-how-to-group-boxplots-on-matlab


%%
cd /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/Publish_figures
