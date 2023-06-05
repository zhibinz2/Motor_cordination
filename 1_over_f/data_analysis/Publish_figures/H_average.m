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


%%
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);

H=zeros(2,12,numSes);
cd /ssd/zhibin/1overf/

for r=1:numSes
    clear intervals conditions sortorders
    runid = num2str(seeds(r,:));
    % load unmatched int.
    path = [runid '_2P/Cleaned_data/'];
    load([path  'clean_' runid '.mat'],'intervals','conditions');
    % sort order and plot results
    % [x,sortorder]=sort(condiSeq);
    [x,sortorder]=sort(conditions);
    for j = 1:12
        % For unmatched int.
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
 

%% average all trials of the same condition in synch/o (2 subplots)
H_L_uncouple_synch = reshape(H(1,[1:3],[1:2:11]),1,[])  % L,uncouple,synch
H_L_Llead_synch    = reshape(H(1,[4:6],[1:2:11]),1,[])  % L,L lead,synch
H_L_Rlead_synch    = reshape(H(1,[7:9],[1:2:11]),1,[])  % L,R lead,synch
H_L_mutual_synch   = reshape(H(1,[10:12],[1:2:11]),1,[])% L, mutual,synch
H_R_uncouple_synch = reshape(H(2,[1:3],[1:2:11]),1,[])  % R,uncouple,synch
H_R_Llead_synch    = reshape(H(2,[4:6],[1:2:11]),1,[])  % R,L lead,synch
H_R_Rlead_synch    = reshape(H(2,[7:9],[1:2:11]),1,[])  % R,R lead,synch
H_R_mutual_synch   = reshape(H(2,[10:12],[1:2:11]),1,[])% R, mutual,synch
H_LR_synch_mean = [mean(H_L_uncouple_synch) mean(H_R_uncouple_synch);...
    mean(H_L_Llead_synch) mean(H_R_Llead_synch);...
    mean(H_L_Rlead_synch) mean(H_R_Rlead_synch);...
    mean(H_L_mutual_synch) mean(H_R_mutual_synch)];

H_L_uncouple_synco = reshape(H(1,[1:3],[2:2:12]),1,[])  % R,uncouple,synco
H_L_Llead_synco    = reshape(H(1,[4:6],[2:2:12]),1,[])  % R,L lead,synco
H_L_Rlead_synco    = reshape(H(1,[7:9],[2:2:12]),1,[])  % R,R lead,synco
H_L_mutual_synco   = reshape(H(1,[10:12],[2:2:12]),1,[])% R, mutual,synco
H_R_uncouple_synco = reshape(H(2,[1:3],[2:2:12]),1,[])  % R,uncouple,synco
H_R_Llead_synco    = reshape(H(2,[4:6],[2:2:12]),1,[])  % R,L lead,synco
H_R_Rlead_synco    = reshape(H(2,[7:9],[2:2:12]),1,[]) % R,R lead,synco
H_R_mutual_synco   = reshape(H(2,[10:12],[2:2:12]),1,[])% R, mutual,synco
H_LR_synco_mean = [mean(H_L_uncouple_synco) mean(H_R_uncouple_synco);...
    mean(H_L_Llead_synco) mean(H_R_Llead_synco);...
    mean(H_L_Rlead_synco) mean(H_R_Rlead_synco);...
    mean(H_L_mutual_synco) mean(H_R_mutual_synco)];

%% separate synch and synco with error bar
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

%% grouped barplot with errorbar (4 subplots) and statistical tests
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


