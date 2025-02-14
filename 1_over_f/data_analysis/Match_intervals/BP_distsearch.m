clear all
close all

%%
seed=20220713;
seed=20220721;
seed=20220804;
seed=20220808;
seed=20220810;
seed=20220811;
seed=20220815;
seed=20220816;
seed=20221003;
%%
clear
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
    2022100402;20221005];
sessions={'synch','synco','synch','synco','synch','synco','synch','synco','synch','synco',...
    'synch','synco'};
cd /ssd/zhibin/1overf/

%% make folder 
mkdir Cleaned_data

mkdir(['/ssd/zhibin/1overf/' num2str(seed) '_2P/Cleaned_data/']);

%% 
% addpath /ssd/zhibin/1overf/20220808_2P/Copy_of_Ramesh20220909
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/Match_intervals

%% run match intervals
tic
clear
for i = 1:12
    cd /ssd/zhibin/1overf
    varlist1 = who();
    seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
    tooshort = 700;
    runid = num2str(seeds(i,:));
    path = [runid '_2P/Segmented_data/'];
    outpath = [runid '_2P/Cleaned_data/'];
    load([path 'data_variables' runid '.mat'])
    

    bpchan = [32]; %for python
    conditions = condiSeq;
    channels = chaninfo;

    for j = 1:12
	    BPL = BP(1).BP{j};
	    BPR = BP(2).BP{j};
	    [samples{j}, intervals{j}, BPL, BPR, intL{j}, intR{j}] = BP_nearest(BPL,BPR, tooshort);
	    dataL{j} = [EEG(1).EEG{j} BPL];
	    dataR{j} = [EEG(2).EEG{j} BPR];
    end;
    eval(['save ' outpath 'clean_' runid '.mat dataL dataR conditions channels sr bpchan samples intervals tooshort'...
        ' seed session sessionTypes labels conditionNames'])
    display(['Completed session ' runid]);

    varlist2 = who();
    varlist_diff = setdiff(varlist2,varlist1);
    % clear all new added variables
    clear(varlist_diff{:})
    clear varlist2 varlist_diff
end
toc





