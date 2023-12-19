

seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);

for n=[1 2 4 6 8 9 11 12]
    ses_dir=['/ssd/zhibin/1overf/' num2str(seeds(n,:)) '_2P/Segmented_data/'];
    cd(ses_dir)
    clear EMG
    load(['data_variables' num2str(seeds(n,:)) '.mat'], 'EMG');
    EMG_filtered_L=EMG(1).EMG;
    EMG_filtered_R=EMG(2).EMG;
    cd /ssd/zhibin/1overf/EMG
    save(['EMG_' num2str(seeds(n,:)) '.mat'], 'EMG_filtered_L','EMG_filtered_R');
end


% n==3 no segment 
% 5 no segment 
% 7 no segment 
% 10 load failed
% EMG_filtered_L=EMG(1).filtered_EMG;
% EMG_filtered_R=EMG(2).filtered_EMG;

n=3
ses_dir=['/ssd/zhibin/1overf/' num2str(seeds(n,:)) '_2P/Segmented_data/'];
cd(ses_dir)
% load(['data_variables' num2str(seeds(n,:)) '.mat'], 'EMG');
load('data_structure20220804.mat', 'data') % CHECK IF NOT ENVELOP
EMG_filtered_L=data(1).EMG;
EMG_filtered_R=data(2).EMG;
cd /ssd/zhibin/1overf/EMG
save(['EMG_' num2str(seeds(n,:)) '.mat'], 'EMG_filtered_L','EMG_filtered_R');

n=5
ses_dir=['/ssd/zhibin/1overf/' num2str(seeds(n,:)) '_2P/Segmented_data/'];
cd(ses_dir)
% load(['data_variables' num2str(seeds(n,:)) '.mat'], 'EMG'); 
load('data_structure20220810.mat') % Failed
load('data_variables20220810.mat')
EMG_filtered_L=data(1).EMG;
EMG_filtered_R=data(2).EMG;
cd /ssd/zhibin/1overf/EMG
save(['EMG_' num2str(seeds(n,:)) '.mat'], 'EMG_filtered_L','EMG_filtered_R');

n=7
ses_dir=['/ssd/zhibin/1overf/' num2str(seeds(n,:)) '_2P/Segmented_data/'];
cd(ses_dir)
load(['data_structure' num2str(seeds(n,:)) '.mat'], 'data'); % CHECK IF NOT ENVELOP
EMG_filtered_L=data(1).EMG;
EMG_filtered_R=data(2).EMG;
cd /ssd/zhibin/1overf/EMG
save(['EMG_' num2str(seeds(n,:)) '.mat'], 'EMG_filtered_L','EMG_filtered_R');

n=10
ses_dir=['/ssd/zhibin/1overf/' num2str(seeds(n,:)) '_2P/Segmented_data/'];
load('data_variables2022100401.mat') % Failed to Read
EMG_filtered_L=data(1).EMG;
EMG_filtered_R=data(2).EMG;
cd /ssd/zhibin/1overf/EMG
save(['EMG_' num2str(seeds(n,:)) '.mat'], 'EMG_filtered_L','EMG_filtered_R');