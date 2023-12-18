

seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);

for n=11:12
    ses_dir=['/ssd/zhibin/1overf/' num2str(seeds(n,:)) '_2P/Segmented_data/'];
    cd(ses_dir)
    clear EMG
    load(['data_variables' num2str(seeds(n,:)) '.mat'], 'EMG');
    EMG_filtered_L=EMG(1).EMG;
    EMG_filtered_R=EMG(2).EMG;
    cd /ssd/zhibin/1overf/EMG
    save(['EMG_' num2str(seeds(n,:)) '.mat'], 'EMG_filtered_L','EMG_filtered_R');
end


if n==3 5 7 10
    EMG_filtered_L=EMG(1).filtered_EMG;
    EMG_filtered_R=EMG(2).filtered_EMG;
end