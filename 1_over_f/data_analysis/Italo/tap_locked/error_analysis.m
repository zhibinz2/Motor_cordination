cd /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/Italo/tap_locked
load('ses12_errors.mat')
ses12_errors=cell_data;
load('ses12_seeds_cell_array.mat');
ses12_seeds_cell_array=cell_data;
load('ses12_start_stop_cell_array.mat')
ses12_start_stop_cell_array=cell_data;
load('ses12_joint_seq_cell_array.mat');
ses12_joint_seq_cell_array=cell_data;

%% symbol and mean error reorganzied with the same sequense
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);
% Convert each number to a string and store it in a cell array
seeds_cell = arrayfun(@num2str, seeds, 'UniformOutput', false)

ses12_errors{13} % sequense same as seeds confirmed

ses12_error_reorg=cell(12,12);
for ses=1:12
    ses12_error_reorg{ses,1}=ses12_errors{ses}.tr1
    ses12_error_reorg{ses,2}=ses12_errors{ses}.tr2
    ses12_error_reorg{ses,3}=ses12_errors{ses}.tr3
    ses12_error_reorg{ses,4}=ses12_errors{ses}.tr4
    ses12_error_reorg{ses,5}=ses12_errors{ses}.tr5
    ses12_error_reorg{ses,6}=ses12_errors{ses}.tr6
    ses12_error_reorg{ses,7}=ses12_errors{ses}.tr7
    ses12_error_reorg{ses,8}=ses12_errors{ses}.tr8
    ses12_error_reorg{ses,9}=ses12_errors{ses}.tr9
    ses12_error_reorg{ses,10}=ses12_errors{ses}.tr10
    ses12_error_reorg{ses,11}=ses12_errors{ses}.tr11
    ses12_error_reorg{ses,12}=ses12_errors{ses}.tr12
end

error_reorg=cell(12,1);
for ses=1:12
    error_reorg{ses}=[ses12_error_reorg{ses,:}];
end

start_stop_reorg=cell(12,1);
joint_seq_reorg=cell(12,1);
for ses=1:12
    ind_reorg=find(strcmp(ses12_seeds_cell_array,seeds_cell{ses}));
    start_stop_reorg{ses}=ses12_start_stop_cell_array{ind_reorg};
    joint_seq_reorg{ses}=ses12_joint_seq_cell_array{ind_reorg};
end

UniqJtS=cell(12,1)
numUniqJtS=nan(12,1);
for ses=1:12
    UniqJtS{ses}=unique(joint_seq_reorg{ses})
    numUniqJtS(ses)=length(UniqJtS{ses});
end

meanError=cell(12,1); 
seError=cell(12,1);
for ses=1:12
    meanError{ses}=zeros(numUniqJtS(ses),1);
    seError{ses}=zeros(numUniqJtS(ses),1);
    for JtS=1:numUniqJtS(ses)
        Inds_js=find(joint_seq_reorg{ses}==UniqJtS{ses}(JtS));
        meanError{ses}(JtS)=mean(error_reorg{ses}(Inds_js));
        seError{ses}(JtS)=std(error_reorg{ses}(Inds_js))/sqrt(length(Inds_js));
    end
end

% standard offset
S_offset_synch=0;
S_offset_synco=2000/1.3/2/2;

figure
clf
for ses=1:12
    subplot(3,4,ses)
    % X values (optional)
    x = 1:length(meanError{ses});  % X-axis positions for each group
    
    % Plot mean with error bars (mean ± standard error)
    errorbar(x, meanError{ses}./2, seError{ses}./2, 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red');
    xlabel('Joint Symbols');
    ylabel('Mean Error (ms)');
    title(['Mean with Standard Error: ses' num2str(ses)]);
    grid on;
    ylim([-10 400])

    if mod(ses, 2) == 0 % even
        yline(S_offset_synco,'m-')
    else
        yline(S_offset_synch,'m-')
    end
end



%% load the conditions
