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

ses12_errors{13} % sequence same as seeds confirmed

% reorganize errors from each trials
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

% concatenate errors in 12 trials 
error_reorg=cell(12,1);
for ses=1:12
    error_reorg{ses}=[ses12_error_reorg{ses,:}];
end

% reorganize start-stop time and joint symbols in session sequence
start_stop_reorg=cell(12,1);
joint_seq_reorg=cell(12,1);
for ses=1:12
    ind_reorg=find(strcmp(ses12_seeds_cell_array,seeds_cell{ses}));
    start_stop_reorg{ses}=ses12_start_stop_cell_array{ind_reorg}+1;
    joint_seq_reorg{ses}=ses12_joint_seq_cell_array{ind_reorg};
end

% unique joint symbols and their counts
UniqJtS=cell(12,1)
numUniqJtS=nan(12,1);
for ses=1:12
    UniqJtS{ses}=unique(joint_seq_reorg{ses})
    numUniqJtS(ses)=length(UniqJtS{ses});
end

%% mean error and standard error of joint symbols
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
        subtitle('synco')
    else
        yline(S_offset_synch,'m-')
        subtitle('synch')
    end
end



%% load the conditions
data_path = '/home/zhibinz2/Documents/GitHub/finger_tapping_behavioral_data/';
data_path = 'C:\Users\zhouz\GitHub\finger_tapping_behavioral_data/'
% load conditions
condition_all=[];
for s=1:numSes
    clear conditions
    runid=num2str(seeds(s,:));
    load([data_path 'clean_' runid '_bpchan.mat'],'conditions');
    condition_all(s,:)=conditions; % 12 session x 12 trials
end

% use start_stop_reorg to assign the condition to each error in time sequeence
condition_reorg=cell(12,1);
for s=1:numSes
    condition_reorg{s}=nan(length(joint_seq_reorg{s}),1);
end
for s=1:numSes
    for tr=1:12
        tr_condi=condition_all(s,tr);
        condition_reorg{s}(start_stop_reorg{s}(tr,1):start_stop_reorg{s}(tr,2))=tr_condi;
    end
end
% replace condition 2 and 3 with 2, and 4 into 3.
for s=1:numSes
    condition_reorg{s}(find(condition_reorg{s}==3))=2;
    condition_reorg{s}(find(condition_reorg{s}==4))=3;
end
%%  Examine relationship btw error and condi
err3m12=cell(12,1);
err3se12=cell(12,1);
for ses=1:12; % pick a session
    subplot(3,4,ses)
    % compute mean and standard error for each condition
    err3m=nan(3,1);err3se=nan(3,1);
    for c=1:3
        err3m(c)=mean(error_reorg{ses}(find(condition_reorg{ses}==c)));
        err3se(c)=std(error_reorg{ses}(find(condition_reorg{ses}==c)))/sqrt(length(error_reorg{ses}(find(condition_reorg{ses}==c))));
    end
    err3m12{ses}=err3m;
    err3se12{ses}=err3se;
end

figure
clf
for ses=1:12; % pick a session
    subplot(3,4,ses)
    % X values (optional)
    x = 1:length(err3m);  % X-axis positions for each group
    % Plot mean with error bars (mean ± standard error)
    errorbar(x, err3m12{ses}./2, err3se12{ses}./2, 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red');
    xlabel('Conditions');
    ylabel('Mean Error (ms)');
    title(['Mean with Standard Error: ses' num2str(ses)]);
    grid on;
    % Set custom x-axis tick positions
    xticks([1 2 3]);xlim([0.5 3.5])
    % Set custom x-axis tick labels
    xticklabels({'uncouple', 'unidir', 'bidir'});
    % ylim([-10 400])
    if mod(ses, 2) == 0 % even
        % yline(S_offset_synco,'m-')
        subtitle('synco')
    else
        % yline(S_offset_synch,'m-')
        subtitle('synch')
    end
end


%% In one session, examine relationship btw error condi and symbol in 3d 


figure
clf
for ses=1:12
    subplot(1,12,ses)
    colors = repmat([nan nan nan], length(condition_reorg{ses}), 1);
    idx = condition_reorg{ses} == 1;
    colors(idx, :) = repmat([1 0 0], sum(idx), 1);
    idx = condition_reorg{ses} == 2;
    colors(idx, :) = repmat([0 1 0], sum(idx), 1);
    idx = condition_reorg{ses} == 3;
    colors(idx, :) = repmat([0 0 1], sum(idx), 1);
    
    scatter3(condition_reorg{ses},error_reorg{ses},joint_seq_reorg{ses}, 36, colors, 'filled')
    
    % Set custom x-axis tick positions
    xticks([1 2 3]);xlim([0.5 3.5])
    % Set custom x-axis tick labels
    xticklabels({'uncouple', 'unidir', 'bidir'});
    
    ylabel('error'); zlabel('jt synmbols')
    
    % view([0 0 1]) % view error and condi
    % view([0 1 0]) % view jt sym and condi
    view([1 0 0]) % view error and jt sym


    if mod(ses, 2) == 0 % even
        % yline(S_offset_synco,'m-')
        subtitle('synco')
    else
        % yline(S_offset_synch,'m-')
        subtitle('synch')
    end
end

% Axes handle for the legend (place it at the bottom)
axes('Position', [0.4, 0.01, 0.3, 0.08], 'Visible', 'off');

% Plot colored markers for the legend
hold on;
plot(0.30, 0.41, 'go', 'MarkerFaceColor', 'g');  
plot(0.31, 0.41, 'ro', 'MarkerFaceColor', 'r');  
plot(0.32, 0.41, 'bo', 'MarkerFaceColor', 'b');  

% Add text labels for the legend
text(0.301, 0.41, 'uncouple', 'Color', 'g', 'VerticalAlignment', 'middle');
text(0.311, 0.41, 'unidir', 'Color', 'r', 'VerticalAlignment', 'middle');
text(0.321, 0.41, 'bidir', 'Color', 'b', 'VerticalAlignment', 'middle');

%% 3d cluster
% Generate synthetic 3D data
rng(1); % For reproducibility
data1 = mvnrnd([1 2 3], eye(3), 100);
data2 = mvnrnd([5 5 5], eye(3), 100);
data3 = mvnrnd([8 1 8], eye(3), 100);
data = [data1; data2; data3];

% Number of clusters
k = 3;

% Perform k-means clustering
[idx, centroids] = kmeans(data, k);

% Create a 3D scatter plot
figure;
scatter3(data(:,1), data(:,2), data(:,3), 36, idx, 'filled');
hold on;
scatter3(centroids(:,1), centroids(:,2), centroids(:,3), 100, 'kx');
title('3D Cluster Visualization');
xlabel('X');
ylabel('Y');
zlabel('Z');
legend('Cluster 1', 'Cluster 2', 'Cluster 3', 'Centroids');
grid on;