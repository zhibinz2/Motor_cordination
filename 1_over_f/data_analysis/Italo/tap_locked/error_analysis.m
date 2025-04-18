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

% reorganize offsets from each trials
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

% concatenate offsets in 12 trials  (in samples)
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

% standard offset (in ms)
S_offset_synch=0;
S_offset_synco=2000/1.3/2/2;
% convert to relative error offset to target (in ms)
rErr=cell(12,1);
for ses=1:12
    if mod(ses, 2) == 0 % syncopation
        rErr{ses}=S_offset_synco-error_reorg{ses}./2;
    else
        rErr{ses}=error_reorg{ses}./2;
    end
end

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
set(gcf,'color','w'); % set backg

%% In one session, examine relationship btw error condi and symbol in 3d 
figure('Position',[108         101         754        1170]);  
% get(gcf, 'outerposition'); figure('Position',[ans]);  
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
    
    scatter3(condition_reorg{ses},error_reorg{ses}/2,joint_seq_reorg{ses}, 36, colors, 'o')
    
    % Set custom x-axis tick positions
    xticks([1 2 3]);xlim([0.5 3.5])
    % Set custom x-axis tick labels
    xticklabels({'uncouple', 'unidir', 'bidir'});
    
    ylabel('error (ms)'); zlabel('joint synmbols')
    
    % view([0 0 1]) % view error and condi
    % view([0 1 0]) % view jt sym and condi
    view([1 0 0]) % view error and jt sym

    title(['ses' num2str(ses)]);
    if mod(ses, 2) == 0 % even
        % yline(S_offset_synco,'m-')
        subtitle('synco')
    else
        % yline(S_offset_synch,'m-')
        subtitle('synch')
    end
end

% Axes handle for the legend (place it at the bottom)
axes('Position', [0.4, 0.01, 0.3, 0.05], 'Visible', 'off');

% Plot colored markers for the legend
hold on;
plot(0.30, 0.2, 'go', 'MarkerFaceColor', 'g');  
plot(0.31, 0.2, 'ro', 'MarkerFaceColor', 'r');  
plot(0.32, 0.2, 'bo', 'MarkerFaceColor', 'b');  

% Add text labels for the legend
text(0.301, 0.2, 'uncouple', 'Color', 'g', 'VerticalAlignment', 'middle');
text(0.311, 0.2, 'unidir', 'Color', 'r', 'VerticalAlignment', 'middle');
text(0.321, 0.2, 'bidir', 'Color', 'b', 'VerticalAlignment', 'middle');

set(gcf,'color','w'); % set backg

%% 3d cluster (example code)
ses=3
data = [condition_reorg{ses} (error_reorg{ses}/2)' (joint_seq_reorg{ses})'];

% Number of clusters
k = 3;

% Perform k-means clustering
[idx, centroids] = kmeans(double(data), k);

% Create a 3D scatter plot
figure;
clf
scatter3(data(:,1), data(:,2), data(:,3), 36, idx, 'o');
hold on;
scatter3(centroids(:,1), centroids(:,2), centroids(:,3), 100, 'ko','filled','LineWidth',3);
title('3D Cluster Visualization');
% Set custom x-axis tick positions
xticks([1 2 3]);xlim([0.5 3.5])
% Set custom x-axis tick labels
xticklabels({'uncouple', 'unidir', 'bidir'});  
ylabel('error (ms)'); zlabel('joint synmbols')
% view([0 0 1]) % view error and condi
% view([0 1 0]) % view jt sym and condi
view([1 0 0]) % view error and jt sym
title({['3D Cluster Visualization'],['ses' num2str(ses)]});
if mod(ses, 2) == 0 % even
    % yline(S_offset_synco,'m-')
    subtitle('synco')
else
    % yline(S_offset_synch,'m-')
    subtitle('synch')
end
legend('', 'Cluster Centroids','Location','north');
grid on;

[caz,cel] = view
view(caz,cel);

set(gcf,'color','w'); % set backg

% view([0 0 1]) % view error and condi
% view([0 1 0]) % view jt sym and condi
view([1 0 0]) % view error and jt sym

%% session_clusterings
% load('session_clusterings.mat');
UniqS=cell(12,2);
for ses=1:12
    varname=['ses' num2str(seeds(ses))];
    loadtmpstruct=load('session_clusterings.mat',varname);
    loadtmpfile=loadtmpstruct.(varname);
    for subj=1:2
        UniqS{ses,subj}=unique(loadtmpfile(subj,:));
    end
end

% count num of JtS in ses 12 in matrix
ses=12;
mat_JtS=zeros(length(UniqS{ses,1}),length(UniqS{ses,2}));
% assign relative offset error to JtS in matrix
ses_er=rErr{ses};
mat_sumer=zeros(length(UniqS{ses,1}),length(UniqS{ses,2}));
for indxy=1:size(loadtmpfile,2)
    mat_JtS(loadtmpfile(1,indxy)+1,loadtmpfile(2,indxy)+1)=mat_JtS(loadtmpfile(1,indxy)+1,loadtmpfile(2,indxy)+1)+1;
    % assign offset to JtS in matrix
    mat_sumer(loadtmpfile(1,indxy)+1,loadtmpfile(2,indxy)+1)=mat_sumer(loadtmpfile(1,indxy)+1,loadtmpfile(2,indxy)+1)+ses_er(indxy);
end

% jet colormap with zero in white
jet0=[[1 1 1];jet];

figure;clf
subplot(311);imagesc(mat_JtS);colorbar;title('symbols');ylabel('subj A'); xlabel('subj B');colormap(jet0);clim([0 max(mat_JtS,[],'all')])
subplot(312);imagesc(mat_sumer);colorbar;title('sum of offset');ylabel('subj A'); xlabel('subj B');colormap(jet0);clim([0 max(mat_sumer,[],'all')])
subplot(313);imagesc(mat_sumer./mat_JtS);colorbar;title('average offset (ms)');ylabel('subj A'); xlabel('subj B');colormap(jet0);clim([0 max(mat_sumer./mat_JtS,[],'all')])
sgtitle(['ses' num2str(ses)])

%  Examine relationship btw error and condi
mat4_JtS12ses=cell(12,1);
mat4_sumer12ses=cell(12,1);
cell_er12s4c=cell(12,1);
for ses=1:12; % pick a session
    % mat_JtS4c=cell(12,1);
    % mat_sumer4c=cell(12,1);
    mat_JtS=zeros(4,length(UniqS{ses,1}),length(UniqS{ses,2}));
    mat_sumer=zeros(4,length(UniqS{ses,1}),length(UniqS{ses,2}));
    cell_er4c=cell(4,1);
    ses_er=rErr{ses};
    % err3m=nan(3,length(UniqS{ses,1}),length(UniqS{ses,2}));
    % err3se=nan(3,length(UniqS{ses,1}),length(UniqS{ses,2}));
    varname=['ses' num2str(seeds(ses))];
    loadtmpstruct=load('session_clusterings.mat',varname);
    loadtmpfile=loadtmpstruct.(varname);
    for c=1:4
        indc=find(condition_reorg{ses}==c);
        cell_er=cell(length(UniqS{ses,1}),length(UniqS{ses,2}));
        for indxy=1:length(indc)
            mat_JtS(c,loadtmpfile(1,indc(indxy))+1,loadtmpfile(2,indc(indxy))+1)=...
                mat_JtS(c,loadtmpfile(1,indc(indxy))+1,loadtmpfile(2,indc(indxy))+1)+1;
            mat_sumer(c,loadtmpfile(1,indc(indxy))+1,loadtmpfile(2,indc(indxy))+1)=...
                mat_sumer(c,loadtmpfile(1,indc(indxy))+1,loadtmpfile(2,indc(indxy))+1)+ses_er(indc(indxy));
            cell_er{loadtmpfile(1,indc(indxy))+1,loadtmpfile(2,indc(indxy))+1}=...
                [cell_er{loadtmpfile(1,indc(indxy))+1,loadtmpfile(2,indc(indxy))+1} ses_er(indc(indxy))];
            % err3m(c)=mean(error_reorg{ses}(indc);
            % err3se(c)=std(error_reorg{ses}(indc))/sqrt(length(error_reorg{ses}(indc)));
        end
        cell_er4c{c}=cell_er;
    end
    mat4_JtS12ses{ses}=mat_JtS;
    mat4_sumer12ses{ses}=mat_sumer;
    cell_er12s4c{ses}=cell_er4c;
    % err3m12{ses}=err3m;
    % err3se12{ses}=err3se;
end
% examine 
ses=3; % ses=5;
mat_JtS=mat4_JtS12ses{ses};
mat_sumer=mat4_sumer12ses{ses};
cell_er4c=cell_er12s4c{ses};
c=1;
cell_er=cell_er4c{c};
for row=1:size(cell_er,1)
    row_means(row)=mean([cell_er{row,:}]);
    row_sems(row)=std([cell_er{row,:}])/sqrt(length([cell_er{row,:}]));
end
for col=1:size(cell_er,2)
    col_means(col)=mean([cell_er{:,col}]);
    col_sems(col)=std([cell_er{:,col}])/sqrt(length([cell_er{:,col}]));
end
figure;
bar(row_means);hold on;errorbar(row_means,row_sems, 'k', 'LineStyle', 'none', 'CapSize', 8);


figure;
for c=1:4
    subplot(3,4,c);imagesc(squeeze(mat_JtS(c,:,:)));colorbar;title('JtS');clim([0 50]);colormap(jet0);
    subplot(3,4,4+c);imagesc(squeeze(mat_sumer(c,:,:)));colorbar; title('sum of err');clim([0 10000]);colormap(jet0);
    subplot(3,4,8+c);imagesc(squeeze(mat_sumer(c,:,:))./squeeze(mat_JtS(c,:,:)));colorbar; title('mean of error (ms)');clim([0 S_offset_synco]);colormap(jet0);
end

for ses=1:12
    mat_JtS=mat4_JtS12ses{ses};
    mat_sumer=mat4_sumer12ses{ses};
    figure;clf
    condi4names={'uncouple','A leading','B leading','bidirectional'}
    for c=1:4
        subplot(2,4,c);imagesc(squeeze(mat_JtS(c,:,:)));colorbar;title('Symbols');clim([0 50]);
        subtitle(condi4names{c});colormap(jet0)
        if c==2
            xlabel('leader');ylabel('follower');
        elseif c==3
            ylabel('leader');xlabel('follower');
        else
            xlabel('participant A'); ylabel('participant B');
        end  
        subplot(2,4,4+c);imagesc(squeeze(mat_sumer(c,:,:))./squeeze(mat_JtS(c,:,:)));colorbar; title('mean of error (ms)');clim([0 S_offset_synco]);
        subtitle(condi4names{c});colormap(jet0);
        if c==2
            xlabel('leader');ylabel('follower');
        elseif c==3
            ylabel('leader');xlabel('follower');
        else
            xlabel('participant A'); ylabel('participant B');
        end  
    end
    sgtitle(['session ' num2str(ses) ])
    set(gcf,'color','w'); % set backg
end

% combine leaders and followers not possible
mat3_JtS12ses=cell(12,1);
mat3_sumer12ses=cell(12,1);
for ses=1:12; % pick a session
    for c=1:4
        if c==3
            mat_JtS=squeeze(mat4_JtS12ses{ses}(c,:,:))';
            mat_sumer=squeeze(mat4_sumer12ses{ses}(c,:,:))';
        else
            mat_JtS=squeeze(mat4_JtS12ses{ses}(c,:,:));
            mat_sumer=squeeze(mat4_sumer12ses{ses}(c,:,:));
        end
    end
end


%% EEG complexity
cd /ssd/zhibin/1overf/Cleaned_sourcedata/cortical_source_data/taplocked/corr_matrices/
% covert to mat files
tic
load('dyad_20220713_matstack.mat')
toc % 1min
subj1tapnum'
subj2tapnum'
subj1filenames_cell_array'
subj2filenames_cell_array'
size(subj1mat)
size(subj2mat)

% test H
i=180;j=300;taps=1:917;
corr_series=subj1mat(taps,i,j);

% add this repo to your path (https://github.com/wodeyara/AdaptiveGraphicalLassoforParCoh)
% addpath /home/zhibinz2/Documents/GitHub/AdaptiveGraphicalLassoforParCoh/Simulations/util
addpath(genpath('/home/zhibinz2/Documents/GitHub/AdaptiveGraphicalLassoforParCoh'))
mat=squeeze(subj1mat(1,:,:));
% convert real value covariance to complex matrix (448x448) 
% then compute coherence
Pcoh=normalizeCSD(mat);

tic
coh=nan(size(subj1mat));
for pc=1:size(subj1mat)
    coh(pc,:,:)=normalizeCSD(squeeze(subj1mat(pc,:,:)));
end
toc % 20s 
corr_series=coh(taps,i,j); % lots of inf , 

pc=150;
imagesc(squeeze(coh(pc,:,:)));colorbar

addpath(genpath('/home/zhibinz2/Documents/GitHub/Motor_coordination_code/util'))
[~,H]=DFA_main(corr_series);



%% imagesc with distribution of x y summation
ses=1
c=1


% Create a figure
for ses =1:12
    mat_JtS=mat4_JtS12ses{ses};
    
    for c=1:4
        mat_temp=squeeze(mat_JtS(c,:,:));
        row_sum = sum(mat_temp, 2);  % Sum along each row 
        col_sum = sum(mat_temp, 1);  % Sum along each column 
        display([ condi4names{c} '  session ' num2str(ses)]);

        figure('Position',[2789         390        1480         796]); % get(gcf, 'outerposition')
        clf

        % sort and plot error first
        subplot(3,6,[10,11,16,17]);

        cell_er4c=cell_er12s4c{ses};
        cell_er=cell_er4c{c};
        clear row_means row_sems col_means col_sems
        for row=1:size(cell_er,1)
            row_means(row)=mean([cell_er{row,:}]);
            row_sems(row)=std([cell_er{row,:}])/sqrt(length([cell_er{row,:}]));
        end
        for col=1:size(cell_er,2)
            col_means(col)=mean([cell_er{:,col}]);
            col_sems(col)=std([cell_er{:,col}])/sqrt(length([cell_er{:,col}]));
        end

        % Step 2: Sort rows and columns based on their sums
        [~, rowOrder] = sort(row_means); % Indices for sorted rows
        [~, colOrder] = sort(col_means); % Indices for sorted columns

        mat_sumer=mat4_sumer12ses{ses};
        mater_temp=squeeze(mat_sumer(c,:,:))./mat_temp;

        % Reorganize the matrix based on sorted sums
        sortedMatrix = mater_temp(rowOrder, colOrder);

        imagesc(sortedMatrix);
        yticks([1:size(sortedMatrix,1)]);xticks([1:size(sortedMatrix,2)]);
        yticklabels(rowOrder-1);xticklabels(colOrder-1);

        clim([0 S_offset_synco]);colormap(jet0)
        % colorbar('westoutside');
        % set(gca, 'XAxisLocation', 'top', 'YAxisLocation', 'left');  % Place x-axis on top
        title({['mean Error (ms)'],['(' condi4names{c} ')'], ['session ' num2str(ses)]},'Rotation',0,'Position',[size(mat_temp,2)+3,-2,0]);
        cb = colorbar('southoutside');cb.Position = [0.54, 0.025, 0.23, 0.01];
        
        if c==2
            xlabel('mean Error (leader)');ylabel(' mean Error (follower)');
        elseif c==3
            ylabel('mean Error (leader)');xlabel('mean Error (follower)');
        else
            xlabel('mean Error (participant A)'); ylabel('mean Error (participant B)');
        end  

        
        subplot(3,6,[12,18]);  % left col
        bar(1:size(mat_temp,1),row_means(rowOrder), 'FaceColor', grey);  % Horizontal bar plot for row sums;
        xlim([0.5 size(mat_temp,1)+0.5]); title(['Mean of Rows']); ylabel([]);xticks([]);
        hold on; 
        % set(gca, 'XDir', 'reverse');
        view(90, 90);
        errorbar(1:size(mat_temp,1),row_means(rowOrder), row_sems(rowOrder), 'k', 'LineStyle', 'none','CapSize', 0);
        hold off;
        % set(gca, 'YAxisLocation', 'right', 'XAxisLocation', 'top');  % Move X-axis to the top
        % set(gca, 'XAxisLocation', 'left', 'YAxisLocation', 'origin');  % Place x-axis on top
        ylim([0 300])
        
        subplot(3,6,[4,5]);  % row bottom
        bar(1:size(mat_temp,2),col_means(colOrder), 'FaceColor', grey);  % Bar plot for column sums
        hold on; errorbar(1:size(mat_temp,2),col_means(colOrder), col_sems(colOrder), 'k', 'LineStyle', 'none','CapSize', 0);
        hold off;
        xlim([0.5 size(mat_temp,2)+0.5]);
        xticks([]);xlabel([]);
        title(['Mean of Columns']);
        set(gca,'XAxisLocation', 'top', 'YAxisLocation', 'left');  % Place x-axis on top
        ylim([0 300])
        
        % same order for the symbols
        subplot(3,6,[7,8,13,14]);
        % Reorganize the matrix based on sorted sums
        sortedMatrix = mat_temp(rowOrder, colOrder);
        imagesc(sortedMatrix);clim([0 50]);colormap(jet0)
        yticks([1:size(sortedMatrix,1)]);xticks([1:size(sortedMatrix,2)]);
        yticklabels(rowOrder-1);xticklabels(colOrder-1);
        % colorbar('westoutside');
        % set(gca, 'XAxisLocation', 'top', 'YAxisLocation', 'left');  % Place x-axis on top
        title({['Symbols'],['(' condi4names{c} ')'], ['session ' num2str(ses)]},'Rotation',0,'Position',[size(mat_temp,2)+3,-2,0]);
        cb = colorbar('southoutside');cb.Position = [0.13, 0.025, 0.23, 0.01];
        
        if c==2
            xlabel('Symbols (leader)');ylabel('Symbols (follower)');
        elseif c==3
            ylabel('Symbols (leader)');xlabel('Symbols (follower)');
        else
            xlabel('Symbols (participant A)'); ylabel('Symbols (participant B)');
        end  
        
        subplot(3,6,[9,15]);  % left col
        barh(1:size(mat_temp,1),row_sum(rowOrder), 'FaceColor', grey);  % Horizontal bar plot for row sums;
        ylim([0.5 size(mat_temp,1)+0.5]);title(['Sum of Rows']);ylabel([]);yticks([]);
        set(gca, 'YDir', 'reverse');
        % set(gca, 'YAxisLocation', 'right', 'XAxisLocation', 'top');  % Move X-axis to the top
        set(gca, 'XAxisLocation', 'bottom', 'YAxisLocation', 'right');  % Place x-axis on top
        xlim([0 250])
        
        subplot(3,6,[1,2]);  % row bottom
        bar(1:size(mat_temp,2),col_sum(colOrder), 'FaceColor', grey);  % Bar plot for column sums
        xlim([0.5 size(mat_temp,2)+0.5]);xticks([]);xlabel([]);
        title(['Sum of Columns']);
        set(gca,'XAxisLocation', 'top', 'YAxisLocation', 'left');  % Place x-axis on top
        ylim([0 250])

        set(gcf,'color','w'); % set backg
        
        exportgraphics(gcf, 'multiple_figures.pdf', 'Append', true); % Append to PDF
        close all;
    end
end

%% sort order example
% Example matrix
A = [9 2 1; 
     8 7 6; 
     4 3 5];
figure;
subplot(121)
imagesc(A);colorbar;colormap(jet)
yticks([1:size(A,1)]);xticks([1:size(A,2)]);
yticklabels(yticks-1);xticklabels(xticks-1);

% Step 1: Compute row and column sums
rowSums = sum(A, 2);  % Sum of each row
colSums = sum(A, 1);  % Sum of each column

% Step 2: Sort rows and columns based on their sums
[~, rowOrder] = sort(rowSums); % Indices for sorted rows
[~, colOrder] = sort(colSums); % Indices for sorted columns

% Step 3: Reorganize the matrix based on sorted sums
sortedMatrix = A(rowOrder, colOrder);
subplot(122)
imagesc(sortedMatrix);colorbar;colormap(jet)
yticks([1:size(A,1)]);xticks([1:size(A,2)]);
yticklabels(rowOrder-1);xticklabels(colOrder-1);

% Display results
disp('Original Matrix:');
disp(A);
disp('Reorganized Matrix:');
disp(sortedMatrix);



%% give up
figure;
clf
ax1 = axes('Position', [0.1, 0.1, 0.4, 0.8]);imagesc(ax1, rand(10, 10));
ax2 = axes('Position', [0.55, 0.1, 0.4, 0.8]);imagesc(ax2, rand(10, 10));

for ses=1:12
    mat_JtS=mat4_JtS12ses{ses};
    mat_sumer=mat4_sumer12ses{ses};
    figure;clf
    condi4names={'uncouple','A leading','B leading','bidirectional'}
    for c=1:4
        mat_temp=squeeze(mat_JtS(c,:,:));
        row_sum = sum(mat_temp, 2);  % Sum along each row 
        col_sum = sum(mat_temp, 1);  % Sum along each column 
        
        imagesc(axes('Position', [0.05, 0.05, 0.15, 0.3]),mat_temp);colorbar('southoutside');clim([0 50]);
        title(['Symbols (' condi4names{c} ')'],'Rotation',90,'Position',[-0.75,size(mat_temp,1)*4/5,0]);
        colormap(jet0)
        if c==2
            xlabel('leader');ylabel('follower');
        elseif c==3
            ylabel('leader');xlabel('follower');
        else
            xlabel('participant A'); ylabel('participant B');
        end  

        barh(axes('Position', [0.2, 0.18, 0.05, 0.163]),row_sum, 'FaceColor', 'r'); 
        ylim([0.5 size(mat_temp,1)+0.5]);yticks([]);


        subplot(2,4,4+c);imagesc(squeeze(mat_sumer(c,:,:))./squeeze(mat_JtS(c,:,:)));colorbar; title('mean of error (ms)');clim([0 S_offset_synco]);
        subtitle(condi4names{c});colormap(jet0);
        if c==2
            xlabel('leader');ylabel('follower');
        elseif c==3
            ylabel('leader');xlabel('follower');
        else
            xlabel('participant A'); ylabel('participant B');
        end  
    end
    sgtitle(['session ' num2str(ses) ])
    set(gcf,'color','w'); % set backg
    
end

%% load conditions and H
open /home/zhibinz2/Documents/GitHub/Motor_coordination_code/Fig3_H_matching.m
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
H_all;
% compute H difference
Hdiff=nan(12,12);
for ses=1:12
    for tr=1:12
        Hdiff(ses,tr)=abs(H_all(1,ses,tr)-H_all(2,ses,tr));
    end
end


condition_all;


% load average dwell time and motif length
load('dict_result.mat');
Average_Dwell_Time_unsorted=nan(144,1);
Average_Motif_Length_unsorted=nan(144,1);
% Condition_unsorted=cell(144,1);
% Session_unsorted=cell(144,1);
% Session_Type_unsorted=cell(144,1);
for i=1:144
    Average_Dwell_Time_unsorted(i)=cell_data.Average_Dwell_Time{i};
    Average_Motif_Length_unsorted(i)=cell_data.Average_Motif_Length{i};
end
Condition_unsorted=cell_data.Condition;
Session_unsorted=cell_data.Session;
Session_Type_unsorted=cell_data.Session_Type;


% sort in time order
Average_Dwell_Time=nan(12,12);
Average_Motif_Length=nan(12,12);
Conditions=cell(12,12);
condition_id=nan(12,12);
sessions_seeds=cell(12,12);
sessions_type=cell(12,12);
sessions_id=nan(12,12);

for ses=1:12
    idx=find(strcmp(Session_unsorted,num2str(seeds(ses,:))));
    Average_Dwell_Time(ses,:)=Average_Dwell_Time_unsorted(idx)';
    Average_Motif_Length(ses,:)=Average_Motif_Length_unsorted(idx,:)';
    condition_id(ses,:)=condition_all(12*(ses-1)+(1:12))';
    for tr=1:12
        Conditions{ses,tr}=Condition_unsorted{idx(tr)};
        sessions_seeds{ses,tr}=Session_unsorted{idx(tr)};
        sessions_type{ses,tr}=Session_Type_unsorted{idx(tr)};
        if strcmp(Session_Type_unsorted{idx(tr)},'Syncopation')
            sessions_id(ses,tr)=2;
        else
            sessions_id(ses,tr)=1;
        end
    end
end

sessions_type144=cell(144,1);
Conditions144=cell(144,1);
Average_Dwell_Time144=nan(144,1);
Average_Motif_Length144=nan(144,1);
for ses=1:12
    for tr=1:12
        sessions_type144{12*(ses-1)+tr}=sessions_type{ses,tr};
        Conditions144{12*(ses-1)+tr}=Conditions{ses,tr};
        Average_Dwell_Time144(12*(ses-1)+tr)=Average_Dwell_Time(ses,tr);
        Average_Motif_Length144(12*(ses-1)+tr)=Average_Motif_Length(ses,tr);
    end
end

% find index test
indx=intersect(find(condition_id==4),find(sessions_id==2));
% 'Uncoupled','Leader-Follower','Mutual'
indx=intersect(find(strcmpi(Conditions,'Mutual')),find(strcmpi(sessions_type,'Syncopation')));

% compute mean and sem of Average_Dwell_Time and Average_Motif_Length and Hdif 2x3
dwell_means=nan(2,3);dwell_sems=nan(2,3);
motif_means=nan(2,3);motif_sems=nan(2,3);

indx=intersect(find(strcmpi(Conditions,'Uncoupled')),find(strcmpi(sessions_type,'Synchronization')));
dw_gp11=[Average_Dwell_Time(indx)];
dwell_means(1,1)=mean(dw_gp11);
dwell_sems(1,1)=std(dw_gp11)/sqrt(length(dw_gp11));
mt_gp11=[Average_Motif_Length(indx)];
motif_means(1,1)=mean(mt_gp11);
motif_sems(1,1)=std(mt_gp11)/sqrt(length(mt_gp11));

indx=intersect(find(strcmpi(Conditions,'Leader-Follower')),find(strcmpi(sessions_type,'Synchronization')));
dw_gp12=[Average_Dwell_Time(indx)];
dwell_means(1,2)=mean(dw_gp12);
dwell_sems(1,2)=std(dw_gp12)/sqrt(length(dw_gp12));
mt_gp12=[Average_Motif_Length(indx)];
motif_means(1,2)=mean(mt_gp12);
motif_sems(1,2)=std(mt_gp12)/sqrt(length(mt_gp12));

indx=intersect(find(strcmpi(Conditions,'Mutual')),find(strcmpi(sessions_type,'Synchronization')));
dw_gp13=[Average_Dwell_Time(indx)];
dwell_means(1,3)=mean(dw_gp13);
dwell_sems(1,3)=std(dw_gp13)/sqrt(length(dw_gp13));
mt_gp13=[Average_Motif_Length(indx)];
motif_means(1,3)=mean(mt_gp13);
motif_sems(1,3)=std(mt_gp13)/sqrt(length(mt_gp13));

indx=intersect(find(strcmpi(Conditions,'Uncoupled')),find(strcmpi(sessions_type,'Syncopation')));
dw_gp21=[Average_Dwell_Time(indx)];
dwell_means(2,1)=mean(dw_gp21);
dwell_sems(2,1)=std(dw_gp21)/sqrt(length(dw_gp21));
mt_gp21=[Average_Motif_Length(indx)];
motif_means(2,1)=mean(mt_gp21);
motif_sems(2,1)=std(mt_gp21)/sqrt(length(mt_gp21));

indx=intersect(find(strcmpi(Conditions,'Leader-Follower')),find(strcmpi(sessions_type,'Syncopation')));
dw_gp22=[Average_Dwell_Time(indx)];
dwell_means(2,2)=mean(dw_gp22);
dwell_sems(2,2)=std(dw_gp22)/sqrt(length(dw_gp22));
mt_gp22=[Average_Motif_Length(indx)];
motif_means(2,2)=mean(mt_gp22);
motif_sems(2,2)=std(mt_gp22)/sqrt(length(mt_gp22));

indx=intersect(find(strcmpi(Conditions,'Mutual')),find(strcmpi(sessions_type,'Syncopation')));
dw_gp23=[Average_Dwell_Time(indx)];
dwell_means(2,3)=mean(dw_gp23);
dwell_sems(2,3)=std(dw_gp23)/sqrt(length(dw_gp23));
mt_gp23=[Average_Motif_Length(indx)];
motif_means(2,3)=mean(mt_gp23);
motif_sems(2,3)=std(mt_gp23)/sqrt(length(mt_gp23));

anova_dwell_synch={}

figure
subplot(221)
bar(dwell_means(1,:));hold on;errorbar(dwell_means(1,:),dwell_sems(1,:), 'k', 'LineStyle', 'none', 'CapSize', 8);
title('synch');ylabel('mean dwell time');xlabel('conditions');xticklabels({'uncouple', 'unidir', 'bidir'});
subplot(222)
bar(dwell_means(2,:));hold on;errorbar(dwell_means(2,:),dwell_sems(2,:), 'k', 'LineStyle', 'none', 'CapSize', 8);
title('synco');ylabel('mean dwell time');xlabel('conditions');xticklabels({'uncouple', 'unidir', 'bidir'});
subplot(223)
bar(motif_means(1,:));hold on;errorbar(motif_means(1,:),motif_sems(1,:), 'k', 'LineStyle', 'none', 'CapSize', 8);
title('synch');ylabel('mean motif length');xlabel('conditions');xticklabels({'uncouple', 'unidir', 'bidir'});
subplot(224)
bar(motif_means(2,:));hold on;errorbar(motif_means(2,:),motif_sems(2,:), 'k', 'LineStyle', 'none', 'CapSize', 8);
title('synco');ylabel('mean motif length');xlabel('conditions');xticklabels({'uncouple', 'unidir', 'bidir'});

% anova statistics and plot
group1=dw_gp21'; group2=dw_gp22'; group3=dw_gp23';
group1=mt_gp11'; group2=mt_gp12'; group3=mt_gp13';
% Combine the data into a single vector
data = [group1, group2, group3];
% Create a grouping variable
group = [ones(size(group1)), 2 * ones(size(group2)), 3 * ones(size(group3))];

open /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/Italo/tap_locked/stat_test.m
sgtitle('mean $\pm$ standard error', 'Interpreter', 'latex')

%% paired t-test

% compute mean and sem of Average_Dwell_Time and Average_Motif_Length and Hdif 2x3
dwell_means=nan(2,3);dwell_sems=nan(2,3);
motif_means=nan(2,3);motif_sems=nan(2,3);

indx=intersect(find(strcmpi(Conditions144,'Uncoupled')),find(strcmpi(sessions_type144,'Synchronization')));
dw_gp11=[Average_Dwell_Time144(indx)];
dwell_means(1,1)=mean(dw_gp11);
dwell_sems(1,1)=std(dw_gp11)/sqrt(length(dw_gp11));
mt_gp11=[Average_Motif_Length144(indx)];
motif_means(1,1)=mean(mt_gp11);
motif_sems(1,1)=std(mt_gp11)/sqrt(length(mt_gp11));

indx=intersect(find(strcmpi(Conditions144,'Leader-Follower')),find(strcmpi(sessions_type144,'Synchronization')));
dw_gp12=[Average_Dwell_Time144(indx)];
dwell_means(1,2)=mean(dw_gp12);
dwell_sems(1,2)=std(dw_gp12)/sqrt(length(dw_gp12));
mt_gp12=[Average_Motif_Length144(indx)];
motif_means(1,2)=mean(mt_gp12);
motif_sems(1,2)=std(mt_gp12)/sqrt(length(mt_gp12));

indx=intersect(find(strcmpi(Conditions144,'Mutual')),find(strcmpi(sessions_type144,'Synchronization')));
dw_gp13=[Average_Dwell_Time144(indx)];
dwell_means(1,3)=mean(dw_gp13);
dwell_sems(1,3)=std(dw_gp13)/sqrt(length(dw_gp13));
mt_gp13=[Average_Motif_Length144(indx)];
motif_means(1,3)=mean(mt_gp13);
motif_sems(1,3)=std(mt_gp13)/sqrt(length(mt_gp13));

indx=intersect(find(strcmpi(Conditions144,'Uncoupled')),find(strcmpi(sessions_type144,'Syncopation')));
dw_gp21=[Average_Dwell_Time144(indx)];
dwell_means(2,1)=mean(dw_gp21);
dwell_sems(2,1)=std(dw_gp21)/sqrt(length(dw_gp21));
mt_gp21=[Average_Motif_Length144(indx)];
motif_means(2,1)=mean(mt_gp21);
motif_sems(2,1)=std(mt_gp21)/sqrt(length(mt_gp21));

indx=intersect(find(strcmpi(Conditions144,'Leader-Follower')),find(strcmpi(sessions_type144,'Syncopation')));
dw_gp22=[Average_Dwell_Time144(indx)];
dwell_means(2,2)=mean(dw_gp22);
dwell_sems(2,2)=std(dw_gp22)/sqrt(length(dw_gp22));
mt_gp22=[Average_Motif_Length144(indx)];
motif_means(2,2)=mean(mt_gp22);
motif_sems(2,2)=std(mt_gp22)/sqrt(length(mt_gp22));

indx=intersect(find(strcmpi(Conditions144,'Mutual')),find(strcmpi(sessions_type144,'Syncopation')));
dw_gp23=[Average_Dwell_Time144(indx)];
dwell_means(2,3)=mean(dw_gp23);
dwell_sems(2,3)=std(dw_gp23)/sqrt(length(dw_gp23));
mt_gp23=[Average_Motif_Length144(indx)];
motif_means(2,3)=mean(mt_gp23);
motif_sems(2,3)=std(mt_gp23)/sqrt(length(mt_gp23));

% paired t-test
group1=dw_gp21'; group2=dw_gp23';
group1=mt_gp11'; group2=mt_gp13';