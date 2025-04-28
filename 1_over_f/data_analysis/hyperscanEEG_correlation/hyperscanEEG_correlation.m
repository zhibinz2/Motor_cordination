open /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/Publish_figures/H_matching.m
% run session 1 up to Bonaza session


%% load EEG pcorr

cd /home/zhibinz2/Documents/GitHub/Motor_coordination_code/
open Fig6_networkx.ipynb
% using partial correlation from Corr_eeg144
load('Corr_eeg144.mat') 
% originally from 
% cd /home/zhibinz2/Documents/GitHub/1overf/Corr_amplitude
% /home/zhibinz2/Documents/GitHub/1overf/Corr_amplitude/step1_ouput_Corr_filter_hilbert_fast_copy.ipynb
eeg_empirical_correlation144x2; % 144x2x7x32x32
eeg_partial_correlation144x2;
intervals144;
samples144;
session144;

eeg_partial_correlation144x2=logical(eeg_partial_correlation144x2);

% open/home/zhibinz2/Documents/GitHub/1overf/Corr_amplitude/step1_ouput_Corr_filter_hilbert_fast_copy.ipynb
% open /home/zhibinz2/Documents/GitHub/Motor_coordination_code/Fig3_H_matching.m
% /home/zhibinz2/Documents/GitHub/1overf/Corr_amplitude/step3_organize_indicies.ipynb
load('Indicies.mat') % this is in time sequence
Uncoupled_synch_ind =   Uncoupled_synch_ind+1;
L_lead_synch_ind    =      L_lead_synch_ind+1;
R_lead_synch_ind    =      R_lead_synch_ind+1;
Mutual_synch_ind    =      Mutual_synch_ind+1;
synch_condi4_ind    = [synch_condi4_ind{1}+1;synch_condi4_ind{2}+1;synch_condi4_ind{3}+1;synch_condi4_ind{4}+1];

Uncoupled_synco_ind =   Uncoupled_synco_ind+1;
L_lead_synco_ind    =      L_lead_synco_ind+1;
R_lead_synco_ind    =      R_lead_synco_ind+1;
Mutual_synco_ind    =      Mutual_synco_ind+1;
synco_condi4_ind    = [synco_condi4_ind{1}+1;synco_condi4_ind{2}+1;synco_condi4_ind{3}+1;synco_condi4_ind{4}+1];

condi4_ind      = [condi4_ind{1}{1}+1;condi4_ind{2}{1}+1;condi4_ind{3}{1}+1;condi4_ind{4}{1}+1];

syn2_condi4_ind = cat(3,synch_condi4_ind,synco_condi4_ind)

Uncoupled_ind   = Uncoupled_ind{1}+1;
L_lead_ind      = L_lead_ind{1}+1;
R_lead_ind      = R_lead_ind{1}+1;
Mutual_ind      = Mutual_ind{1}+1;

synch_ind       = synch_ind{1}+1;
synco_ind       = synco_ind{1}+1;
syn_ind         = [syn_ind{1}{1}+1; syn_ind{2}{1}+1];

% conditional organization
A_follow_synco = zeros(36,7,32,32);
A_lead_synco = zeros(36,7,32,32);
A_ind_synco = zeros(36,7,32,32);
A_mutual_synco = zeros(36,7,32,32);
A_follow_synch = zeros(36,7,32,32);
A_lead_synch = zeros(36,7,32,32);
A_ind_synch = zeros(36,7,32,32);
A_mutual_synch = zeros(36,7,32,32);

% Assuming indices is a struct and eeg_partial_correlation144x2subj is a 5-D array in MATLAB
A_follow_synco(1:18,:,:,:)  = eeg_partial_correlation144x2(R_lead_synco_ind,1,:,:,:);
A_follow_synco(19:36,:,:,:) = eeg_partial_correlation144x2(L_lead_synco_ind,2,:,:,:);
A_lead_synco(1:18,:,:,:)    = eeg_partial_correlation144x2(L_lead_synco_ind,1,:,:,:);
A_lead_synco(19:36,:,:,:)   = eeg_partial_correlation144x2(R_lead_synco_ind,2,:,:,:);
A_mutual_synco(1:18,:,:,:)  = eeg_partial_correlation144x2(Mutual_synco_ind,1,:,:,:);
A_mutual_synco(19:36,:,:,:) = eeg_partial_correlation144x2(Mutual_synco_ind,2,:,:,:);
A_ind_synco(1:18,:,:,:)     = eeg_partial_correlation144x2(Uncoupled_synco_ind,1,:,:,:);
A_ind_synco(19:36,:,:,:)    = eeg_partial_correlation144x2(Uncoupled_synco_ind,2,:,:,:);
A_follow_synch(1:18,:,:,:)  = eeg_partial_correlation144x2(R_lead_synch_ind,1,:,:,:);
A_follow_synch(19:36,:,:,:) = eeg_partial_correlation144x2(L_lead_synch_ind,2,:,:,:);
A_lead_synch(1:18,:,:,:)    = eeg_partial_correlation144x2(L_lead_synch_ind,1,:,:,:);
A_lead_synch(19:36,:,:,:)   = eeg_partial_correlation144x2(R_lead_synch_ind,2,:,:,:);
A_mutual_synch(1:18,:,:,:)  = eeg_partial_correlation144x2(Mutual_synch_ind,1,:,:,:);
A_mutual_synch(19:36,:,:,:) = eeg_partial_correlation144x2(Mutual_synch_ind,2,:,:,:);
A_ind_synch(1:18,:,:,:)     = eeg_partial_correlation144x2(Uncoupled_synch_ind,1,:,:,:);
A_ind_synch(19:36,:,:,:)    = eeg_partial_correlation144x2(Uncoupled_synch_ind,2,:,:,:);


A_syn=nan(2,4,36,7,32,32);
% organized into A_syn: 2 syn x 4 condi x 36 tr x 7 freq x 32 chan x 32 chan
A_syn(1,1,:,:,:,:)=A_ind_synch;
A_syn(1,2,:,:,:,:)=A_lead_synch;
A_syn(1,3,:,:,:,:)=A_follow_synch;
A_syn(1,4,:,:,:,:)=A_mutual_synch;
A_syn(2,1,:,:,:,:)=A_ind_synco;
A_syn(2,2,:,:,:,:)=A_lead_synco;
A_syn(2,3,:,:,:,:)=A_follow_synco;
A_syn(2,4,:,:,:,:)=A_mutual_synco;

% load Pcorr with conditional organization
% /home/zhibinz2/Documents/GitHub/Motor_coordination_code/Fig6_networkx.ipynb
load('Asyn.mat')
A_syn_test=nan(2,4,36,7,32,32);
A_synch=Asyn{1};A_synco=Asyn{2}; % same sequeence 
for syn = 1:2
    for st = 1:4
        A_syn_test(syn,st,:,:,:,:)=Asyn{syn}{st}; % same sequeence, checked
    end
end
% check
for syn=1:2
    for st =1:4
        display(['syn ' num2str(syn) ' st ' num2str(st)])
        sum(squeeze(A_syn_test(syn,st,:,:,:,:)-A_syn(syn,st,:,:,:,:)),'all')
    end
end
sum(squeeze(A_syn_test(1,1,:,:,:,:)-A_syn(1,1,:,:,:,:)),'all')

%% reorganzie (degree ctr, eff, bw)
cd /home/zhibinz2/Documents/GitHub/Motor_coordination_code/
load('nx144.mat'); % 144x2x7x32
dg_ctr144;
efficiency144;
bw_ctr144;

Uncoupled_synco_ind;
L_lead_synco_ind   ;
R_lead_synco_ind   ;
Mutual_synco_ind   ;

Uncoupled_synch_ind;
L_lead_synch_ind   ;
R_lead_synch_ind   ;
Mutual_synch_ind   ;

% organize from 144 tr x 2 subj x 7 freq x 32 chan to 2 syn x 4 condi x 36 tr x 7fre x 32 chan
dgctr_follow_synco  = zeros(36,7,32);
dgctr_lead_synco    = zeros(36,7,32);
dgctr_ind_synco     = zeros(36,7,32);
dgctr_mutual_synco  = zeros(36,7,32);
dgctr_follow_synch  = zeros(36,7,32);
dgctr_lead_synch    = zeros(36,7,32);
dgctr_ind_synch     = zeros(36,7,32);
dgctr_mutual_synch  = zeros(36,7,32);
dgctr_follow_synco(1:18,:,:)  = squeeze(dg_ctr144(R_lead_synco_ind,1,:,:));
dgctr_follow_synco(19:36,:,:) = squeeze(dg_ctr144(L_lead_synco_ind,2,:,:));
dgctr_lead_synco(1:18,:,:)    = squeeze(dg_ctr144(L_lead_synco_ind,1,:,:));
dgctr_lead_synco(19:36,:,:)   = squeeze(dg_ctr144(R_lead_synco_ind,2,:,:));
dgctr_mutual_synco(1:18,:,:)  = squeeze(dg_ctr144(Mutual_synco_ind,1,:,:));
dgctr_mutual_synco(19:36,:,:) = squeeze(dg_ctr144(Mutual_synco_ind,2,:,:));
dgctr_ind_synco(1:18,:,:)     = squeeze(dg_ctr144(Uncoupled_synco_ind,1,:,:));
dgctr_ind_synco(19:36,:,:)    = squeeze(dg_ctr144(Uncoupled_synco_ind,2,:,:));
dgctr_follow_synch(1:18,:,:)  = squeeze(dg_ctr144(R_lead_synch_ind,1,:,:));
dgctr_follow_synch(19:36,:,:) = squeeze(dg_ctr144(L_lead_synch_ind,2,:,:));
dgctr_lead_synch(1:18,:,:)    = squeeze(dg_ctr144(L_lead_synch_ind,1,:,:));
dgctr_lead_synch(19:36,:,:)   = squeeze(dg_ctr144(R_lead_synch_ind,2,:,:));
dgctr_mutual_synch(1:18,:,:)  = squeeze(dg_ctr144(Mutual_synch_ind,1,:,:));
dgctr_mutual_synch(19:36,:,:) = squeeze(dg_ctr144(Mutual_synch_ind,2,:,:));
dgctr_ind_synch(1:18,:,:)     = squeeze(dg_ctr144(Uncoupled_synch_ind,1,:,:));
dgctr_ind_synch(19:36,:,:)    = squeeze(dg_ctr144(Uncoupled_synch_ind,2,:,:));
dgctr_synch=zeros(4,36,7,32); % 4condi x 36tr x 7freq x 32chan
dgctr_synco=zeros(4,36,7,32);
dgctr_synch(1,:,:,:)=dgctr_ind_synch;
dgctr_synch(2,:,:,:)=dgctr_lead_synch;
dgctr_synch(3,:,:,:)=dgctr_follow_synch;
dgctr_synch(4,:,:,:)=dgctr_mutual_synch;
dgctr_synco(1,:,:,:)=dgctr_ind_synco;
dgctr_synco(2,:,:,:)=dgctr_lead_synco;
dgctr_synco(3,:,:,:)=dgctr_follow_synco;
dgctr_synco(4,:,:,:)=dgctr_mutual_synco;
dgctr_syn=zeros(2,4,36,7,32); % 2syn x 4condi x 36tr x 7freq x 32chan
dgctr_syn(1,:,:,:,:)=dgctr_synch;
dgctr_syn(2,:,:,:,:)=dgctr_synco;

effi_follow_synco  = zeros(36,7,32);
effi_lead_synco    = zeros(36,7,32);
effi_ind_synco     = zeros(36,7,32);
effi_mutual_synco  = zeros(36,7,32);
effi_follow_synch  = zeros(36,7,32);
effi_lead_synch    = zeros(36,7,32);
effi_ind_synch     = zeros(36,7,32);
effi_mutual_synch  = zeros(36,7,32);
effi_follow_synco(1:18,:,:)  = squeeze(efficiency144(R_lead_synco_ind,1,:,:));
effi_follow_synco(19:36,:,:) = squeeze(efficiency144(L_lead_synco_ind,2,:,:));
effi_lead_synco(1:18,:,:)    = squeeze(efficiency144(L_lead_synco_ind,1,:,:));
effi_lead_synco(19:36,:,:)   = squeeze(efficiency144(R_lead_synco_ind,2,:,:));
effi_mutual_synco(1:18,:,:)  = squeeze(efficiency144(Mutual_synco_ind,1,:,:));
effi_mutual_synco(19:36,:,:) = squeeze(efficiency144(Mutual_synco_ind,2,:,:));
effi_ind_synco(1:18,:,:)     = squeeze(efficiency144(Uncoupled_synco_ind,1,:,:));
effi_ind_synco(19:36,:,:)    = squeeze(efficiency144(Uncoupled_synco_ind,2,:,:));
effi_follow_synch(1:18,:,:)  = squeeze(efficiency144(R_lead_synch_ind,1,:,:));
effi_follow_synch(19:36,:,:) = squeeze(efficiency144(L_lead_synch_ind,2,:,:));
effi_lead_synch(1:18,:,:)    = squeeze(efficiency144(L_lead_synch_ind,1,:,:));
effi_lead_synch(19:36,:,:)   = squeeze(efficiency144(R_lead_synch_ind,2,:,:));
effi_mutual_synch(1:18,:,:)  = squeeze(efficiency144(Mutual_synch_ind,1,:,:));
effi_mutual_synch(19:36,:,:) = squeeze(efficiency144(Mutual_synch_ind,2,:,:));
effi_ind_synch(1:18,:,:)     = squeeze(efficiency144(Uncoupled_synch_ind,1,:,:));
effi_ind_synch(19:36,:,:)    = squeeze(efficiency144(Uncoupled_synch_ind,2,:,:));
effi_synch=zeros(4,36,7,32); % 4condi x 36tr x 7freq x 32chan
effi_synco=zeros(4,36,7,32);
effi_synch(1,:,:,:)=effi_ind_synch;
effi_synch(2,:,:,:)=effi_lead_synch;
effi_synch(3,:,:,:)=effi_follow_synch;
effi_synch(4,:,:,:)=effi_mutual_synch;
effi_synco(1,:,:,:)=effi_ind_synco;
effi_synco(2,:,:,:)=effi_lead_synco;
effi_synco(3,:,:,:)=effi_follow_synco;
effi_synco(4,:,:,:)=effi_mutual_synco;
effi_syn=zeros(2,4,36,7,32); % 2syn x 4condi x 36tr x 7freq x 32chan
effi_syn(1,:,:,:,:)=effi_synch;
effi_syn(2,:,:,:,:)=effi_synco;

bwctr_follow_synco  = zeros(36,7,32);
bwctr_lead_synco    = zeros(36,7,32);
bwctr_ind_synco     = zeros(36,7,32);
bwctr_mutual_synco  = zeros(36,7,32);
bwctr_follow_synch  = zeros(36,7,32);
bwctr_lead_synch    = zeros(36,7,32);
bwctr_ind_synch     = zeros(36,7,32);
bwctr_mutual_synch  = zeros(36,7,32);
bwctr_follow_synco(1:18,:,:)  = squeeze(bw_ctr144(R_lead_synco_ind,1,:,:));
bwctr_follow_synco(19:36,:,:) = squeeze(bw_ctr144(L_lead_synco_ind,2,:,:));
bwctr_lead_synco(1:18,:,:)    = squeeze(bw_ctr144(L_lead_synco_ind,1,:,:));
bwctr_lead_synco(19:36,:,:)   = squeeze(bw_ctr144(R_lead_synco_ind,2,:,:));
bwctr_mutual_synco(1:18,:,:)  = squeeze(bw_ctr144(Mutual_synco_ind,1,:,:));
bwctr_mutual_synco(19:36,:,:) = squeeze(bw_ctr144(Mutual_synco_ind,2,:,:));
bwctr_ind_synco(1:18,:,:)     = squeeze(bw_ctr144(Uncoupled_synco_ind,1,:,:));
bwctr_ind_synco(19:36,:,:)    = squeeze(bw_ctr144(Uncoupled_synco_ind,2,:,:));
bwctr_follow_synch(1:18,:,:)  = squeeze(bw_ctr144(R_lead_synch_ind,1,:,:));
bwctr_follow_synch(19:36,:,:) = squeeze(bw_ctr144(L_lead_synch_ind,2,:,:));
bwctr_lead_synch(1:18,:,:)    = squeeze(bw_ctr144(L_lead_synch_ind,1,:,:));
bwctr_lead_synch(19:36,:,:)   = squeeze(bw_ctr144(R_lead_synch_ind,2,:,:));
bwctr_mutual_synch(1:18,:,:)  = squeeze(bw_ctr144(Mutual_synch_ind,1,:,:));
bwctr_mutual_synch(19:36,:,:) = squeeze(bw_ctr144(Mutual_synch_ind,2,:,:));
bwctr_ind_synch(1:18,:,:)     = squeeze(bw_ctr144(Uncoupled_synch_ind,1,:,:));
bwctr_ind_synch(19:36,:,:)    = squeeze(bw_ctr144(Uncoupled_synch_ind,2,:,:));
bwctr_synch=zeros(4,36,7,32); % 4condi x 36tr x 7freq x 32chan
bwctr_synco=zeros(4,36,7,32);
bwctr_synch(1,:,:,:)=bwctr_ind_synch;
bwctr_synch(2,:,:,:)=bwctr_lead_synch;
bwctr_synch(3,:,:,:)=bwctr_follow_synch;
bwctr_synch(4,:,:,:)=bwctr_mutual_synch;
bwctr_synco(1,:,:,:)=bwctr_ind_synco;
bwctr_synco(2,:,:,:)=bwctr_lead_synco;
bwctr_synco(3,:,:,:)=bwctr_follow_synco;
bwctr_synco(4,:,:,:)=bwctr_mutual_synco;
bwctr_syn=zeros(2,4,36,7,32); % 2syn x 4condi x 36tr x 7freq x 32chan
bwctr_syn(1,:,:,:,:)=bwctr_synch;
bwctr_syn(2,:,:,:,:)=bwctr_synco;

%% plot 
% 32 chan x 7 freq in syn and synco (degree ctr, eff, bw)
syn2names={'synch','synco'};
nx3names={'degree ctr','efficiency','betweenness ctr'};
condi4names={'Independent','Leading','Following','Mutual'};

for condi=1:4
    figure;
    for syn = 1:2
        subplot(3,2,syn)
        imagesc(squeeze(mean(squeeze(dgctr_syn(syn,condi,:,:,:)),1)));colorbar;
        ylabel('freq');xlabel('chan');title([syn2names{syn} ' ' nx3names{1} ' mean'])
        clim([0.1 0.6])
        subplot(3,2,2+syn)
        imagesc(squeeze(mean(squeeze(effi_syn(syn,condi,:,:,:)),1)));colorbar;
        ylabel('freq');xlabel('chan');title([syn2names{syn} ' ' nx3names{2} ' mean'])
        clim([0.4 0.8])
        subplot(3,2,4+syn)
        imagesc(squeeze(mean(squeeze(bwctr_syn(syn,condi,:,:,:)),1)));colorbar;
        ylabel('freq');xlabel('chan');title([syn2names{syn} ' ' nx3names{3} ' mean'])
        clim([0.01 0.1])
    end
    sgtitle(condi4names{condi})
end


% barplot of 32 chan each with syn and synco (degree ctr, eff, bw) x 7 freq
% pick mutual condi=4 and betweenness ctr for display
syn2colors=[darkgreen;pink];
condi=4; 
mat7freq=squeeze(mean(squeeze(bwctr_syn(:,condi,:,:,:)),2)); % 2syn x 7 freq x 32 chan
band7labels = {'Delta','Theta', 'Alpha', 'Mu', 'Beta1', 'Beta2','Gamma'};
figure;
clf
for freq=1:7
    disp_mat=squeeze(mat7freq(:,freq,:));
    % Transpose to 32x2 so that each row corresponds to a channel
    data = disp_mat';
    subplot(7,1,freq)
    % Create bar plot
    b = bar(data, 'grouped');
    % Set colors: green for synch, pink for synco
    b(1).FaceColor = syn2colors(1,:);% green
    b(2).FaceColor = syn2colors(2,:);% pink
    % Set x-axis labels
    xticks(1:32);
    xticklabels(arrayfun(@(x) sprintf('Ch%d', x), 1:32, 'UniformOutput', false));
    xlabel('Channel');
    ylabel('Value');
    legend({'synch', 'synco'});
    title(band7labels{freq});
    xtickangle(45); % Optional: angle the labels if too crowded
    ylim([0 0.1])
end
sgtitle([condi4names{condi} ' ' nx3names{3}]);

%% compute (1/abs(Ha-Hb))
open /home/zhibinz2/Documents/GitHub/Motor_coordination_code/Fig3_H_matching.m

% compute 1/abs(Ha-Hb)
% H1over=round(1./abs(H_Lall-H_Rall),1); % 144x1
% H1over=round(1./log(abs(H_Lall-H_Rall)),1); % 144x1
% H1over= min([H_Lall H_Rall],[],2) ./ max([H_Lall H_Rall],[],2); 

H1over=1-abs(H_Lall-H_Rall)./(H_Lall+H_Rall);

H1over=1./abs(H_Lall-H_Rall);

H1over=1./(H_Lall-H_Rall);

H1over=abs(H_Lall-H_Rall); % figure;plot(abs(H_Lall-H_Rall));hold on; plot(1./abs(H_Lall-H_Rall))

H1over=H_Lall-H_Rall;

H1over=[H1over;H1over]; %stacked 288x1

H1over=[H_Lall;H_Rall]; %stacked 288x1

% average value of 1-abs(Hl-hr)/(Hl+HR) in each condition 
condi4_ind; 
syn2_condi4_ind;

H1over;

meanH1over=nan(2,3);
stdH1over=nan(2,3);
semH1over=nan(2,3);
maxH1over=nan(2,3);
minH1over=nan(2,3);
for syn=1:2
    meanH1over(syn,1)=mean(squeeze(H1over(syn2_condi4_ind(1,:,syn))));
    meanH1over(syn,2)=mean(squeeze(H1over(syn2_condi4_ind(2:3,:,syn))),'all');
    meanH1over(syn,3)=mean(squeeze(H1over(syn2_condi4_ind(4,:,syn))));
    stdH1over(syn,1)=std(squeeze(H1over(syn2_condi4_ind(1,:,syn))));
    stdH1over(syn,2)=std([squeeze(H1over(syn2_condi4_ind(2,:,syn))); squeeze(H1over(syn2_condi4_ind(3,:,syn)))]');
    stdH1over(syn,3)=std(squeeze(H1over(syn2_condi4_ind(4,:,syn))));
    semH1over(syn,1)=stdH1over(syn,1)/length(squeeze(H1over(syn2_condi4_ind(1,:,syn))));
    semH1over(syn,2)=stdH1over(syn,2)/length([squeeze(H1over(syn2_condi4_ind(2,:,syn))); squeeze(H1over(syn2_condi4_ind(3,:,syn)))]')
    semH1over(syn,3)=stdH1over(syn,3)/length(squeeze(H1over(syn2_condi4_ind(4,:,syn))));
    maxH1over(syn,1)=max(squeeze(H1over(syn2_condi4_ind(1,:,syn))));
    maxH1over(syn,2)=max([squeeze(H1over(syn2_condi4_ind(2,:,syn))); squeeze(H1over(syn2_condi4_ind(3,:,syn)))]');
    maxH1over(syn,3)=max(squeeze(H1over(syn2_condi4_ind(4,:,syn))));
    minH1over(syn,1)=min(squeeze(H1over(syn2_condi4_ind(1,:,syn))));
    minH1over(syn,2)=min([squeeze(H1over(syn2_condi4_ind(2,:,syn))); squeeze(H1over(syn2_condi4_ind(3,:,syn)))]');
    minH1over(syn,3)=min(squeeze(H1over(syn2_condi4_ind(4,:,syn))));
end
max(maxH1over,[],'all')
min(minH1over,[],'all')

st3names;syn2names;

% Create grouped bar graph
figure;
clf
% Sample data
means = meanH1over; % 2 states x 3 conditions
stds = semH1over;
% Labels
state_labels = {'Synchronization', 'Syncopation'};
cond_labels = {'Independent', 'Unidirectional', 'Bidirectional'};
hold on;
hb = bar(means'); % Transpose to group by condition
% Set colors for each state
colors = [0.2 0.6 0.8; 0.8 0.4 0.4]; % You can change this
for i = 1:numel(hb)
    hb(i).FaceColor = colors(i,:);
end
% Add error bars
ngroups = size(means, 2); % number of conditions
nbars = size(means, 1);   % number of states
% Calculate bar centers for error bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = hb(i).XEndPoints;
end
% Plot error bars
errorbar(x', means', stds', 'k', 'linestyle', 'none', 'LineWidth', 1.2);
% Annotate min and max
for i = 1:nbars
    for j = 1:ngroups
        xpos = x(i, j);
        % Annotate min below bar
        text(xpos, means(i,j)-stds(i,j), sprintf('Min: %.2f', minH1over(i,j)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 9);
        % Annotate max above bar
        text(xpos, means(i,j)+stds(i,j), sprintf('Max: %.2f', maxH1over(i,j)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 9);
    end
end
% Formatting
xticks([1 2 3])
set(gca, 'XTickLabel', cond_labels, 'FontSize', 12);
ylabel('Value with standard deviation');
ylabel('Value with standard error')
legend(state_labels, 'Location', 'northwest');
title('Values of 1-abs(H_L-H_R)./(H_L+H_R)');
title('Values of 1./abs(H_L-H_R)');
title('Values of H_L-H_R');
title('Values of abs(H_L-H_R)');
title('Values of 1./(H_L-H_R)');
box on;
ylim([0.7 1.05])


% open /home/zhibinz2/Documents/GitHub/Motor_coordination_code/Fig3_H_matching.m
H1overcell=cell(2,3);
H1overmean=nan(2,3);
H1overstd=nan(2,3);
H1oversem=nan(2,3);
H1overmax=nan(2,3);
H1overmin=nan(2,3);
% uncouple
H1overcell{1,1}=1-abs(H_Lall(uncoupleInd_LR(synchind(1:18)))      -H_Rall(uncoupleInd_LR(synchind(1:18))))./(H_Lall(uncoupleInd_LR(synchind(1:18)))      +H_Rall(uncoupleInd_LR(synchind(1:18))));

H1overcell{2,1}=1-abs(H_Lall(uncoupleInd_LR(syncoind(1:18)))      -H_Rall(uncoupleInd_LR(syncoind(1:18))))./(H_Lall(uncoupleInd_LR(syncoind(1:18)))      +H_Rall(uncoupleInd_LR(syncoind(1:18))));

% undir
A=1-abs(H_Lall(leadingInd_LR(synchind(1:18)))   -H_Rall(leadingInd_LR(synchind(1:18))))   ./(H_Lall(leadingInd_LR(synchind(1:18)))   +H_Rall(leadingInd_LR(synchind(1:18))));
B=1-abs(H_Rall(followingInd_LR(synchind(1:18))) -H_Lall(followingInd_LR(synchind(1:18)))) ./(H_Rall(followingInd_LR(synchind(1:18))) +H_Lall(followingInd_LR(synchind(1:18))));
H1overcell{1,2}=[A;B];

A=1-abs(H_Lall(leadingInd_LR(syncoind(1:18)))   -H_Rall(leadingInd_LR(syncoind(1:18))))   ./(H_Lall(leadingInd_LR(syncoind(1:18)))   +H_Rall(leadingInd_LR(syncoind(1:18))));
B=1-abs(H_Rall(followingInd_LR(syncoind(1:18))) -H_Lall(followingInd_LR(syncoind(1:18)))) ./(H_Rall(followingInd_LR(syncoind(1:18))) +H_Lall(followingInd_LR(syncoind(1:18))));
H1overcell{2,2}=[A;B];

% bidir
H1overcell{1,3}=1-abs(H_Lall(mutualInd_LR(synchind(1:18)))-H_Rall(mutualInd_LR(synchind(1:18))))./(H_Lall(mutualInd_LR(synchind(1:18)))+H_Rall(mutualInd_LR(synchind(1:18))));
H1overcell{2,3}=1-abs(H_Lall(mutualInd_LR(syncoind(1:18)))-H_Rall(mutualInd_LR(syncoind(1:18))))./(H_Lall(mutualInd_LR(syncoind(1:18)))+H_Rall(mutualInd_LR(syncoind(1:18))));

for syn=1:2
    for st=1:3
        H1overmean(syn,st)=mean(H1overcell{syn,st});
        H1overstd(syn,st)=std(H1overcell{syn,st});
        H1oversem(syn,st)=H1overstd(syn,st)/sqrt(length(H1overcell{syn,st}));
        H1overmax(syn,st)=max(H1overcell{syn,st});
        H1overmin(syn,st)=min(H1overcell{syn,st});
    end
end
max(H1overmax,[],'all')
min(H1overmin,[],'all')

% Create grouped bar graph
figure;
clf
% Sample data
means = H1overmean; % 2 states x 3 conditions
stds = H1oversem;
% Labels
state_labels = {'Synchronization', 'Syncopation'};
cond_labels = {'Independent', 'Unidirectional', 'Bidirectional'};
hold on;
hb = bar(means'); % Transpose to group by condition
% Set colors for each state
colors = [0.2 0.6 0.8; 0.8 0.4 0.4]; % You can change this
for i = 1:numel(hb)
    hb(i).FaceColor = colors(i,:);
end
% Add error bars
ngroups = size(means, 2); % number of conditions
nbars = size(means, 1);   % number of states
% Calculate bar centers for error bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = hb(i).XEndPoints;
end
% Plot error bars
errorbar(x', means', stds', 'k', 'linestyle', 'none', 'LineWidth', 1.2);
% Annotate min and max
for i = 1:nbars
    for j = 1:ngroups
        xpos = x(i, j);
        % Annotate min below bar
        text(xpos, means(i,j)-stds(i,j), sprintf('Min: %.2f', H1overmin(i,j)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 9);
        % Annotate max above bar
        text(xpos, means(i,j)+stds(i,j), sprintf('Max: %.2f', H1overmax(i,j)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 9);
    end
end
% Formatting
xticks([1 2 3])
set(gca, 'XTickLabel', cond_labels, 'FontSize', 12);
ylabel('Value with standard error');
legend(state_labels, 'Location', 'northwest');
title('Values of 1-abs(H_L-H_R)./(H_L+H_R)');
box on;
ylim([0.75 0.95])
%% compute Degree Centrality (1/abs(DCa-DCb)


% compute 1/abs(NXa-NXb)
cd /home/zhibinz2/Documents/GitHub/Motor_coordination_code/
load('nx144.mat'); % 144x2subjectsx7x32
dg_ctr144;
efficiency144;
bw_ctr144;

% dgctr1over  =squeeze(round(1./abs(dg_ctr144(:,1,:,:)     - dg_ctr144(:,2,:,:)    ))); % 144x7x32
% effi1over   =squeeze(round(1./abs(efficiency144(:,1,:,:) - efficiency144(:,2,:,:))));
% bwctr1over  =squeeze(round(1./abs(bw_ctr144(:,1,:,:)     - bw_ctr144(:,2,:,:)    )));

% for i=1:144
%     for freq=1:7
%         for chan=1:32
%             dgctr1over(i,freq,chan)  =min(dg_ctr144(i,1,freq,chan),dg_ctr144(i,2,freq,chan))/max(dg_ctr144(i,1,freq,chan),dg_ctr144(i,2,freq,chan)); % 144x7x32
%             effi1over(i,freq,chan)   =min(efficiency144(i,1,freq,chan),efficiency144(i,2,freq,chan))/max(efficiency144(i,1,freq,chan),efficiency144(i,2,freq,chan));
%             bwctr1over(i,freq,chan)  =min(bw_ctr144(i,1,freq,chan),bw_ctr144(i,2,freq,chan))/max(bw_ctr144(i,1,freq,chan),bw_ctr144(i,2,freq,chan));
%         end
%     end
% end

% stack 2 subjects
dgctr1over  =nan(288,7,32);
dgctr1over(1:144,:,:)=squeeze(dg_ctr144(:,1,:,:));       
dgctr1over(145:288,:,:)=squeeze(dg_ctr144(:,2,:,:));
effi1over   =nan(288,7,32);
effi1over(1:144,:,:) =squeeze(efficiency144(:,1,:,:));   
effi1over(145:288,:,:)=squeeze(efficiency144(:,2,:,:));
bwctr1over  =nan(288,7,32);
bwctr1over(1:144,:,:)=squeeze(bw_ctr144(:,1,:,:));       
bwctr1over(145:288,:,:)=squeeze(bw_ctr144(:,2,:,:));

% combined 3nx into one nx3over for looping
% nx3over=nan(3,144,7,32);
% nx3over(1,:,:,:)=dgctr1over;
% nx3over(2,:,:,:)=effi1over;
% nx3over(3,:,:,:)=bwctr1over;

nx3over=nan(3,288,7,32);
nx3over(1,:,:,:)=dgctr1over;
nx3over(2,:,:,:)=effi1over;
nx3over(3,:,:,:)=bwctr1over;

% % 3 states (combine 2syn)
% condi4_ind;
% independent_ind     =condi4_ind(1,:);
% unidir_ind          =[condi4_ind(2,:) condi4_ind(3,:)];
% bidir_ind           =condi4_ind(4,:);

%% 2 syn 3 states
synch_condi4_ind; synco_condi4_ind;
syn2st3over=cell(2,3); % 2syn x 3 state
syn2st3over{1,1}=synch_condi4_ind(1,:);
syn2st3over{1,2}=[synch_condi4_ind(2,:) synch_condi4_ind(3,:)];
syn2st3over{1,3}=synch_condi4_ind(4,:);
syn2st3over{2,1}=synco_condi4_ind(1,:);
syn2st3over{2,2}=[synco_condi4_ind(2,:) synco_condi4_ind(3,:)];
syn2st3over{2,3}=synco_condi4_ind(4,:);

% stack the indicies
% for syn=1:2
%     for st=1:3
%         syn2st3over{syn,st}
%     end
% end

% correlation betwee 1over (2 syn x 3 states )
coef1over = nan(2,3,3,7,32); % 2syn x 3 state (independent/uni/bi) x 3 nx x 7freq x32chan
for syn=1:2
    for st=1:3
        ind_select=syn2st3over{syn,st};
        % stack the indicies
        ind2_select=[ind_select 144+ind_select];
        H_corr=H1over(ind2_select);
        for nx=1:3
            for freq=1:7
                for chan=1:32
                    nxfreqchan_corr=squeeze(nx3over(nx,ind2_select,freq,chan))';
                    % there are inf in nxfreqchan_corr
                    coef1over(syn,st,nx,freq,chan)=...
                        corr(H_corr,nxfreqchan_corr);
                end
            end
        end
    end
end

% fgure 2 syn x 3 states (each figure with 3nx x (7freq x 32chan))
st3names={'Independent','Unidir','Bidir'};
for syn=1:2
    figure;
    for st=1:3
        for nx=1:3
            subplot(3,3,(nx-1)*3+st) % 3 nx x 3 st
            mat_display=squeeze(coef1over(syn,st,nx,:,:));
            % mat_display(isnan(mat_display))=0; % replace NaN with 0
            imagesc(mat_display);colorbar;colormap('jet');clim([-1 1]);
            xlabel('chan');ylabel('freq');title([nx3names{nx} ' ' st3names{st}])
        end
        % sgtitle([syn2names{syn} ' correlation with 1-abs(H_L-H_R)./(H_L+H_R)'],'Color',syn2colors(syn,:))
        sgtitle([syn2names{syn} ' correlation with H'],'Color',syn2colors(syn,:))
    end
end

%% 2 syn 4 condi
synch_condi4_ind; synco_condi4_ind;
syn2st3over=cell(2,4); % 2syn x 4 condi
for condi=1:4
    syn2st3over{1,condi}=synch_condi4_ind(condi,:)
    syn2st3over{2,condi}=synco_condi4_ind(condi,:)
end

dgctr_syn;
effi_syn; % 2syn x 4condi x 36tr x 7freq x 32chan
bwctr_syn;

nx3over=nan(3,2,4,36,7,32); % 3nx x 2syn x 4condi x 36tr x 7freq x 32chan
nx3over(1,:,:,:,:,:)=dgctr_syn;
nx3over(2,:,:,:,:,:)=effi_syn;
nx3over(3,:,:,:,:,:)=bwctr_syn;
% just plot mean of nx
for syn=1:2
    figure;
    for condi=1:4
        for nx=1:3
            subplot(3,4,(nx-1)*4+condi) % 3 nx x 4 condi
            mat_tmp=squeeze(nx3over(nx,syn,condi,:,:,:));
            mat_display=squeeze(mean(mat_tmp,1));
            % mat_display(isnan(mat_display))=0; % replace NaN with 0
            imagesc(mat_display);colorbar;colormap('jet');
            if nx==1
                clim([0.2 0.6]);
            elseif nx==2
                clim([0.4 0.8]);
            else
                clim([0.01 0.1]);
            end
            xlabel('chan');ylabel('freq');title([nx3names{nx} ' ' condi4names{condi}])
        end
        sgtitle([syn2names{syn} ' mean'],'Color',syn2colors(syn,:))
    end
end
% make independent the reference condition and contrast each value against independent
nx3contrast=nan(3,2,3,7,32); % 3nx x 2syn x 3 contrast condi x 36tr x 7freq x 32chan
for syn=1:2
    figure;
    for nx=1:3
        for condi=2:4
            mat_indi_tmp=squeeze(nx3over(nx,syn,1,:,:,:));
            mat_indi=squeeze(mean(mat_indi_tmp,1));
            mat_tmp=squeeze(nx3over(nx,syn,condi,:,:,:));
            mat_display=squeeze(mean(mat_tmp,1))-mat_indi;
            nx3contrast(nx,syn,condi-1,:,:)=mat_display;
            subplot(3,3,(nx-1)*3+condi-1);
            imagesc(mat_display);colorbar;colormap('jet');
            xlabel('chan');ylabel('freq');title([nx3names{nx} ' ' condi4names{condi}]);
            if nx==1
                clim([-0.15 0.15])
            elseif nx==2
                clim([-0.15 0.15])
            else
                clim([-0.06 0.06])
            end
        end
    end
    sgtitle([syn2names{syn} ' constrast mean - independent'],'Color',syn2colors(syn,:))
end


% correlation betwee 1over (2 syn x 3 states )
coef1over = nan(2,4,3,7,32); % 2syn x 4condi (independent/leading/following/bi) x 3 nx x 7freq x32chan
for syn=1:2
    for condi=1:4
        ind_select=syn2st3over{syn,condi};
        % stack the indicies
        ind2_select=[ind_select 144+ind_select];
        H_corr=H1over(ind2_select);
        for nx=1:3
            for freq=1:7
                for chan=1:32
                    nxfreqchan_corr=squeeze(nx3over(nx,syn,condi,:,freq,chan));
                    % there are inf in nxfreqchan_corr
                    coef1over(syn,condi,nx,freq,chan)=...
                        corr(H_corr,nxfreqchan_corr);
                end
            end
        end
    end
end

% fgure 2 syn x 4 condi (each figure with 3nx x (7freq x 32chan))
for syn=1:2
    figure;
    for condi=1:4
        for nx=1:3
            subplot(3,4,(nx-1)*4+condi) % 3 nx x 4 condi
            mat_display=squeeze(coef1over(syn,condi,nx,:,:));
            % mat_display(isnan(mat_display))=0; % replace NaN with 0
            imagesc(mat_display);colorbar;colormap('jet');clim([-1 1]);
            xlabel('chan');ylabel('freq');title([nx3names{nx} ' ' condi4names{condi}])
        end
        % sgtitle([syn2names{syn} ' correlation with 1-abs(H_L-H_R)./(H_L+H_R)'],'Color',syn2colors(syn,:))
        % sgtitle([syn2names{syn} ' correlation with H'],'Color',syn2colors(syn,:))
        sgtitle([syn2names{syn} ' correlation with abs(H_L-H_R)'],'Color',syn2colors(syn,:))
        sgtitle([syn2names{syn} ' correlation with H_L-H_R'],'Color',syn2colors(syn,:))
    end
end


%% zscore
% size(r100_4crH1over)
load('zscore_crH1over.mat')
size(zscore_crH)

addpath(genpath('/home/zhibinz2/Documents/GitHub/eeglab'))
eeglab

cd /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/channels_info
load('chaninfo.mat')

for ccc=1:3
    figure('Position',[1084           1        1196        1413]);  
    % get(gcf, 'outerposition');figure('Position',[ans]); 
    
    % Define subplot grid size
    nRows = 8; % 2 syn × 4 condi
    nCols = 7; % 7 freqs

    clf
    for syn=1:2
        for condi=1:4
            for freq=1:7
                subplot(8,7,(syn-1)*4*7+7*(condi-1)+freq)
                array_display=squeeze(zscore_crH(ccc,syn,condi,freq,:));
                topoplot(array_display,chaninfo','nosedir','+X'); %colorbar;
                clim([-3 3])
                % title([condi4names{condi}]);
                % subtitle([syn2names{syn} ' freq: ' band7labels{freq}],'Color',syn2colors(syn,:))
            end
        end
    end
    % sgtitle([nx3names{ccc}])
    % display(num2str(ccc))

    % Add row (left side) labels
    for syn = 1:2
        for condi = 1:4
            row = (syn-1)*4 + condi;
            % normalized position for text
            y_pos = 1-0.06 - 0.85*(row - 0.5)/nRows;
            annotation('textbox', [0.01, y_pos, 0.1, 0.03], ...
                'String', {[syn2names{syn}], [condi4names{condi}]}, ...
                'EdgeColor', 'none', ...
                'HorizontalAlignment', 'left', ...
                'FontWeight', 'bold', ...
                'Color', syn2colors(syn,:));
        end
    end

    % delete(findall(gcf, 'Type', 'annotation'))
    
    % Add column (top) labels for frequencies
    for freq = 1:7
        x_pos = 0.07+0.8*(freq - 0.5) / nCols;
        annotation('textbox', [x_pos, 0.93, 0.1, 0.03], ...
            'String', band7labels{freq}, ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center', ...
            'FontWeight', 'bold');
    end

    % delete(findall(gcf, 'Type', 'annotation'))
    annotation('textbox',[0.3, 0.97, 0.4, 0.03], ...
        'String',{[nx3names{ccc}],'Zscore of correlation with 1-abs(H_L-H_R)./(H_L+H_R)'}, ...
        'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', ...
        'FontWeight', 'bold', ...
        'FontSize',12);

    % Create an invisible axes to attach the colorbar
    cb_ax = axes('Position', [0.2 0.05 0.6 0.02], 'Visible', 'off');
    % Create colorbar on that axes
    cbar = colorbar(cb_ax, 'Location', 'southoutside');
    % Adjust colormap and limits (optional)
    colormap(jet)
    caxis([-3 3])  % Match your topoplot clim

end

% reorganzie into matrix
zscore1over = nan(2,4,3,7,32); % 2syn x 4condi (independent/leading/following/bi) x 3 nx x 7freq x32chan
for ccc=1:3
    for syn=1:2
        for condi=1:4
            for freq=1:7
                zscore1over(syn,condi,ccc,freq,:)=squeeze(zscore_crH(ccc,syn,condi,freq,:));
            end
        end
    end
end
for syn=1:2
    figure;
    for condi=1:4
        for nx=1:3
            subplot(3,4,(nx-1)*4+condi) % 3 nx x 4 condi
            mat_display=squeeze(zscore1over(syn,condi,nx,:,:));
            % mat_display(isnan(mat_display))=0; % replace NaN with 0
            imagesc(mat_display);colorbar;colormap('jet');clim([-3 3]);
            xlabel('chan');ylabel('freq');title([nx3names{nx} ' ' condi4names{condi}])
        end
        sgtitle([syn2names{syn} ' zscore correlation with 1-abs(H_L-H_R)./(H_L+H_R)'],'Color',syn2colors(syn,:))
        % sgtitle([syn2names{syn} ' correlation with H'],'Color',syn2colors(syn,:))
    end
end


%% Correlation with Xcorr
open /home/zhibinz2/Documents/GitHub/Motor_coordination_code/Fig4_Cross_Corr.m
% sort order for EEG in the same order for Xcorr10Lag (maybe not)
% reorganzie (degree ctr, eff, bw)
cd /home/zhibinz2/Documents/GitHub/Motor_coordination_code/
load('nx144.mat'); % 144x2x7x32 in time sequeence

data_path='/home/zhibinz2/Documents/GitHub/Cleaned_data/';
Xcorr144Lag=nan(144,21);
n=1;
for r=1:numSes;
    clear intervals conditions sortorders
    runid = num2str(seeds(r,:));
    load([data_path  'clean_' runid '.mat'],'conditions','intervals');
    for j = 1:12
        clear intL_good_dmean intR_good_dmean
        % remove the mean
        intL_good_dmean=intervals{j}(:,1)-mean(intervals{j}(:,1));
        intR_good_dmean=intervals{j}(:,2)-mean(intervals{j}(:,2));
        % Xcorr based on int_dmean_drm (before d removal)********
        r12=[];lags12=[];
        [r12,lags12]=xcorr(intL_good_dmean,intR_good_dmean,10,'normalized');
        Xcorr144Lag(n,:)=r12;
        n=n+1;
    end
end

syn2Ind={[1:2:11],[2:2:12]};syn2names={'Synch','Synco'};

% extract lag -1 0 +1 
Xcorr3lg=nan(144,3);
for n=1:144
    Xcorr3lg(n,:)=Xcorr144Lag(n,10:12);
end

% organize Xcorr3lg from 144 tr x 3 to 2 syn x 4 condi x 36 tr x 3? (maybe no)
meanXcorr3lg=nan(2,3,3);
stdXcorr3lg=nan(2,3,3);
semXcorr3lg=nan(2,3,3);
maxXcorr3lg=nan(2,3,3);
minXcorr3lg=nan(2,3,3);
for syn=1:2
    meanXcorr3lg(syn,1,:)=mean(squeeze(Xcorr3lg(syn2_condi4_ind(1,:,syn),:)));
    meanXcorr3lg(syn,2,:)=mean(squeeze(Xcorr3lg([syn2_condi4_ind(2,:,syn) syn2_condi4_ind(3,:,syn)],:)),1);
    meanXcorr3lg(syn,3,:)=mean(squeeze(Xcorr3lg(syn2_condi4_ind(4,:,syn),:)));
    stdXcorr3lg(syn,1,:)=std(squeeze(Xcorr3lg(syn2_condi4_ind(1,:,syn),:)));
    stdXcorr3lg(syn,2,:)=std([squeeze(Xcorr3lg(syn2_condi4_ind(2,:,syn),:)); squeeze(Xcorr3lg(syn2_condi4_ind(3,:,syn),:))]);
    stdXcorr3lg(syn,3,:)=std(squeeze(Xcorr3lg(syn2_condi4_ind(4,:,syn),:)));
    semXcorr3lg(syn,1,:)=stdXcorr3lg(syn,1)/length(squeeze(Xcorr3lg(syn2_condi4_ind(1,:,syn),:)));
    semXcorr3lg(syn,2,:)=stdXcorr3lg(syn,2)/length([squeeze(Xcorr3lg(syn2_condi4_ind(2,:,syn),:)); squeeze(Xcorr3lg(syn2_condi4_ind(3,:,syn),:))]);
    semXcorr3lg(syn,3,:)=stdXcorr3lg(syn,3)/length(squeeze(Xcorr3lg(syn2_condi4_ind(4,:,syn),:)));
    maxXcorr3lg(syn,1,:)=max(squeeze(Xcorr3lg(syn2_condi4_ind(1,:,syn),:)));
    maxXcorr3lg(syn,2,:)=max([squeeze(Xcorr3lg(syn2_condi4_ind(2,:,syn),:)); squeeze(Xcorr3lg(syn2_condi4_ind(3,:,syn),:))]);
    maxXcorr3lg(syn,3,:)=max(squeeze(Xcorr3lg(syn2_condi4_ind(4,:,syn),:)));
    minXcorr3lg(syn,1,:)=min(squeeze(Xcorr3lg(syn2_condi4_ind(1,:,syn),:)));
    minXcorr3lg(syn,2,:)=min([squeeze(Xcorr3lg(syn2_condi4_ind(2,:,syn),:)); squeeze(Xcorr3lg(syn2_condi4_ind(3,:,syn),:))]);
    minXcorr3lg(syn,3,:)=min(squeeze(Xcorr3lg(syn2_condi4_ind(4,:,syn),:)));
end

st3names;syn2names;

% Create grouped bar graph
lag3names={'-1','0','+1'};
figure;
clf
% Sample data
for lg=1:3;
    subplot(1,3,lg)
    means = meanXcorr3lg(:,:,lg); % 2 states x 3 conditions
    stds = semXcorr3lg(:,:,lg);
    % Labels
    state_labels = {'Synchronization', 'Syncopation'};
    cond_labels = {'Independent', 'Unidirectional', 'Bidirectional'};
    hold on;
    hb = bar(means'); % Transpose to group by condition
    % Set colors for each state
    colors = [0.2 0.6 0.8; 0.8 0.4 0.4]; % You can change this
    for i = 1:numel(hb)
        hb(i).FaceColor = colors(i,:);
    end
    % Add error bars
    ngroups = size(means, 2); % number of conditions
    nbars = size(means, 1);   % number of states
    % Calculate bar centers for error bars
    x = nan(nbars, ngroups);
    for i = 1:nbars
        x(i,:) = hb(i).XEndPoints;
    end
    % Plot error bars
    errorbar(x', means', stds', 'k', 'linestyle', 'none', 'LineWidth', 1.2);
    % Annotate min and max
    for i = 1:nbars
        for j = 1:ngroups
            xpos = x(i, j);
            % Annotate min below bar
            text(xpos, means(i,j)-stds(i,j), sprintf('Min: %.2f', minXcorr3lg(i,j,lg)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 9);
            % Annotate max above bar
            text(xpos, means(i,j)+stds(i,j), sprintf('Max: %.2f', maxXcorr3lg(i,j,lg)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 9);
        end
    end
    % Formatting
    xticks([1 2 3])
    set(gca, 'XTickLabel', cond_labels, 'FontSize', 12);
    ylabel('Value with standard deviation');
    ylabel('Value with standard error')
    legend(state_labels, 'Location', 'northwest');
    title(['Xcorr at lag ' lag3names{lg}]);
    box on;
    ylim([0 0.7])
end

%% Permutation with contrast to independent
open /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/hyperscanEEG_correlation/permutation20250425.ipynb
cd /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/hyperscanEEG_correlation/
load('nx3permu.mat') %  3nx x 2syn x 3contrast condi x7freq x32chan
load('nx3permu2.mat') %  3nx x 2syn x 3contrast condi x7freq x32chan
nx3gain10   
nx3gain01   
nx3counter10
nx3counter01

for syn=1:2
    figure('Position',[1327         772         884         510]); % get(gcf, 'outerposition') 
    % get(gcf, 'outerposition'); figure('Position',[ans]);  ;
    for nx=1:3
        for condi=2:4
            tmp=squeeze(nx3gain10(nx,syn,condi-1,:,:));
            ind_display=(tmp>97500);
            tmp(ind_display)=1;
            mat_display=tmp;
            subplot(3,3,(nx-1)*3+condi-1);
            imagesc(mat_display);colorbar;colormap('jet');
            xlabel('chan');ylabel('freq');title([nx3names{nx} ' ' condi4names{condi}]);
            if nx==1
                clim([1-0.6 1.6])
            elseif nx==2
                clim([1-0.3 1.3])
            else
                clim([1-5 1+5])
            end
        end
    end
    sgtitle([syn2names{syn} ' permutation (> independent)'],'Color',syn2colors(syn,:))
end
