clear
%% addpath
addpath /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing
addpath /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/DFA

% d estimate and removal
% http://www.lucafaes.net/LMSE-MSE_VARFI.html
addpath /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/MSE-VARFI
% Granger Causality
addpath(genpath('/home/zhibinz2/Documents/GitHub/MVGC1'));
cd /home/zhibinz2/Documents/GitHub/MVGC1
run startup
addpath /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/VAR_Granger


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

%%
% SECT 9 Xcorr before/after d removal and save new intervals (matched int) (sorted order) 
% Trials={'uncouple','uncouple','uncouple','L-lead','L-lead','L-lead',...
%     'R-lead','R-lead','R-lead','mutual','mutual','mutual'}';
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);
% sessions={'synch','synco','synch','synco','synch','synco','synch','synco','synch','synco',...
%     'synch','synco'};

% (Before d removal)********
% d estimation from H (DFA method) based on the good matched intervals
H=zeros(2,12,numSes);XcorrPeakLag=zeros(12,numSes);XcorrPeak=zeros(12,numSes);
Xcorr10Lag=nan(numSes,12,21);
cd /ssd/zhibin/1overf/
tic
for r=1:numSes
    clear intervals conditions sortorders
    runid = num2str(seeds(r,:));
    path = [runid '_2P/Cleaned_data/'];
    load([path  'clean_' runid '.mat'],'intervals','conditions');
    % sort order and plot results
    [x,sortorder]=sort(conditions);
    for j = 1:12
        clear intL_good_dmean intR_good_dmean
        % After aligment
        % intervals{sortorder(j)}(:,1);
        % intervals{sortorder(j)}(:,2);
        % remove the mean
        intL_good_dmean=intervals{sortorder(j)}(:,1)-mean(intervals{sortorder(j)}(:,1));
        intR_good_dmean=intervals{sortorder(j)}(:,2)-mean(intervals{sortorder(j)}(:,2));
        % estimate the d
        [~,H(1,j,r)]=DFA_main(intL_good_dmean);
        [~,H(2,j,r)]=DFA_main(intR_good_dmean);
%         % Xcorr based on int_dmean_drm (before d removal)********
%         r12=[];lags12=[];
%         [r12,lags12]=xcorr(intL_good_dmean,intR_good_dmean,10,'normalized');
%         XcorrPeakLag(j,r)=lags12(find(r12==max(r12)));
%         XcorrPeak(j,r)=max(r12);
%         Xcorr10Lag(r,j,1:21)=r12;
    end
end
toc

d1=H-0.5;

%%

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
% d removal
H_est2=zeros(2,12,numSes);
int_dmean_drm=cell(2,12,numSes);
cd /ssd/zhibin/1overf/
tic
for r=1:numSes;
    clear intervals conditions sortorders
    runid = num2str(seeds(r,:));
    path = [runid '_2P/Cleaned_data/'];
    load([path  'clean_' runid '.mat'],'intervals','conditions');
    % sort order and plot results
    [x,sortorder]=sort(conditions);
    for j = 1:12
        clear intL_good_dmean intR_good_dmean
        % After aligment
        % intervals{sortorder(j)}(:,1);
        % intervals{sortorder(j)}(:,2);
        % remove the mean
        intL_good_dmean=intervals{sortorder(j)}(:,1)-mean(intervals{sortorder(j)}(:,1));
        intR_good_dmean=intervals{sortorder(j)}(:,2)-mean(intervals{sortorder(j)}(:,2));
        % d removal  (method1: based on DFA, d=h-0.5)
        [int_dmean_drm{1,j,r}]=remove_d(intL_good_dmean,d1(1,j,r));
        [int_dmean_drm{2,j,r}]=remove_d(intR_good_dmean,d1(2,j,r));
%         % H estimate again using DFA method
%         [~,H_est2(1,j,r)]=DFA_main(int_dmean_drm{1,j,r});
%         [~,H_est2(2,j,r)]=DFA_main(int_dmean_drm{2,j,r});
    end
end
toc


%%
% orgainized the synchronization trials into 4 conditions
condi4Ind={[1:3],[4:6],[7:9],[10:12]};
condi4names={'Uncouple','L-lead','R-lead','Mutual'};
syn2Ind={[1:2:11],[2:2:12]};syn2names={'Synch','Synco'};

% MVGC
% organize into time x 2 L/R matrix then apply MVGC and plot
Fs=cell(2,4);
for syn=1:2
    for condi=1:4
        test_data=int_dmean_drm(:,condi4Ind{condi},syn2Ind{syn});
        % Method1a: concatenate into one big trial
        Int12LR=[];
        for b=1:size(test_data,2) % block
            for s=1:size(test_data,3) % session
                    Int12LR=[Int12LR; [test_data{1,b,s} test_data{2,b,s}]];
            end
        end
%         % Method1b: randomized the concatenation
%         Int12LR=[];
%         for b=randperm(size(test_data,2))
%             for s=randperm(size(test_data,3))
%                     Int12LR=[Int12LR; [test_data{1,b,s} test_data{2,b,s}]];
%             end
%         end
%         if syn==1; subplot(2,4,condi);end
%         if syn==2; subplot(2,4,4+condi);end
        X=[];
        X=permute(Int12LR,[2,1]);
        % Parameters
        ntrials   = size(X,3);     % number of trials
        nobs      = size(X,2);  % number of observations per trial
        regmode   = 'LWR';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
        icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default) 
        morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
        momax     = 10;     % maximum model order for model order estimation
        % recomended acmaxlags minimum = 105
        acmaxlags = 10;   % maximum autocovariance lags (empty for automatic calculation)
        tstat     = '';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
        alpha     = 0.05;   % significance level for significance test
        mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')
        % for plotting only
        fs        = 1.3;    % sample rate (Hz)
        fres      = [];     % frequency resolution (empty for automatic calculation)
        nvars = size(X,1); % number of variables
        % model order estimation
        [AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X,momax,icregmode);
        % Select model order.
        if     strcmpi(morder,'actual')
            morder = amo;
            fprintf('\nusing actual model order = %d\n',morder);
        elseif strcmpi(morder,'AIC')
            morder = moAIC;
            fprintf('\nusing AIC best model order = %d\n',morder);
        elseif strcmpi(morder,'BIC')
            morder = moBIC;
            fprintf('\nusing BIC best model order = %d\n',morder);
        else
            fprintf('\nusing specified model order = %d\n',morder);
        end
        % VAR model estimation 
        [A,SIG] = tsdata_to_var(X,morder,regmode);
        % Check for failed regression
        assert(~isbad(A),'VAR estimation failed');
        % Autocovariance calculation
        [G,info] = var_to_autocov(A,SIG,acmaxlags);
        var_acinfo(info,true); % report results (and bail out on error)
        % Granger causality calculation: time domain 
        F = autocov_to_pwcgc(G);
        % save F
        % F(2,1): L->R; F(1,2): R->L;
        Fs{syn,condi}=[F(2,1) F(1,2)];
    end
end

%%
% organize into time x 2 L/R matrix then apply MVGC and collect p values (2 suplots separate synch/o)
Fs=nan(2,4);Ps=nan(2,4); % 2 syn type x 4 conditions
if true % calculation of GC and P values
for condi=1
    for t=1:2
    test_data=[];
    test_data=int_dmean_drm(:,condi4Ind{condi},syn2Ind{t});
    L_mat=[];R_mat=[];
    L_mat=cell2mat(reshape(squeeze(test_data(1,:,:)),18,1));
    R_mat=cell2mat(reshape(squeeze(test_data(2,:,:)),18,1));
    X=[];
    X=[L_mat R_mat; R_mat L_mat];
    X=permute(X,[2,1]);
    [F,pval,sig] = myGCautocov(X);
    Fs(t,1)=F(2,1);
    Ps(t,1)=pval(2,1);
    end
end

for condi=[2 3]
    for t=1:2
    test_data=[];
    test_data=int_dmean_drm(:,condi4Ind{2},syn2Ind{t});
    L_mat=[];R_mat=[];
    L_mat=cell2mat(reshape(squeeze(test_data(1,:,:)),18,1));
    R_mat=cell2mat(reshape(squeeze(test_data(2,:,:)),18,1));
    X1=[];
    X1=[L_mat R_mat];

    test_data=[];
    test_data=int_dmean_drm(:,condi4Ind{3},syn2Ind{t});
    L_mat=[];R_mat=[];
    L_mat=cell2mat(reshape(squeeze(test_data(1,:,:)),18,1));
    R_mat=cell2mat(reshape(squeeze(test_data(2,:,:)),18,1));
    X2=[];
    X2=[R_mat L_mat];
    
    X=[];
    X=[X1; X2];
    X=permute(X,[2,1]);
    [F,pval,sig] = myGCautocov(X);
    Fs(t,2)=F(2,1);Fs(t,3)=F(1,2);
    Ps(t,2)=pval(2,1);Ps(t,3)=pval(1,2);
    end
end

for condi=4
    for t=1:2
    test_data=[];
    test_data=int_dmean_drm(:,condi4Ind{condi},syn2Ind{t});
    L_mat=[];R_mat=[];
    L_mat=cell2mat(reshape(squeeze(test_data(1,:,:)),18,1));
    R_mat=cell2mat(reshape(squeeze(test_data(2,:,:)),18,1));
    X=[];
    X=[L_mat R_mat; R_mat L_mat];
    X=permute(X,[2,1]);
    [F,pval,sig] = myGCautocov(X);
    Fs(t,4)=F(2,1);
    Ps(t,4)=pval(2,1);
    end
end
end

%%

if true % barplot 4groups x 2types
canvas(0.23, 0.4);
model_series = Fs';
clear b
b = bar(model_series, 'FaceColor','flat');
b(1).FaceColor=darkgreen;b(2).FaceColor=pink;
xticks(1:4);xticklabels({'Uncoupled','Leader -> Follower','Follower -> Leader','Mutual'});
xl = get(gca,'XTickLabel');  
set(gca,'XTickLabel',xl,'fontsize',20,'FontWeight','bold')
% set(gca,'XTickLabel',{'Uncoupled','Leading','Following','Mutual'},'fontsize',12,'FontWeight','bold');
% label p values
% delete(findall(gcf,'type','text'))
% for t=1:2
%     for p=1:4
%         if Ps(t,p) <0.01; signifColor=red; else; signifColor=black; end
%         if t==1
%             text(p-0.27, Fs(t,p)+0.005, sprintf('*','Color',signifColor,'FontSize',16))
%             % text(p-0.27, Fs(t,p)+0.005, sprintf('p=%.2f',Ps(t,p)),'Color',signifColor,'FontSize',16);
%         else
%             text(p+0.01, Fs(t,p)+0.005, sprintf('*','Color',signifColor,'FontSize',16))
%             % text(p+0.01, Fs(t,p)+0.005, sprintf('p=%.2f',Ps(t,p)),'Color',signifColor,'FontSize',16);
%         end
%     end
% end
legend({'Synch','Synco'},'location','northwest','FontSize', 30);
ylabel('GC','FontSize',17); 
yl = get(gca,'YTickLabel');  
set(gca,'YTickLabel',yl,'fontsize',20,'FontWeight','bold');
% delete(sg)
% sg=annotation('textbox',[0.2 0.01 0.7 0.07],'string',...
%     {['MVGC (concatenate data) H-int ^{* SECT 21}'],char(datetime('now'))});
set(gcf,'color','w'); % set background white for copying in ubuntu
end



%%
cd /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/Publish_figures
