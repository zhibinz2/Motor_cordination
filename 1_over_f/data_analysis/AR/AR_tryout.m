%% AR toolboxestryout
cd /home/zhibin/Documents/GitHub/matlab/external/arfit
open arfit
open arconf
open ardem

cd /home/zhibin/Documents/GitHub/matlab/external/MVGC/core
open tsdata_to_var.m % A
open tsdata_to_autocov.m % G

cd /home/zhibin/Documents/GitHub/matlab/external/MSE-VARFI
open arfi2ar2.m
open analysis_script.m

%% select the subject data
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);
sessions={'synch','synco','synch','synco','synch','synco','synch','synco','synch','synco',...
        'synch','synco'};
cd /ssd/zhibin/1overf/
s=1
runid=num2str(seeds(s,:));
load(['/ssd/zhibin/1overf/' runid '_2P/Cleaned_data/clean_' runid '.mat' ]) 


cd /ssd/zhibin/1overf/all_session20220713_1005
load('int_dmean_drm.mat')

%% try out reorganizing methods

% Method1a: concatenate into one big trial
Int12LR=[];
for b=1:size(test_data,2)
    for s=1:size(test_data,3)
            Int12LR=[Int12LR; [test_data{1,b,s} test_data{2,b,s}]];
    end
end
% Method1b: randomized the concatenation
Int12LR=[];
for b=randperm(size(test_data,2))
    for s=randperm(size(test_data,3))
            Int12LR=[Int12LR; [test_data{1,b,s} test_data{2,b,s}]];
    end
end


% Method2: truncate to 100 samples
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing
% [Int12LR] = reorg4arfit(test_data)
Int12LR=nan(100,2,size(test_data,2)*size(test_data,3));
nL=1;nR=1;
for p=1:size(test_data,1)
    for b=1:size(test_data,2)
        for s=1:size(test_data,3)
            if p==1;
                Int12LR(1:100,1,nL)=test_data{1,b,s}(1:100,:);
                nL=nL+1;
            end
            if p==2
                Int12LR(1:100,2,nR)=test_data{2,b,s}(1:100,:);
                nR=nR+1;
            end
        end
    end
end
clc
Int12LR
% Method3: resample to 200 samples
Maxlength=200;
Int12LR=nan(Maxlength,2,size(test_data,2)*size(test_data,3));
nL=1;nR=1;
for p=1:size(test_data,1)
    for b=1:size(test_data,2)
        for s=1:size(test_data,3)
            if p==1;
                Int12LR(1:Maxlength,1,nL)=resample(test_data{1,b,s},Maxlength,length(test_data{1,b,s}));
                nL=nL+1;
            end
            if p==2
                Int12LR(1:Maxlength,2,nR)=resample(test_data{2,b,s},Maxlength,length(test_data{2,b,s}));
                nR=nR+1;
            end
        end
    end
end
clc
Int12LR
%% try out arfit and arconf
cd /home/zhibin/Documents/GitHub/matlab/external/arfit
% demostration example
open ardem

% % one tapping int series from L
% Int1L=[int_dmean_drm{1,1,1}];
% % 1 int series from L & R
% Int1LR=[int_dmean_drm{1:2,1,1}];
% % with 3rd dimension as trials
% Int12LR=nan(100,2,12);
% for p=1:2
%     for b=1:12
%         for s=1%:12
%             Int12LR(1:100,p,b)=int_dmean_drm{p,b,s}(1:100,:);
%         end
%     end
% end


% try arfit
% select v to try out
% v=Int1L; 
% v=Int12LR(1:100,1,12);
v_arfit=Int12LR;
% v_arfit=Int12LR+1500;

% v; % time series data x observations x trials

addpath cd /home/zhibin/Documents/GitHub/matlab/external/arfit
pmin=1; % minimal order
pmax=10; % maximal order
[w, A, C, sbc, fpe, th]=arfit(v_arfit, pmin, pmax) % , 'sbc', 'zero')
% v(k,:)' = w' + A1*v(k-1,:)' +...+ Ap*v(k-p,:)' + noise(C)

%  For a bivariate AR(2process)
%  v(k,:)' = w' + A1*v(k-1,:)' + A2*v(k-2,:)' + eta(k,:)',
w; % least squares estimates of the intercept vector w (m x 1)
A; % coefficient matrices A=[A1 ... Ap]
C; % noise covariance matrix (m x m)
sbc; % Schwarz's Bayesian Criterion
fpe; % logarithm of Akaike's Final Prediction Error
th; % contains information needed for the computation of confidence intervals

%  Next it is necessary to check whether the fitted model is adequate
%  to represent the time series v. A necessary condition for model
%  adequacy is the uncorrelatedness of the residuals. The module ARRES
[siglev,res] = arres(w,A,v_arfit);
%  returns the time series of residuals res(:,:,itr) for each
%  realization itr=1,...,ntr of the original time series and the
%  significance level siglev with which a modified Li-McLeod
%  portmanteau test rejects the null hypothesis that the residuals
%  (using all available realizations) are uncorrelated. A model
%  passes this test if, say, siglev > 0.05.  In our example, the
%  significance level of the modified Li-McLeod portmanteau
%  statistic is 
disp(siglev); 

% locate p for the minimal value
[minSBC,p_indSBC]=min(sbc)
[minFPE,p_indFPE]=min(fpe)

% try arconf
[Aerr, werr]=arconf(A, C, w, th)
% (A +/- Aerr) and (w +/- werr) are approximate 95% confidence intervals for the elements of the coefficient matrix A
Aerr; % margins of error
werr; % and for the components of the intercept vector w. The input

% See if the intercept cross zero
figure;errorbar(1:length(w),w,werr,'ko');
xlim([0 length(w)+1]);yline(0,'color',[1 0.8 0.2]);

% Granger Causality
% https://www.youtube.com/watch?v=XqsSB_vpHLs
% GC of L->R
GC_L2R=log(C(2,2)/C(1,2))
% GC of R->L
GC_R2L=log(C(1,1)/C(1,2))
% plot
figure;bar([GC_L2R GC_R2L]);xticks([1 2]);xticklabels({'L->R','R->L'});

%% try out tsdata_to_var and tsdata_to_autocov (MVGC)
v_tsdata=permute(Int12LR,[2,1,3]);
X=v_tsdata; % multi-trial time series data
p=10; % model order (number of lags)
regmode='LWR' ; % regression mode: 'LWR' (default) or 'OLS'
[A,SIG,E] = tsdata_to_var(X,p,regmode);
A; % VAR coefficients matrix
SIG; % residuals covariance matrix
E; % residuals time series

cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/VAR_Granger
open mvgc_demo_autocov_tryout.m


X=v_tsdata; % multi-trial time series data
q=5; % number of lags
G = tsdata_to_autocov(X,q)
G; % sample autocovariance sequence

% Parameters
ntrials   = size(X,3);     % number of trials
nobs      = size(X,2);  % number of observations per trial

regmode   = 'LWR';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
momax     = 20;     % maximum model order for model order estimation

% recomended acmaxlags minimum = 105
acmaxlags = 20;   % maximum autocovariance lags (empty for automatic calculation)

tstat     = '';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')

% for plotting only
fs        = 1.3;    % sample rate (Hz)
fres      = [];     % frequency resolution (empty for automatic calculation)

nvars = size(X,1); % number of variables

figure;bar([F(2,1) F(1,2)]);;xticks([1 2]);xticklabels({'L->R','R->L'});
%% Issues
% each trial has to be the same length for the input as a 3d matrix (1:100)

% arfit: A size varies with number of trials

% arfit: A=[] 100x2x12. 
% arconf: Index in position 1 exceeds array bounds.

% arconf: intercef+error cross zero bad?

% The A output (coefficient matrices) of arfit variea


%%
