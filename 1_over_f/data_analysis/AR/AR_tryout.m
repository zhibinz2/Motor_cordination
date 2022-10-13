%% AR tryout

cd /home/zhibin/Documents/GitHub/matlab/external/arfit
open arfit
open arconf


cd /home/zhibin/Documents/GitHub/matlab/external/MVGC/core
open tsdata_to_var.m % A
open tsdata_to_autocov.m % G


cd /home/zhibin/Documents/GitHub/matlab/external/MSE-VARFI
open arfi2ar2.m

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

%% try out arfit and arconf
cd /home/zhibin/Documents/GitHub/matlab/external/arfit
open ardem


% one tapping int series from L
Int1L=[int_dmean_drm{1,1,1}];

% 1 int series from L & R
Int1LR=[int_dmean_drm{1:2,1,1}];

% add 3rd dimension as trials
% replicate 12 trials
% Int12LR=repmat(Int1LR,1,1,12);
Int12LR=nan(100,2,12);
for p=1:2
    for b=1:12
        for s=1%:12
            Int12LR(1:100,p,b)=int_dmean_drm{p,b,s}(1:100,:);
%             Int12LR(1:100,3,b)=int_dmean_drm{2,b,s}(1:100,:);
%             Int12LR(1:100,4,b)=int_dmean_drm{2,b,s}(1:100,:);
        end
    end
end


% try arfit
% select v to try out
% v=Int1L; 
% v=Int12LR(1:100,1,12);
v=Int12LR;

v; % time series data x observations x trials

addpath cd /home/zhibin/Documents/GitHub/matlab/external/arfit
pmin=1; % minimal order
pmax=10; % maximal order
[w, A, C, sbc, fpe, th]=arfit(v, pmin, pmax) % , 'sbc', 'zero')
% v(k,:)' = w' + A1*v(k-1,:)' +...+ Ap*v(k-p,:)' + noise(C)

%  For a bivariate AR(2process)
%  v(k,:)' = w' + A1*v(k-1,:)' + A2*v(k-2,:)' + eta(k,:)',
w; % least squares estimates of the intercept vector w (m x 1)
A; % coefficient matrices A=[A1 ... Ap]
C; % noise covariance matrix (m x m)
sbc; % Schwarz's Bayesian Criterion
fpe; % logarithm of Akaike's Final Prediction Error
th; % contains information needed for the computation of confidence intervals

% locate p for the minimal value
[minSBC,indSBC]=min(sbc)
[minFPE,indFPE]=min(fpe)

% try arconf
[Aerr, werr]=arconf(A, C, w, th)
Aerr; % margins of error
werr; % (A +/- Aerr) and (w +/- werr) are approximate 95% confidence intervals for the elements of the coefficient matrix A

% See if the intercept cross zero
errorbar(1:length(w),w,werr,'ko');
xlim([0 length(w)+1])

% Granger Causality
% https://www.youtube.com/watch?v=XqsSB_vpHLs
GC=

%% try out tsdata_to_var and tsdata_to_autocov
v=permute(Int12LR,[2,1,3]);
X=v; % multi-trial time series data
p=5; % model order (number of lags)
regmode='LWR' ; % regression mode: 'LWR' (default) or 'OLS'
[A,SIG,E] = tsdata_to_var(X,p,regmode);
A; % VAR coefficients matrix
SIG; % residuals covariance matrix
E; % residuals time series

cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/VAR_Granger
open mvgc_demo_autocov_tryout.m


X=v; % multi-trial time series data
q=5; % number of lags
G = tsdata_to_autocov(X,q)
G; % sample autocovariance sequence

%% Issues
% each trial has to be the same length for the input as a 3d matrix (1:100)

% arfit: Time series too short 2x100x12. 

% arfit: A=[] 100x2x12. 
% arconf: Index in position 1 exceeds array bounds.



