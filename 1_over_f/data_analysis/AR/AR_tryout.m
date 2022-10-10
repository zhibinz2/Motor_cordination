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
runid=num2str(seeds(s,:));
load(['/ssd/zhibin/1overf/' runid '_2P/Cleaned_data/clean_' runid '.mat' ]) 
intervals_H_removed

%% try out arfit and arconf
% one tapping int series from L
b=1;
Int1L=[intervals_H_removed{b}(:,1)];

% 1 int series from L & R
Int1LR=intervals_H_removed{b}(:,1:2);

% add 3rd dimension as trials
% replicate 12 trials
% Int12LR=repmat(Int1LR,1,1,12);
Int12LR=nan(100,2,12);
for b=1:12
    Int12LR(:,:,b)=[intervals_H_removed{b}(1:100,1:2)];
end

% try arfit
% select v to try out
v=Int1L; 
v=Int1LR;
v=Int12LR;
v; % time series data x observations x trials

pmin=0; % minimal order
pmax=4; % maximal order
[w, A, C, sbc, fpe, th]=arfit(v, pmin, pmax, 'sbc', 'zero')
% v(k,:)' = w' + A1*v(k-1,:)' +...+ Ap*v(k-p,:)' + noise(C)
w; % least squares estimates of the intercept vector w (m x 1)
A; % coefficient matrices A=[A1 ... Ap]
C; % noise covariance matrix (m x m)
sbc; % Schwarz's Bayesian Criterion
fpe; % logarithm of Akaike's Final Prediction Error
th; % contains information needed for the computation of confidence intervals

% try arconf
[Aerr, werr]=arconf(A, C, w, th)
Aerr; % margins of error
werr; % (A +/- Aerr) and (w +/- werr) are approximate 95% confidence intervals for the elements of the coefficient matrix A

%% try out tsdata_to_var and tsdata_to_autocov

X=v; % multi-trial time series data
p=5; % model order (number of lags)
regmode='LWR' ; % regression mode: 'LWR' (default) or 'OLS'
[A,SIG,E] = tsdata_to_var(X,p,regmode)
A; % VAR coefficients matrix
SIG; % residuals covariance matrix
E; % residuals time series

cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/VAR_Granger
open mvgc_demo_autocov_tryout.m



X=v; % multi-trial time series data
q=5; % number of lags
G = tsdata_to_autocov(X,q)
G; % sample autocovariance sequence



