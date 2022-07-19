function [F,pval,sig] = mygranger(y1,y2);
% This function estimate and remove d in y1 and y2
% Then compute Granger Causality

% d estimate
d_min=-0.5; d_max=1;
est_d1=d_estimation(y1,d_min,d_max);
est_d2=d_estimation(y2,d_min,d_max);
% d filter
[Ytmp1,~]=remove_d(y1',est_d1); %filtered data
[Ytmp2,~]=remove_d(y2',est_d2); %filtered data

% Parameters

ntrials   = 1;     % number of trials
nobs      =  min([length(Ytmp1) length(Ytmp2)]);  % number of observations per trial

regmode   = 'LWR';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
momax     = 20;     % maximum model order for model order estimation

acmaxlags = 100;   % maximum autocovariance lags (empty for automatic calculation)

tstat     = 'F';    % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDRD'; % multiple hypothesis test correction (see routine 'significance')

fs        = 1.3;    % sample rate (Hz)
fres      = [];     % frequency resolution (empty for automatic calculation)

X=[];
X=cat(1,Ytmp1(1:nobs)',Ytmp2(1:nobs)');
X=cat(3,X); % over n trials
seed      = 0;      % random seed (0 for unseeded)

% Model order estimation (<mvgc_schema.html#3 |A2|>)
[AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X,momax,icregmode);

% % Plot information criteria.
% figure(1); clf;
% plot_tsdata([AIC BIC]',{'AIC','BIC'},1/fs);
% title('Model order estimation');
% fprintf('\nbest model order (AIC) = %d\n',moAIC);
% fprintf('best model order (BIC) = %d\n',moBIC);

% Select model order.
if     strcmpi(morder,'actual')
    morder = amo;
%     fprintf('\nusing actual model order = %d\n',morder);
elseif strcmpi(morder,'AIC')
    morder = moAIC;
%     fprintf('\nusing AIC best model order = %d\n',morder);
elseif strcmpi(morder,'BIC')
    morder = moBIC;
%     fprintf('\nusing BIC best model order = %d\n',morder);
else
%     fprintf('\nusing specified model order = %d\n',morder);
end

% VAR model estimation (<mvgc_schema.html#3 |A2|>)

% Estimate VAR model of selected order from data.
[A,SIG] = tsdata_to_var(X,morder,regmode);
assert(~isbad(A),'VAR estimation failed - bailing out');
% Report information on the estimated VAR, and check for errors.
% _IMPORTANT:_ We check the VAR model for stability and symmetric
% positive-definite residuals covariance matrix. _THIS CHECK SHOULD ALWAYS BE
% PERFORMED!_ - subsequent routines may fail if there are errors here. If there
% are problems with the data (e.g. non-stationarity, colinearity, etc.) there's
% also a good chance they'll show up at this point - and the diagnostics may
% supply useful information as to what went wrong.
info = var_info(A,SIG);
assert(~info.error,'VAR error(s) found - bailing out');

% Granger causality calculation: time domain  (<mvgc_schema.html#3 |A13|>)

% Calculate time-domain pairwise-conditional causalities from VAR model parameters
% by state-space method [4]. The VAR model is transformed into an equivalent state-
% space model for computation. Also return p-values for specified test (F-test or
% likelihood-ratio test; this is optional - if p-values are not required, then it
% is not necessary to supply time series |X|, regression mode |regmode|, or test
% specification |tstat|).
[F,pval] = var_to_pwcgc(A,SIG,X,regmode,tstat);
% Check for failed GC calculation
assert(~isbad(F,false),'GC calculation failed - bailing out');
% Significance-test p-values, correcting for multiple hypotheses.
sig = significance(pval,alpha,mhtc);
end

