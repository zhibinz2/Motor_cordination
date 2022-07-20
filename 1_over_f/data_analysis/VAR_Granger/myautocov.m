function [F,pval,sig] = myautocov(y1,y2);
% This function estimate and remove d in y1 and y2
% Then compute Granger Causality
% MVGC demo: autocov method.

%% d estimate
d_min=-0.5; d_max=1;
est_d1=d_estimation(y1,d_min,d_max);
est_d2=d_estimation(y2,d_min,d_max);
%d filter
[Ytmp1,~]=remove_d(y1',est_d1); %filtered data
[Ytmp2,~]=remove_d(y2',est_d2); %filtered data

%% Parameters
ntrials   = 1;     % number of trials
nobs      = min([length(Ytmp1) length(Ytmp2)]);  % number of observations per trial

regmode   = 'LWR';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
momax     = 20;     % maximum model order for model order estimation

% recomended acmaxlags minimum = 122 for some trials
acmaxlags = 30;   % maximum autocovariance lags (empty for automatic calculation)

tstat     = '';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')

% for plotting only
fs        = 1.3;    % sample rate (Hz)
fres      = [];     % frequency resolution (empty for automatic calculation)

% combine the time series
X=[];
X=cat(1,Ytmp1(1:nobs)',Ytmp2(1:nobs)');
X=cat(3,X);         % over n trials

nvars = size(X,1); % number of variables
%% Model order estimation (<mvgc_schema.html#3 |A2|>)
[AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X,momax,icregmode);

% % Plot information criteria.
% figure(1); clf;
% plot_tsdata([AIC BIC]',{'AIC','BIC'},1/fs);
% title('Model order estimation');

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

%% VAR model estimation (<mvgc_schema.html#3 |A2|>)
[A,SIG] = tsdata_to_var(X,morder,regmode);
% Check for failed regression
assert(~isbad(A),'VAR estimation failed');
% NOTE: at this point we have a model and are finished with the data! - all
% subsequent calculations work from the estimated VAR parameters A and SIG.

% Autocovariance calculation (<mvgc_schema.html#3 |A5|>)
% The autocovariance sequence drives many Granger causality calculations (see
% next section). Now we calculate the autocovariance sequence G according to the
% VAR model, to as many lags as it takes to decay to below the numerical
% tolerance level, or to acmaxlags lags if specified (i.e. non-empty).
[G,info] = var_to_autocov(A,SIG,acmaxlags);
% The above routine does a LOT of error checking and issues useful diagnostics.
% If there are problems with your data (e.g. non-stationarity, colinearity,
% etc.) there's a good chance it'll show up at this point - and the diagnostics
% may supply useful information as to what went wrong. It is thus essential to
% report and check for errors here.
var_acinfo(info,true); % report results (and bail out on error)

%% Granger causality calculation: time domain  (<mvgc_schema.html#3 |A13|>)
% Calculate time-domain pairwise-conditional causalities - this just requires
% the autocovariance sequence.
F = autocov_to_pwcgc(G);
% Check for failed GC calculation
assert(~isbad(F,false),'GC calculation failed');
% Significance test using theoretical null distribution, adjusting for multiple
% hypotheses.
pval = mvgc_pval(F,morder,nobs,ntrials,1,1,nvars-2,tstat); % take careful note of arguments!
sig  = significance(pval,alpha,mhtc);

end