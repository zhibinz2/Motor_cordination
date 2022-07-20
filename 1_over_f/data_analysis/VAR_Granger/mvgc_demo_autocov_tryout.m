%% load some data to test
cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/VAR_Granger
clear;
cd /ssd/zhibin/1overf/20220713_2P/Segmented_data/

open /home/zhibin/Documents/GitHub/MVGC1/demo/mvgc_demo_autocov.m
load('data_variables20220713.mat')

y1=Calinterval(BP(1).BP{10}')./sr;
y2=Calinterval(BP(2).BP{10}')./sr;

%% estimate and remove d
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/MSE-VARFI
d_min=-0.5; d_max=1;
est_d1=d_estimation(y1,d_min,d_max);
est_d2=d_estimation(y2,d_min,d_max);
%d filter
[Ytmp1,~]=remove_d(y1',est_d1); %filtered data
[Ytmp2,~]=remove_d(y2',est_d2); %filtered data
%d estimation again (Whittle semiparametric estimator)
est_d12=d_estimation(Ytmp1,d_min,d_max);
est_d22=d_estimation(Ytmp2,d_min,d_max);

figure;
subplot(2,2,1);plot(y1); title(['before d removal: est-d = ' num2str(est_d1)]);
subplot(2,2,3);plot(Ytmp1);title(['removed est-d: Now d= ' num2str(est_d12)]);
subplot(2,2,2);plot(y2); title(['before d removal: est-d = ' num2str(est_d2)]);
subplot(2,2,4);plot(Ytmp2);title(['removed est-d: Now d= ' num2str(est_d22)]);

%% try out granger causality on Ytemp1 and Ytemp2
addpath(genpath('/home/zhibin/Documents/GitHub/MVGC1'));
cd /home/zhibin/Documents/GitHub/MVGC1
run startup
cd /home/zhibin/Documents/GitHub/MVGC1/demo
open mvgc_demo_autocov.m

%% Parameters

ntrials   = 1;     % number of trials
nobs      = min([length(Ytmp1) length(Ytmp2)]);  % number of observations per trial

regmode   = 'LWR';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
momax     = 20;     % maximum model order for model order estimation

% recomended acmaxlags minimum = 105
acmaxlags = 105;   % maximum autocovariance lags (empty for automatic calculation)

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

% % for VAR test data generation only
% seed      = 0;      % random seed (0 for unseeded)

%% Model order estimation (<mvgc_schema.html#3 |A2|>)

% Calculate information criteria up to specified maximum model order.

ptic('\n*** tsdata_to_infocrit\n');
[AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X,momax,icregmode);
ptoc('*** tsdata_to_infocrit took ');

% Plot information criteria.

figure(1); clf;
plot_tsdata([AIC BIC]',{'AIC','BIC'},1/fs);
title('Model order estimation');

% We don't know the acutal model order
% amo = size(AT,3); % actual model order

fprintf('\nbest model order (AIC) = %d\n',moAIC);
fprintf('best model order (BIC) = %d\n',moBIC);
% fprintf('actual model order     = %d\n',amo);

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

%% VAR model estimation (<mvgc_schema.html#3 |A2|>)

% Estimate VAR model of selected order from data.

ptic('\n*** tsdata_to_var... ');
[A,SIG] = tsdata_to_var(X,morder,regmode);
ptoc;

% Check for failed regression

assert(~isbad(A),'VAR estimation failed');

% NOTE: at this point we have a model and are finished with the data! - all
% subsequent calculations work from the estimated VAR parameters A and SIG.

%% Autocovariance calculation (<mvgc_schema.html#3 |A5|>)

% The autocovariance sequence drives many Granger causality calculations (see
% next section). Now we calculate the autocovariance sequence G according to the
% VAR model, to as many lags as it takes to decay to below the numerical
% tolerance level, or to acmaxlags lags if specified (i.e. non-empty).

ptic('*** var_to_autocov... ');
[G,info] = var_to_autocov(A,SIG,acmaxlags);
ptoc;

% The above routine does a LOT of error checking and issues useful diagnostics.
% If there are problems with your data (e.g. non-stationarity, colinearity,
% etc.) there's a good chance it'll show up at this point - and the diagnostics
% may supply useful information as to what went wrong. It is thus essential to
% report and check for errors here.

var_acinfo(info,true); % report results (and bail out on error)

%% Granger causality calculation: time domain  (<mvgc_schema.html#3 |A13|>)

% Calculate time-domain pairwise-conditional causalities - this just requires
% the autocovariance sequence.

ptic('*** autocov_to_pwcgc... ');
F = autocov_to_pwcgc(G);
ptoc;

% Check for failed GC calculation

assert(~isbad(F,false),'GC calculation failed');

% Significance test using theoretical null distribution, adjusting for multiple
% hypotheses.

pval = mvgc_pval(F,morder,nobs,ntrials,1,1,nvars-2,tstat); % take careful note of arguments!
sig  = significance(pval,alpha,mhtc);

% Plot time-domain causal graph, p-values and significance.

figure(2); clf;
sgtitlex('Pairwise-conditional Granger causality - time domain');
subplot(1,3,1);
plot_pw(F);
title('Pairwise-conditional GC');
subplot(1,3,2);
plot_pw(pval);
title('p-values');
subplot(1,3,3);
plot_pw(sig);
title(['Significant at p = ' num2str(alpha)])

% For good measure we calculate Seth's causal density (cd) measure - the mean
% pairwise-conditional causality. We don't have a theoretical sampling
% distribution for this.

cd = mean(F(~isnan(F)));

fprintf('\ncausal density = %f\n',cd);

%% Granger causality calculation: frequency domain  (<mvgc_schema.html#3 |A14|>)

% Calculate spectral pairwise-conditional causalities at given frequency
% resolution - again, this only requires the autocovariance sequence.

ptic('\n*** autocov_to_spwcgc... ');
f = autocov_to_spwcgc(G,fres);
ptoc;

% Check for failed spectral GC calculation

assert(~isbad(f,false),'spectral GC calculation failed');

% Plot spectral causal graph.

figure(3); clf;
sgtitlex('Pairwise-conditional Granger causality - frequency domain');
plot_spw(f,fs);

%% Granger causality calculation: frequency domain -> time-domain  (<mvgc_schema.html#3 |A15|>)

% Check that spectral causalities average (integrate) to time-domain
% causalities, as they should according to theory.

fprintf('\nfrequency-domain GC integration check... ');
Fint = smvgc_to_mvgc(f); % integrate spectral MVGCs
amax = maxabs(F+Fint)/2;
if amax < 1e-5; amax = 1; end % in case all GCs very small
mre = maxabs(F-Fint)/amax;
if mre < 1e-5
    fprintf('OK (maximum relative error ~ %.0e)\n',mre);
else
    fprintf(2,'WARNING: high maximum relative error ~ %.0e\n',mre);
end

%%
% <mvgc_demo.html back to top>
