cd /ssd/zhibin/1overf/20220610_2P/Segmented_data/1_50Hz_ICAautomized
clear
load('BP20220610.mat')
load('FB20220610.mat')

condi=1; BP_L=BPCondi1L; FB_L=FBCondi1L; BP_R=BPCondi1R;FB_R=FBCondi1R;
condi=2; BP_L=BPCondi2L; FB_L=FBCondi2L; BP_R=BPCondi2R;FB_R=FBCondi2R;
condi=3; BP_L=BPCondi3L; FB_L=FBCondi3L; BP_R=BPCondi3R;FB_R=FBCondi3R;
condi=4; BP_L=BPCondi4L; FB_L=FBCondi4L; BP_R=BPCondi4R;FB_R=FBCondi4R;

% syncopation expt
conditionNames={'uncoupled' 'L-lead' 'R-lead' 'mutual-1.3Hz'};
% cut off from the last bottom press
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing
BP_L=cutofflast(BP_L);
BP_R=cutofflast(BP_R);

y1=Calinterval(BP_L')./2;
y2=Calinterval(BP_R')./2;

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
open mvgc_demo_statespace.m

%% Parameters

ntrials   = 1;     % number of trials
nobs      = 575;   % number of observations per trial

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
X=cat(1,Ytmp1',Ytmp2(1:575)');
X=cat(3,X); % over n trials
seed      = 0;      % random seed (0 for unseeded)

%% Model order estimation (<mvgc_schema.html#3 |A2|>)

% Calculate information criteria up to specified maximum model order.

ptic('\n*** tsdata_to_infocrit\n');
[AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X,momax,icregmode);
ptoc('*** tsdata_to_infocrit took ');

% Plot information criteria.

figure(1); clf;
plot_tsdata([AIC BIC]',{'AIC','BIC'},1/fs);
title('Model order estimation');

fprintf('\nbest model order (AIC) = %d\n',moAIC);
fprintf('best model order (BIC) = %d\n',moBIC);

% We don't know the acutal model order

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
assert(~isbad(A),'VAR estimation failed - bailing out');
ptoc;

% Report information on the estimated VAR, and check for errors.
%
% _IMPORTANT:_ We check the VAR model for stability and symmetric
% positive-definite residuals covariance matrix. _THIS CHECK SHOULD ALWAYS BE
% PERFORMED!_ - subsequent routines may fail if there are errors here. If there
% are problems with the data (e.g. non-stationarity, colinearity, etc.) there's
% also a good chance they'll show up at this point - and the diagnostics may
% supply useful information as to what went wrong.

info = var_info(A,SIG);
assert(~info.error,'VAR error(s) found - bailing out');


%% Granger causality calculation: time domain  (<mvgc_schema.html#3 |A13|>)

% Calculate time-domain pairwise-conditional causalities from VAR model parameters
% by state-space method [4]. The VAR model is transformed into an equivalent state-
% space model for computation. Also return p-values for specified test (F-test or
% likelihood-ratio test; this is optional - if p-values are not required, then it
% is not necessary to supply time series |X|, regression mode |regmode|, or test
% specification |tstat|).

ptic('*** var_to_pwcgc... ');
[F,pval] = var_to_pwcgc(A,SIG,X,regmode,tstat);
ptoc;

% Check for failed GC calculation

assert(~isbad(F,false),'GC calculation failed - bailing out');

% Significance-test p-values, correcting for multiple hypotheses.

sig = significance(pval,alpha,mhtc);

% Plot time-domain causal graph, p-values and significance.

figure(2); clf;
sgtitlex('Pairwise-conditional Granger causality - time domain');
subplot(1,3,1);
plot_pw(F);
title('Pairwise-conditional GC');
subplot(1,3,2);
plot_pw(pval);
title(['p-values (' tstat '-test)']);
subplot(1,3,3);
plot_pw(sig);
title(['Significant at \alpha = ' num2str(alpha)]);

%% Granger causality calculation: frequency domain  (<mvgc_schema.html#3 |A14|>)

% If not specified, we set the frequency resolution to something sensible. Warn if
% resolution is very large, as this may lead to excessively long computation times,
% and/or out-of-memory issues.

if isempty(fres)
    fres = 2^nextpow2(info.acdec); % based on autocorrelation decay; alternatively, you could try fres = 2^nextpow2(nobs);
	fprintf('\nfrequency resolution auto-calculated as %d (increments ~ %.2gHz)\n',fres,fs/2/fres);
end
if fres > 20000 % adjust to taste
	fprintf(2,'\nWARNING: large frequency resolution = %d - may cause computation time/memory usage problems\nAre you sure you wish to continue [y/n]? ',fres);
	istr = input(' ','s'); if isempty(istr) || ~strcmpi(istr,'y'); fprintf(2,'Aborting...\n'); return; end
end

% Calculate spectral pairwise-conditional causalities at given frequency
% resolution by state-space method.

ptic('\n*** var_to_spwcgc... ');
f = var_to_spwcgc(A,SIG,fres);
assert(~isbad(f,false),'spectral GC calculation failed - bailing out');
ptoc;

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
