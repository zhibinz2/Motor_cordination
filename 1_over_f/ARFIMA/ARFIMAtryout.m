%% simulation function 1 
% https://www.mathworks.com/matlabcentral/fileexchange/25611-arfima-simulations
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/ARFIMA
cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/ARFIMA

% [Z] = ARFIMA_SIM(N,F,O,d,stdx,er)
%%% INPUTS
%%%->N =  # % Length of the time series we would like to generate  
N=240;
%%%->F = [ F1 F2 F3 .... ] % Parameters of the AR model, length(F) is the order p. Default p = 0
%%%->O = [ O1 O2 O3 .... ] % Parameters of the MA model, length(O) is the order q. Default q = 0   
%%%->d = # ; % Fractionally differencing parameter, default d = 0 
%%%->stdx = % Optional input: parameter to force the standard deviation of the
%output time series. Impose std(Z)==stdx    
%%%-->er = % Optional input: predefined time seres of white noise  

%%%%%%%%%% EXAMPLES 
%%% White noise 
[Z] = ARFIMA_SIM(N); 
plot(Z,'r.');

%%% AR(1) model 
F1=[3]; % p
[Z] = ARFIMA_SIM(N,[F1]); 
plot(Z,'r.');

%%% MA(1) model 
O1 = [1];% q
[Z] = ARFIMA_SIM(N,[],[O1]);
plot(Z,'r.');

%%% ARMA(2,2) model 
F1=0;F2=1;O1=0;O2=2;
[Z] = ARFIMA_SIM(N,[F1,F2],[O1,O2]);
plot(Z,'r.');

%%% ARFIMA(0,d,0)
d=0.2;
[Z] = ARFIMA_SIM(N,[],[],d);
plot(Z,'r.');

%%% ARFIMA(1,d,1)
F1=1;O1=1;d=0.2;
[Z] = ARFIMA_SIM(N,[F1],[O1],d);
plot(Z,'r.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% simulation function 2 
% https://www.mathworks.com/matlabcentral/fileexchange/53301-arfima-p-d-q
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/ARFIMA/dgp_arfima
cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/ARFIMA/dgp_arfima
open dgp_arfima.m

% sim=dgp_arfima(c,ARo,MAo,T,SD,d,F);
%--------------------------------------------------------------------------
% Input:
% c : Constant term
% AR: Parameters of the AR process (column vector: [AR(1),...AR(p)]')
% MA: Parameters of the MA process (column vector: [MA(1),...MA(q)]')
% T : Simulated sample size
% SD: Standard Deviation of the white noise
% d : Fractional memory parameter. d in [0,inf)
% F : Figures; plot of:
%               1.- the series 
%               2.- the sample FAC
%     [possible values: 0-1]
%               1: Figure
%--------------------------------------------------------------------------
% Output:
% sim           : Simulation of the arfima process 
%--------------------------------------------------------------------------
% EXAMPLES : 
% ARFIMA(p,d,q)
% ARFIMA(2,0.5,1) >>>>  Results=dgp_arfima(1,[0.5 0.1],[0.2],100,1,0.5);
Results=dgp_arfima(1,[0.5 0.1],[0.2],100,1,0.5,1);

% ARFIMA(0,0.8,0) >>>>  Results=dgp_arfima(0,[],[],100,1,0.8);
Results=dgp_arfima(0,[],[],100,1,0.8,1);

%% model estimation (not work)
% https://www.mathworks.com/matlabcentral/fileexchange/30238-arfima-p-d-q-estimator?s_tid=srchtitle_arfima_3
cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/ARFIMA/arfima_est_v2
% [model] = arfima_estimate(Z,modeltype,arma_part)
% Input arguments:
%   Z          : the time series input for the estimation
%   modeltype  : 'FWHI' or 'FML' for arfima(p,d,q) model estimation
%                 NOTE: it is currently implemented assuming that
%                       the noise sequence is gaussian 
%   arma_part  : the AR and MA lag respectively, eg.[1 2]
%                for ARFIMA(1,d,2) 
%
% Output arguments: 
%   model      : parameters and descriptive statistics returned
%                in a data structure 
%                mean : sample mean 
%                   d : the fractional differencing parameter
%                  AR : AR(1 ... p) parameters
%                  MA : MA(1 ... q) parameters
%              sigma2 : the gaussian noise seq.'s variance
%        loglikelihood: the log likelihood function value @params
%              akaike : the Akaike information criterion
%              errors : the estimated noise sequence (for further
%                       investigation)
%            fit_stas : the goodness of fit statistic - see: 
%                       help arfima_gof 

[model] = arfima_estimate(Z,'FML',[1 2])
% Dependent on armaxfilter.m from mfe-toolbox (https://github.com/zhibinz2/mfe-toolbox) 

%% model selection / goodness of fit (not work)
% https://www.mathworks.com/matlabcentral/fileexchange/30355-arfima-p-d-q-goodness-of-fit-test?s_tid=srchtitle_arfima_4
cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/ARFIMA/arfima_gof
open arfima_pr.m
% [params_perreg,params_MLE,arma_part,wn_var] = arfima_pr(Z,esttype,arma_part,uroottest)
% [params_perreg,params_MLE,arma_part,wn_var] = arfima_pr(Z,esttype,arma_part)

% Input     : Z - the long memory time series
%             esttype - estimation type, 'GPH' = Geweke-Porter-Hudak
%                                                periodogram regression only
%                                        'WHI' = executes the Geweke-Porter-Hudak
%                                                method,and then, starting
%                                                on the estimated parameters, the
%                                                Whittle MLE estimation
%             arma_part - [p q] where p is the AR part q is the MA 
%                         part's length. [0 0] for the ARFIMA(0,d,0) case.
%                         Else, [] if it is to be determined by the
%                         user with the help of the ARMA diagnostics. 
%             uroottest - optional: 'AUGDF0'or ... 'AUGDF3'.See augdfautolag's help. 
%                         Include only if you'd like to test the series for
%                         unit root non-stationarity. Default is no unit root
%                         test,estimation goes on anyway                         
% Output    : params_perreg - parameters estimated by the periodogram
%                             regression; standard errors underneath
%             params_MLE - parameters estimated by the Whittle method;
%                          void if method was set to 'GPH'; standard errors 
%                          underneath
%             the estimated parameters: d phi_1 ... phi_p theta_1 ...
%             ... theta_q 
%             Note that in comparison with MATLAB, the estimation 
%             gives the opposite sign to the AR coefficients,namely eg.: 
%                         X(t) = phi_1 * X(t-1) + e(t)  in the AR(1) case 
%             arma_part - the chosen ARMA(p,q) model's AR and MA lags, [p q]
%             wn_var - the variance of the white noise sequence, namely:
%                      (1-phi_1 * L - ... - phi_p * L^p)((1-L)^d) Z 
%                                           =
%                      (1-theta_1*L - ... - theta_q * L^q) e(t);                       
%                            where e(t) ~ (0,wn_var)

[params_perreg,params_MLE,arma_part,wn_var] = arfima_pr(Z,'GPH',[],'AUGDF0')


% **************************************************
open arfima_gof.m
%  [pval,Tn] = arfima_gof(Z,params,arma_part,wn_var)
%  [pval,Tn] = arfima_gof(Z,params,arma_part,wn_var,display)
% Input: Z - the time series
%        params - a size (1,1+p+q) parameter vector with 
%                 d phi_1 ... phi_1 ... phi_p theta_1 ...
%                 ... theta_q ; d =0 for ARMA(p,q) models
%        arma_part - [p q],where p is the AR and q is the MA lag length  
%        wn_var - the estimated white noise sequence's variance 
%        display - optional: whether to display the estimated
%                  autocovariance and autocorrelation of the 
%                  (hopefully) white noise process. 'DISP' if
%                  you'd like to display. Default is: no plot.
% Output: pval - the p-value of the test statistics
%         Tn - the test statistics evaluated @params 

[pval,Tn] = arfima_gof(Z,params,arma_part,wn_var,'DISP')

%%
