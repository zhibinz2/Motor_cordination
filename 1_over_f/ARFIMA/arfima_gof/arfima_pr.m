%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     - - ARFIMA model estimation by periodogram regression. - -
%% Estimates the d fractional differencing parameter by log-periodogram 
%% regression, then calculates the ARMA part by filtering the series.
%%
%% Can be used separately, or - my advice - with the ML estimator; works fine
%% in itself, but in many cases, the Whittle approximate ML ameliorates the 
%% estimation results substantially, particularly if p and / or q large. 
%%
%%     [params_perreg,params_MLE,arma_part,wn_var] = arfima_pr(Z,esttype,arma_part)
%%
%%   [params_perreg,params_MLE,arma_part,wn_var] = arfima_pr(Z,esttype,arma_part,uroottest)
%%
%% Input     : Z - the long memory time series
%%
%%             esttype - estimation type, 'GPH' = Geweke-Porter-Hudak
%%                                                periodogram regression only
%%
%%                                        'WHI' = executes the Geweke-Porter-Hudak
%%                                                method,and then, starting
%%                                                on the estimated parameters, the
%%                                                Whittle MLE estimation
%%
%%             arma_part - [p q] where p is the AR part q is the MA 
%%                         part's length. [0 0] for the ARFIMA(0,d,0) case.
%%                         Else, [] if it is to be determined by the
%%                         user with the help of the ARMA diagnostics. 
%%
%%             uroottest - optional: 'AUGDF0'or ... 'AUGDF3'.See augdfautolag's help. 
%%                         Include only if you'd like to test the series for
%%                         unit root non-stationarity. Default is no unit root
%%                         test,estimation goes on anyway
%%                         
%% Output    : params_perreg - parameters estimated by the periodogram
%%                             regression; standard errors underneath
%%
%%             params_MLE - parameters estimated by the Whittle method;
%%                          void if method was set to 'GPH'; standard errors 
%%                          underneath
%%
%%             the estimated parameters: d phi_1 ... phi_p theta_1 ...
%%             ... theta_q 
%%
%%             Note that in comparison with MATLAB, the estimation 
%%             gives the opposite sign to the AR coefficients,namely eg.: 
%%
%%                         X(t) = phi_1 * X(t-1) + e(t)  in the AR(1) case 
%%
%%             arma_part - the chosen ARMA(p,q) model's AR and MA lags, [p q]
%%
%%             wn_var - the variance of the white noise sequence, namely:
%%                      (1-phi_1 * L - ... - phi_p * L^p)((1-L)^d) Z 
%%                                           =
%%                      (1-theta_1*L - ... - theta_q * L^q) e(t); 
%%                      
%%                            where e(t) ~ (0,wn_var)
%%
%%
%% References: Beran, Jan: Statistics for Long memory processes. 
%%                         Chapman & Hall, New York(1994)
%%
%%             Brockwell, P. - Davis, R.: Time Series Analysis.Ch.12.4.
%%                         Springer Verlag, New York(1993)
%%
%% Copyright: György Inzelt 28-01-2011
%%            inzeltgy@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[params_perreg,params_MLE,arma_part,wn_var] = arfima_pr(Z,esttype,arma_part,uroottest)
clc
switch nargin
    case 3 %if uroottest not included
        cont_ind =1;
    case 4 %if uroottest included
           switch uroottest
               case('AUGDF0')
                 [adfstat,adf_pval,critval] = augdfautolag(Z,0,10,'BIC');
                 if adf_pval < 0.05
                   disp('Series stationary')
                   cont_ind =1;
                    elseif adf_pval > 0.05
                   %if there is a unit root (an eigenvalue of the AR part is on the unit circle),
                   %you might consider integer differencing or to continue 
                     cont_ind = input('DGP contains unit root,algorithm might not work properly due to singularity. \n You should first integer differenciate the series  \n Continue anyway? [1 / 0]');
                 end
                case('AUGDF1')
                 [adfstat,adf_pval,critval] = augdfautolag(Z,1,10,'BIC');
                 if adf_pval < 0.05
                   disp('Series stationary')
                   cont_ind =1;
                    elseif adf_pval > 0.05
                   %if there is a unit root (an eigenvalue of the AR part is on the unit circle),
                   %you might consider integer differencing or to continue 
                     cont_ind = input('DGP contains unit root,algorithm might not work properly due to singularity. \n You should first integer differenciate the series  \n Continue anyway? [1 / 0]');
                 end 
                case('AUGDF2')
                 [adfstat,adf_pval,critval] = augdfautolag(Z,2,10,'BIC');
                 if adf_pval < 0.05
                   disp('Series stationary')
                   cont_ind =1;
                    elseif adf_pval > 0.05
                   %if there is a unit root (an eigenvalue of the AR part is on the unit circle),
                   %you might consider integer differencing or to continue 
                     cont_ind = input('DGP contains unit root,algorithm might not work properly due to singularity. \n You should first integer differenciate the series  \n Continue anyway? [1 / 0]');
                 end
                case('AUGDF3')
                 [adfstat,adf_pval,critval] = augdfautolag(Z,3,10,'BIC');
                 if adf_pval < 0.05
                   disp('Series stationary')
                   cont_ind =1;
                    elseif adf_pval > 0.05
                   %if there is a unit root (an eigenvalue of the AR part is on the unit circle),
                   %you might consider integer differencing or to continue 
                     cont_ind = input('DGP contains unit root,algorithm might not work properly due to singularity. \n You should first integer differenciate the series  \n Continue anyway? [1 / 0]');
                 end
               otherwise
                cont_ind =1;   
           end
end
%continue or not based on input argument numbers and / or user preference
switch cont_ind
    case(1)%continue
[Pxx,w] = periodogram(Z);
Y = log(Pxx);
w(1,1) = w(2,1)/2;
%preliminary calculations
x = exp(-sqrt(-1)*w);  
x = abs(1-x);
%zero should not be there
x(1,1) = abs(1-exp(-sqrt(-1)*w(1,1)));
x(2:length(x),1) = x(2:length(x),1).^2;
x(1,1) = log(x(1,1)^2);
x(2:length(x),1) = log(x(2:length(x),1));
xj = [ones(length(x),1)  x];

switch esttype
%only the periodogram regression - dm is the estimated fdiff
%parameter,depending on the chosen bandwidth
    case('GPH')
        disp('Method chosen: Geweke-Porter-Hudak log-periodogram regression')
%these lines are here to demonstrate the d parameter's sensitivity
dm = zeros(3,(length(Y))-20);
for m = 20:(length(Y))
    B = zeros(2,1);
    [B,BINT,R,RINT,STATS] = regress(Y(1:m,1),xj(1:m,1:2));
    vardhatr = (imag(log(-1))^2)/(6*sum((x(1:m,1) - mean(x(1:m,1))).^2));
    sdhatr = sqrt(vardhatr);
    dm(1:3,m-19) = [B ;sqrt(STATS(1,4))];
end

ind_d = ceil(length(Z)^0.8);
dhat2 = -dm(2,ind_d);
vardhat2 = (imag(log(-1))^2)/(6*sum((x(1:ind_d,1) - mean(x(1:ind_d,1))).^2));
sdhat2 = sqrt(vardhat2);

Jx = fft(Z);
w1 = (0:2*pi/length(Jx):2*pi - 2*pi/length(Jx))';
w1(1,1) = w1(2,1)/2;

%the C++ (1-L)^d filter
%Uhat = dFilter(dhat2,Jx,w1);
%Uhat = Uhat/length(Uhat);

%the MATLAB (1-L)^d filter
Uhat = d_filter(dhat2,Jx,w1);

%estimation of the ARMA part
opt0 = optimset('Display','off');
       if isempty(arma_part) ==1
         
             for p = 0:2%MODIFY HERE up to (L^p)
                 for q = 0:2%MODIFY HERE up to (L^q)
                     [arma_pq,LgL,ers,se_reg, diags, cov_rob, cov] = armaxfilter(Uhat,1,[1:p],[1:q],[],[],opt0);
                     stdP = diag(sqrt(cov_rob));
                     size_arma_pq = size(arma_pq);
                     if size_arma_pq(1) < size_arma_pq(2),arma_pq = arma_pq';end
                     
                     if p~=0 && q~=0
                        s(p+1,q+1) = struct('ARMA',[p q],'Constant',[arma_pq(1);stdP(1)],'AR',[arma_pq(2:2+p-1)';stdP(2:2+p-1)'],'MA',[arma_pq(2+p:2+p+q-1)';stdP(2+p:2+p+q-1)'],'LogL',LgL,'HQC',diags.HQC,'SBIC',diags.SBIC,'WhiteNoiseVar',var(ers));
                     
                        elseif p==0 && q~=0
                        s(p+1,q+1) = struct('ARMA',[p q],'Constant',[arma_pq(1);stdP(1)],'AR',[],'MA',[arma_pq(2+p:length(arma_pq))';stdP(2+p:length(arma_pq))'],'LogL',LgL,'HQC',diags.HQC,'SBIC',diags.SBIC,'WhiteNoiseVar',var(ers));
                         
                        elseif p~=0 && q==0
                        s(p+1,q+1) = struct('ARMA',[p q],'Constant',[arma_pq(1);stdP(1)],'AR',[arma_pq(2:2+p-1)';stdP(2:2+p-1)'],'MA',[],'LogL',LgL,'HQC',diags.HQC,'SBIC',diags.SBIC,'WhiteNoiseVar',var(ers));
                    
                        elseif p==0 && q==0
                        s(p+1,q+1) = struct('ARMA',[p q],'Constant',[arma_pq(1);stdP(1)],'AR',[],'MA',[],'LogL',LgL,'HQC',diags.HQC,'SBIC',diags.SBIC,'WhiteNoiseVar',var(ers));
                     
                     end
                     M = sprintf('%s%2d%2d\n%s%f\n%20f\n%s%f\n%s%f\n%s%f\n','ARMA order:',s(p+1,q+1).ARMA(1),s(p+1,q+1).ARMA(2),'Constant:  ',s(p+1,q+1).Constant(1),s(p+1,q+1).Constant(2),'Log-likelihood:      ',s(p+1,q+1).LogL,'Hannan-QuinnIC:       ',s(p+1,q+1).HQC,'Schwarz IC:           ',s(p+1,q+1).SBIC);
                     disp(M);
                     clear M;
                 end
             end
             Lags = input('Select model ([p q]):');
             arma_part=Lags;
             params_perreg =  [ [dhat2;sdhat2] s(Lags(1)+1,Lags(2)+1).AR s(Lags(1)+1,Lags(2)+1).MA ];
             wn_var = s(Lags(1)+1,Lags(2)+1).WhiteNoiseVar;
             
     elseif sum(arma_part)== 0 && isempty(arma_part) ~=1
             params_perreg = [dhat2; sdhat2];
             wn_var = var(Uhat);
             
     elseif sum(arma_part)~=0 && isempty(arma_part) ~=1
             p = arma_part(1);
             q = arma_part(2);
             [arma_pq,LgL,ers,se_reg, diags, cov_rob, cov] = armaxfilter(Uhat,1,[1:p],[1:q],[],[],opt0);
             stdP = diag(sqrt(cov_rob));
             s = struct('ARMA',[arma_part],'Constant',[arma_pq(1);stdP(1)],'AR',[arma_pq(2:2+p-1)';stdP(2:2+p-1)'],'MA',[arma_pq(2+p:length(arma_pq))';stdP(2+p:length(arma_pq))'],'LogLikelihood',LgL,'HQC',diags.HQC,'SBIC',diags.SBIC,'WhiteNoiseVar',var(ers));
                     
             params_perreg = [ [dhat2;sdhat2] s.AR s.MA ];
             wn_var = s.WhiteNoiseVar;
      end

params_MLE = [];
%quicker periodogram regression,then the MLE
    case('WHI')
         disp('Method chosen: Whittle approximate maximum likelihood estimation')
        
         ind_d = ceil(length(Z)^0.8);
         B = zeros(2,1);
        [B,BINT,R,RINT,STATS] = regress(Y(1:ind_d,1),xj(1:ind_d,1:2));
         dhat2 = -B(2);
         vardhat2 = (imag(log(-1))^2)/(6*sum((x(1:ind_d,1) - mean(x(1:ind_d,1))).^2));
         sdhat2 = sqrt(vardhat2);

         Jx = fft(Z);
         w1 = (0:2*pi/length(Jx):2*pi - 2*pi/length(Jx))';
         w1(1,1) = w1(2,1)/2;

         %the C++ (1-L)^d filter
         %Uhat = dFilter(dhat2,Jx,w1);
         %Uhat = Uhat/length(Uhat);

         %the MATLAB (1-L)^d filter
         Uhat = d_filter(dhat2,Jx,w1);

         %estimation of the ARMA part
         opt0 = optimset('Display','off');
         if isempty(arma_part) ==1
              for p = 0:2%MODIFY HERE up to custom AR lags
                 for q = 0:2%MODIFY HERE up to custom MA lags
                     [arma_pq,LgL,ers,se_reg, diags, cov_rob, cov] = armaxfilter(Uhat,1,[1:p],[1:q],[],[],opt0);                 
                     stdP = diag(sqrt(cov_rob));
                     size_arma_pq = size(arma_pq);
                     if size_arma_pq(1) < size_arma_pq(2),arma_pq = arma_pq';end
                     
                     if p~=0 && q~=0
                        s(p+1,q+1) = struct('ARMA',[p q],'Constant',[arma_pq(1);stdP(1)],'AR',[arma_pq(2:2+p-1)';stdP(2:2+p-1)'],'MA',[arma_pq(2+p:2+p+q-1)';stdP(2+p:2+p+q-1)'],'LogL',LgL,'HQC',diags.HQC,'SBIC',diags.SBIC,'WhiteNoiseVar',var(ers));
                     
                        elseif p==0 && q~=0
                        s(p+1,q+1) = struct('ARMA',[p q],'Constant',[arma_pq(1);stdP(1)],'AR',[],'MA',[arma_pq(2+p:length(arma_pq))';stdP(2+p:length(arma_pq))'],'LogL',LgL,'HQC',diags.HQC,'SBIC',diags.SBIC,'WhiteNoiseVar',var(ers));
                         
                        elseif p~=0 && q==0
                        s(p+1,q+1) = struct('ARMA',[p q],'Constant',[arma_pq(1);stdP(1)],'AR',[arma_pq(2:2+p-1)';stdP(2:2+p-1)'],'MA',[],'LogL',LgL,'HQC',diags.HQC,'SBIC',diags.SBIC,'WhiteNoiseVar',var(ers));
                    
                        elseif p==0 && q==0
                        s(p+1,q+1) = struct('ARMA',[p q],'Constant',[arma_pq(1);stdP(1)],'AR',[],'MA',[],'LogL',LgL,'HQC',diags.HQC,'SBIC',diags.SBIC,'WhiteNoiseVar',var(ers));
                     
                     end
                     M = sprintf('%s%2d%2d\n%s%f\n%20f\n%s%f\n%s%f\n%s%f\n','ARMA order:',s(p+1,q+1).ARMA(1),s(p+1,q+1).ARMA(2),'Constant:  ',s(p+1,q+1).Constant(1),s(p+1,q+1).Constant(2),'Log-likelihood:      ',s(p+1,q+1).LogL,'Hannan-QuinnIC:       ',s(p+1,q+1).HQC,'Schwarz IC:           ',s(p+1,q+1).SBIC);
                     disp(M);
                     clear M;
                 end
             end
             Lags = input('Select ARMA model ([p q]):');
             arma_part = Lags;
             
             if sum(Lags)~=0 && Lags(1)~=0 && Lags(2) ~=0
             params_perreg =  [ [dhat2;sdhat2] s(Lags(1)+1,Lags(2)+1).AR s(Lags(1)+1,Lags(2)+1).MA ];
             startingval = [ s(Lags(1)+1,Lags(2)+1).AR(1,:) s(Lags(1)+1,Lags(2)+1).MA(1,:) ];
             arma_part = Lags;
             elseif sum(Lags)~=0 && Lags(1)~=0 && Lags(2) ==0
             params_perreg =  [ [dhat2;sdhat2]  s(Lags(1)+1,Lags(2)+1).AR ];
             startingval =  s(Lags(1)+1,Lags(2)+1).AR(1,:) ;
             arma_part = Lags;
             elseif sum(Lags)~=0 && Lags(1)==0 && Lags(2) ~=0
             params_perreg =  [ [dhat2;sdhat2]  s(Lags(1)+1,Lags(2)+1).MA ];
             startingval =  s(Lags(1)+1,Lags(2)+1).MA(1,:) ;
             arma_part = Lags;
             elseif sum(Lags)==0
             params_perreg = [dhat2 ; sdhat2 ];
             arma_part = [0 0];
             startingval = [];   
             end
           
         elseif sum(arma_part)== 0 && isempty(arma_part) ~=1
             params_perreg = [dhat2 ; sdhat2 ];
             arma_part = [0 0];
             startingval = [];
             
         elseif sum(arma_part)~=0 && isempty(arma_part) ~=1
             p = arma_part(1);
             q = arma_part(2);
             [arma_pq,LgL,ers,se_reg, diags, cov_rob, cov] = armaxfilter(Uhat,1,[1:p],[1:q],[],[],opt0);
              stdP = diag(sqrt(cov_rob));
              s = struct('ARMA',[arma_part],'Constant',[arma_pq(1);stdP(1)],'AR',[arma_pq(2:2+p-1)';stdP(2:2+p-1)'],'MA',[arma_pq(2+p:length(arma_pq))';stdP(2+p:length(arma_pq))'],'LogLikelihood',LgL,'HQC',diags.HQC,'SBIC',diags.SBIC);
                     
             params_perreg = [ [dhat2;sdhat2] s.AR s.MA ];
             startingval = [ s.AR(1,:) s.MA(1,:) ];
        end
 
lpr = size(params_perreg); 
%changed Lb&Ub to include stationary but unusual AR
%models,eg.(1-1.3L+0.35L^2) X(t) = e(t).
Lb = -1.99*ones(1,lpr(2));
Ub =  1.99*ones(1,lpr(2));%works even with d > 0.5 (the non-stationary case)
opt1 = optimset('Algorithm','active-set');

[params_MLE,LL,EFLAG,OUTPUT,LAMBDA,GRAD,D2LL] = fmincon(@(params) arfima_whittle(params,Z,arma_part),[0 startingval],[],[],[],[],Lb,Ub,[],opt1); 

%if the algorithm got stuck at a suboptimal value
if sum(LAMBDA.upper)~=0 || sum(LAMBDA.lower)~=0    
[params_MLE,LL,EFLAG,OUTPUT,LAMBDA,GRAD,D2LL] = fmincon(@(params) arfima_whittle(params,Z,arma_part),[dhat2 startingval],[],[],[],[],Lb,Ub,[],opt1); 
end
%calculating the white noise variance and the MLE parameter estimates' std errors
params_MLE_stds =  diag(sqrt((1/length(Z))*inv(D2LL)));
wn_var = LL;
params_MLE = [params_MLE; params_MLE_stds'];

    otherwise %if not 'GPH' or 'WHI'
        disp('Unknown method') 
        params_perreg = [];
        params_MLE = [];
        arma_part =[];
        wn_var=[];
end
    otherwise %if terminated due to a unit root
        disp('Estimation terminated due to detected unit root in the data generating process.')
        params_perreg = [];
        params_MLE = [];
        arma_part = [];
        wn_var=[];
end
end