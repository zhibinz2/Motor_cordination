function [reg,ypred_fit,X_MC,Y_MC] = mynpls_pred(X,Y) %,Fac)
% [R2,reg,ypred,Xfactors,Yfactors,Core,B] = npls_pred(X,Y) 
% Drop 1 and cross-validate pls model
%X = EEG measures organized as subjects x channels (or channel pairs)
%Y = Dependent measures organized as subjects x measures
%obviously number of subjects is the same in both matrices and the subjects need to line up. 
%Fac=number of factors to use for prediction.  
%HINT: run it once with Fac = 1.  You will get a table that tells you how the fir went. Then rerun with more factors 
%as appropriate. Remember more factors always gives a better fit but not necessarily a better prediction.  Only add factors that 
%add a lot to the fit. 
%R2 = prediction error  -  a perfect model would return close to 100.  Unlike usual R2 can be rnegative.   
%reg - cell array that contains the regression coefficients between x and y.  ref{Fac} is the one correspondiong to the cross-validation. 
%if your data was at different channels this tells you the weight on each channel.  You can make a map of it. 
% Build PLS model

% Mean center raw data
X_MC=mean_center(X);
Y_MC=mean_center(Y);

disp('/////Calculating fitted model/////')

[Z_fit] = dosc(X_MC,Y_MC,1,1E-3);   % Direct orthogonal signal correction
[Xfactors,Yfactors,Core,B,ypred_fit,ssx,ssy,reg]=npls(Z_fit,Y_MC,10,1);% PLS regression model
% OUTPUTS
% Xfactors; Yfactors; Core; B; ypred_fit; ssx; ssy; reg
%     clear ypred_fit ssx ssy 

% 
% disp('/////Calculating cross-validation model/////')
% 
% if nargin<3
%     Fac=input('# of components to include (80% variance accounted):  ');
% end
% 
% n=size(X,1);
% ypred=zeros(n,1);
% T=zeros(n,Fac);
% 
% for i=1:n;
%     % Remove i from complete data set
%     X_set=X(setdiff(1:n,i),:);
%     Y_set=Y(setdiff(1:n,i),:);
%     
%     % Calculate mean values
%     X_mean=mean(X_set,1);
%     Y_mean=mean(Y_set);
%     
%     % Mean center data subset
%     X_set=mean_center(X_set);
%     Y_set=mean_center(Y_set);
%     
%     % Direct orthogonal signal correction
%     [Z_set,W_dosc,P_dosc] = dosc(X_set,Y_set,1,1E-3);
%     
%     % Build NPLS prediction model
%     [Xfacs,Yfacs,Cr,b] = npls(Z_set,Y_set,10,nan);
%     
%     % Remove X mean for X(i)
%     X_fix=X(i,:)-X_mean;
%     % Remove DOSC component from X(i)
%     [Znew]=dosc_pred(X_fix,W_dosc,P_dosc);
%     
%     % Prediction using NPLS model
%     [ypred_tp,T_tp]=npred(Znew,Fac,Xfacs,Yfacs,Cr,b,1);
%     
%     % Add mean center value to predicted value
%     ypred(i,1)=ypred_tp+Y_mean;
%     T(i,:)=T_tp;
%     
%     clear X_set Y_set X_mean Y_mean X_fix
%     clear Z_set W_dosc P_dosc
%     clear Xfacs Yfacs Cr b
%     clear Znew
%     clear ypred_tp T_tp   
% end
% 
% % Calculate variance accounted for in cross-validated model

% ssEr=sum(sum((Y-ypred_fit).^2))% ssEr=sum(sum((Y-ypred).^2))
% ssOr=sum(sum((Y-mean(Y)).^2))
% R2=100*(1-ssEr/ssOr)

end

