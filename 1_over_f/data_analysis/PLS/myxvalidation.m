function [ssEr,R2] = myxvalidation(X,Y,Fac,nSet);
% Do cross-validation 
% nSet: every nSet elements from the same subject were removed in each loop

disp('/////Calculating cross-validation model/////')

if nargin<3
    Fac=input('# of components to include (80% variance accounted):  ');
end

n=size(X,1);
% ypred=zeros(8,1);
% Y_means=zeros(8,1);
T=zeros(n,Fac);

% % Remove 3 from complete data set (every 3 were from the same subject)
% nSet=6;

for i=1:(n/nSet); %  #(n/nSet) fold cross validation
    X_set=X(setdiff(1:n,((i-1)*nSet+1):((i-1)*nSet+nSet)),:);
    Y_set=Y(setdiff(1:n,((i-1)*nSet+1):((i-1)*nSet+nSet)),:);

    % Calculate mean values
    X_mean=mean(X_set,1);
    Y_mean=mean(Y_set);
    Y_means(i)=Y_mean;

    % Mean center data subset
    X_set=mean_center(X_set);
    Y_set=mean_center(Y_set);
    
    % Direct orthogonal signal correction
    [Z_set,W_dosc,P_dosc] = dosc(X_set,Y_set,1,1E-3);
    
    % Build NPLS prediction model
    [Xfacs,Yfacs,Cr,b] = npls(Z_set,Y_set,10,nan);
    
    % Remove X mean for X(i)
    X_fix=X((((i-1)*nSet+1):((i-1)*nSet+nSet)),:)-X_mean;   
    % Remove DOSC component from X(i)
    [Znew]=dosc_pred(X_fix,W_dosc,P_dosc);
    
    % Prediction using NPLS model
    [ypred_tp,T_tp]=npred(Znew,Fac,Xfacs,Yfacs,Cr,b,1);
    
    % Add mean center value to predicted value
    ypred((((i-1)*nSet+1):((i-1)*nSet+nSet)),1)=ypred_tp+Y_mean;  
    T(((i-1)*nSet+1):((i-1)*nSet+nSet),:)=T_tp;
    
    clear X_set Y_set X_mean Y_mean X_fix
    clear Z_set W_dosc P_dosc
    clear Xfacs Yfacs Cr b
    clear Znew
    clear ypred_tp T_tp   
end

% Calculate variance accounted for in cross-validated model
ssEr=sum(sum((Y-ypred).^2))
ssOr=sum(sum((Y-mean(Y)).^2))
R2=100*(1-ssEr/ssOr)

end

