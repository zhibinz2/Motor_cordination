function [Ytmp1,Ytmp2] = rm_est_d(y1,y2)
% estimate d and remove d from 2 time series
% y1 y2 must be a row vector
% if not row vector, we now convert it
if iscolumn(y1);    
    y1=y1';
end;
if iscolumn(y2);    
    y2=y2';
end;
% d estimate
d_min=-0.5; d_max=1;
est_d1=d_estimation(y1,d_min,d_max);
est_d2=d_estimation(y2,d_min,d_max);
%d filter
[Ytmp1,~]=remove_d(y1',est_d1); %filtered data
[Ytmp2,~]=remove_d(y2',est_d2); %filtered data
end
