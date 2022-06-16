function [IndHigh,IndMedium,IndLow,corrBP] = TriSeg(BP_L,BP_R,sr,win);
% x1 and x2 are two time series (x1 and x2 can be boolean of BP)
% cut x1 and x2 into 10 s per trials
SamplesLength=sr*win; % chunck of 10 seconds
nTrials=min([floor(length(BP_L)/SamplesLength) floor(length(BP_R)/SamplesLength)]);
corrBP=[];
for i=1:nTrials
    xBP_L1=[];xBP_R1=[];
    xBP_L1=BP_L((SamplesLength*(i-1)+1):(SamplesLength*i),:);
    xBP_R1=BP_R((SamplesLength*(i-1)+1):(SamplesLength*i),:);
    xBP_L2=[Calinterval(xBP_L1')'];
    xBP_R2=[Calinterval(xBP_R1')'];
    % Cutting longer one short method
%     minTap=min([length(xBP_L2) length(xBP_R2)]);
%     z = corr(xBP_L2(1:minTap), xBP_R2(1:minTap));
    % Resample method
    Maxlength=max([length(xBP_L2) length(xBP_R2)]);
    xBP_L3=resample(xBP_L2,Maxlength,length(xBP_L2));
    xBP_R3=resample(xBP_R2,Maxlength,length(xBP_R2));
    z = corr(xBP_L3, xBP_R3);
    % Save the corrcoef
    corrBP(i)=z;
    % z = corrcoef(Calinterval(xBP_L1')',Calinterval(xBP_R1')');
    % corrBP(i) = z(1,2) 
end

% sorting corrBP into 3 categories and get indices of trials for the 3 groups
[corrBPvalues,IndDescend]=sort(corrBP,'descend');
% Get the index of the EEG epochs according coorBPvalues
IndHigh=IndDescend(find(corrBPvalues > prctile(corrBPvalues,67)));
IndMedium=IndDescend(find(corrBPvalues <= prctile(corrBPvalues,67) & corrBPvalues >= prctile(corrBPvalues,33)));
IndLow=IndDescend(find(corrBPvalues < prctile(corrBPvalues,33)));

end

