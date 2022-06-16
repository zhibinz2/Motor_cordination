%% TriSegment_EEG
% cut x1 and x2 into 10 s per trials
x1=EEG_L;x2=EEG_R;
x1=zscore(x1);x2=zscore(x2);

win=10; % chunck of 10 seconds
SamplesLength=sr*win; % chunck of 10 seconds
nTrials=min([floor(length(x1)/SamplesLength) floor(length(x2)/SamplesLength)]);

xx1=[];xx2=[];
for i=1:nTrials
    xx1(:,i)=x1((SamplesLength*(i-1)+1):(SamplesLength*i),:);
    xx2(:,i)=x2((SamplesLength*(i-1)+1):(SamplesLength*i),:);
end

%% TriSegment BP same way as EEG
% win=10; % chunck of 10 seconds
% SamplesLength=sr*win; 
% nTrials=min([floor(length(BP_L)/SamplesLength) floor(length(BP_R)/SamplesLength)]);

corrBP=[]
for i=1:nTrials
    xBP_L1=[];xBP_R1=[];
    xBP_L1=BP_L((SamplesLength*(i-1)+1):(SamplesLength*i),:);
    xBP_R1=BP_R((SamplesLength*(i-1)+1):(SamplesLength*i),:);
    xBP_L2=[Calinterval(xBP_L1')'];
    xBP_R2=[Calinterval(xBP_R1')'];
    minTap=min([length(xBP_L2) length(xBP_R2)]);
    z = corr(xBP_L2(1:minTap), xBP_R2(1:minTap));
    corrBP(i)=z;
    % z = corrcoef(Calinterval(xBP_L1')',Calinterval(xBP_R1')');
    % corrBP(i) = z(1,2) 
end

%% sorting corrBP into 3 categories and get indices of trials for the 3 groups
[corrBPvalues,IndDescend]=sort(corrBP,'descend');
% Get the index of the EEG epochs according coorBPvalues
IndHigh=IndDescend(find(corrBPvalues > prctile(corrBPvalues,67)));
IndMedium=IndDescend(find(corrBPvalues <= prctile(corrBPvalues,67) & corrBPvalues >= prctile(corrBPvalues,33)));
IndLow=IndDescend(find(corrBPvalues < prctile(corrBPvalues,33)));

%% Select EEG
xx1(:,IndHigh);
xx1(:,IndMedium);
xx1(:,IndLow);

xx2(:,IndHigh);
xx2(:,IndMedium);
xx2(:,IndLow);
%% find the Botton presses within a EEG time interval
% % example
% x=Calinterval(BP_L'); % x is the intervals
% y=find(BP_L); % y is the indexes of BP
% y(1) % index of the first BP, which is 38
% z=[38; x] % index of each BP in EEG 
% z=[38 x] % index of each BP in EEG 
% zz=cumsum(x) 
% plot(zz)
% y = find(zz < 20000) % BP in the first 10 s
% y = find(zz < 40000 & zz > 20000) % BP in the second 10 s
% a = find(BP_L);
% y = find(a < 40000 & a > 20000)
% y = find(a < 40000 & a > 20000)'
% y = find(a < 60000 & a > 30000)'


%% time series of the corrcoef
% corrSeries=smoothing(corrSeries,win);
% corrSeries
% corrSeries=[];
% for i=1:(Maxlength-win+1)
% corrSeries(i)=corr(resam_BP_inter1(i:i+win-1),resam_BP_inter2(i:i+win-1));
% end
% corrSeries
% figure
% plot(corrSeries)
% hold on
% plot(corrSeries(7:13:end),'ro')
% close
% figure
% plot(corrSeries)
% hold on
% plot(7:13:541,corrSeries(7:13:end),'ro') % plot to examine the time intervals
