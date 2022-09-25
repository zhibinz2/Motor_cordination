function [deltaEEG,thetaEEG,alphaEEG,betaEEG,gammaEEG] = filter5band(EEG,sr)
% This function use makefilter to filter each block of EEG into 5 bands
% my band class
% delta 1-3, theta 4-7, alpha 8-12, beta 13-30, gamma 31-50

% delta
clear filtered_data1 filtered_data2 filtered_data3
% % high pass (no paddings needed)
% Hd = makefilter(sr,0.25,0.01,6,20,0); % for keeping readiness potential
% filtered_data1=filtfilthd(Hd,EEG);
% add padding (this step takes several minutes)
padding=zeros(round(size(EEG,1)/10), size(EEG,2));
filtered_data2=cat(1,padding,EEG,padding);
% low pass 
Hd = makefilter(sr,3,3.5,6,20,0);  
filtered_data3=filtfilthd(Hd,filtered_data2);
% remove padding
deltaEEG=filtered_data3((size(padding,1)+1):(size(padding,1)+size(EEG,1)),:);

% theta
clear filtered_data1 filtered_data2 filtered_data3
% high pass (no paddings needed)
Hd = makefilter(sr,4,3.5,6,20,0); % for keeping readiness potential
filtered_data1=filtfilthd(Hd,EEG);
% add padding (this step takes several minutes)
padding=zeros(round(size(filtered_data1,1)/10), size(filtered_data1,2));
filtered_data2=cat(1,padding,filtered_data1,padding);
% low pass 
Hd = makefilter(sr,7,7.5,6,20,0);  
filtered_data3=filtfilthd(Hd,filtered_data2);
% remove padding
thetaEEG=filtered_data3((size(padding,1)+1):(size(padding,1)+size(EEG,1)),:);


% alpha
clear filtered_data1 filtered_data2 filtered_data3
% high pass (no paddings needed)
Hd = makefilter(sr,8,7.5,6,20,0); % for keeping readiness potential
filtered_data1=filtfilthd(Hd,EEG);
% add padding (this step takes several minutes)
padding=zeros(round(size(filtered_data1,1)/10), size(filtered_data1,2));
filtered_data2=cat(1,padding,filtered_data1,padding);
% low pass 
Hd = makefilter(sr,12,12.5,6,20,0);  
filtered_data3=filtfilthd(Hd,filtered_data2);
% remove padding
alphaEEG=filtered_data3((size(padding,1)+1):(size(padding,1)+size(EEG,1)),:);

% beta
clear filtered_data1 filtered_data2 filtered_data3
% high pass (no paddings needed)
Hd = makefilter(sr,13,12.5,6,20,0); % for keeping readiness potential
filtered_data1=filtfilthd(Hd,EEG);
% add padding (this step takes several minutes)
padding=zeros(round(size(filtered_data1,1)/10), size(filtered_data1,2));
filtered_data2=cat(1,padding,filtered_data1,padding);
% low pass 
Hd = makefilter(sr,30,30.5,6,20,0);  
filtered_data3=filtfilthd(Hd,filtered_data2);
% remove padding
betaEEG=filtered_data3((size(padding,1)+1):(size(padding,1)+size(EEG,1)),:);

% gamma
clear filtered_data1 filtered_data2 filtered_data3
% high pass (no paddings needed)
Hd = makefilter(sr,31,30.5,6,20,0); % for keeping readiness potential
gammaEEG=filtfilthd(Hd,EEG);
% % add padding (this step takes several minutes)
% padding=zeros(round(size(filtered_data1,1)/10), size(filtered_data1,2));
% filtered_data2=cat(1,padding,filtered_data1,padding);
% % low pass 
% Hd = makefilter(sr,49.5,50,6,20,0);  
% filtered_data3=filtfilthd(Hd,filtered_data2);
% % remove padding
% betaEEG=filtered_data3((size(padding,1)+1):(size(padding,1)+size(EEG,1)),:);

end