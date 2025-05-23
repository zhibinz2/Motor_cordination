function [deltaEEG,thetaEEG,alphaEEG,betaEEG,gammaEEG] = filter5band(EEG,sr)
% This function use makefilter to filter each block of EEG into 5 bands
% my band class
% delta 1-3, theta 4-7, alpha 8-12, beta 13-30, gamma 31-50
%% add padding
padding=zeros(round(size(EEG,1)/10), size(EEG,2));
filtered_data1=cat(1,padding,EEG,padding);
%% delta
filtered_data2=[];
% low pass 
Hd = makefilter(sr,3,3.5,6,20,0);  
filtered_data2=filtfilthd(Hd,filtered_data1);
% remove padding
deltaEEG=filtered_data2((size(padding,1)+1):(size(padding,1)+size(EEG,1)),:);
% plot(deltaEEG)

%% theta
filtered_data2=[]; filtered_data3=[];
% high pass 
Hd = makefilter(sr,4,3.5,6,20,0); 
filtered_data2=filtfilthd(Hd,filtered_data1);
% low pass 
Hd = makefilter(sr,7,7.5,6,20,0);  
filtered_data3=filtfilthd(Hd,filtered_data2);
% remove padding
thetaEEG=filtered_data3((size(padding,1)+1):(size(padding,1)+size(EEG,1)),:);
% plot(thetaEEG)

% alpha
filtered_data2=[]; filtered_data3=[];
% high pass 
Hd = makefilter(sr,8,7.5,6,20,0); 
filtered_data2=filtfilthd(Hd,filtered_data1);
% low pass 
Hd = makefilter(sr,12,12.5,6,20,0);  
filtered_data3=filtfilthd(Hd,filtered_data2);
% remove padding
alphaEEG=filtered_data3((size(padding,1)+1):(size(padding,1)+size(EEG,1)),:);
% plot(alphaEEG)

% beta
filtered_data2=[]; filtered_data3=[];
% high pass 
Hd = makefilter(sr,13,12.5,6,20,0); % for keeping readiness potential
filtered_data2=filtfilthd(Hd,filtered_data1);
% low pass 
Hd = makefilter(sr,30,30.5,6,20,0);  
filtered_data3=filtfilthd(Hd,filtered_data2);
% remove padding
betaEEG=filtered_data3((size(padding,1)+1):(size(padding,1)+size(EEG,1)),:);
% plot(betaEEG)

% gamma
filtered_data2=[];
% high pass 
Hd = makefilter(sr,31,30.5,6,20,0); % for keeping readiness potential
filtered_data2=filtfilthd(Hd,filtered_data1);
% % low pass 
% Hd = makefilter(sr,49.5,50,6,20,0);  
% filtered_data3=filtfilthd(Hd,filtered_data2);
% remove padding
gammaEEG=filtered_data2((size(padding,1)+1):(size(padding,1)+size(EEG,1)),:);
% plot(gammaEEG)
end