function [reRef_data] = reRef2(EEG);
% This function pass EEG (sample x chan) and get back average referenced data. 
averageRef=mean(EEG,2); % can we average over the good epochs too?
averageRef=repmat(averageRef,1,size(EEG,1));
reRef_data=EEG-averageRef; % this have all the channels, including VEOG,HEOG,M1,M2
end

