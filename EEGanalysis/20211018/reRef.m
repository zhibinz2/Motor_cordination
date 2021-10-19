function [reRef_data] = reRef(filtered_data,goodchans);
% This function pass filtered_data and get back average referenced data. 
averageRef=mean(filtered_data(:,goodchans,:),2); % can we average over the good epochs too?
averageRef=repmat(averageRef,1,size(filtered_data,2));
reRef_data=filtered_data-averageRef; % this have all the channels, including VEOG,HEOG,M1,M2
end

