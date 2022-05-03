stimTrials_synchronization_r,
stimTrials_syncopation_r,
stimTrials_rt_r,
tapTrials_synchronization_r,
tapTrials_syncopation_r,
tapTrials_rt_r

%%
tapTrials_synchronization_erp=mean(tapTrials_synchronization_r(:,4:32,:),3);
[U,S,V]=svd(tapTrials_synchronization_erp,"econ");
figure;image(S);
image(U);image(V);
diag(S) % the eigen values (a set of weights) for the 32 channels 

% See if U(:,1) looks like the readiness potential?
plot(U(:,1));xlabel('time');title('U(:,1)');xline(2000,'m',{'reponse'});
plot(-1*U(:,1));xlabel('time');title('-1*U(:,1)');xline(2000,'m',{'reponse'});
% topoplot to examine them
cd /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/channels_info
load('chaninfo.mat')
chaninfo29=chaninfo(4:32);
figure;
topoplot(V(:,1),chaninfo,'nosedir','+X');
topoplot(V(:,1),chaninfo29,'nosedir','+X');title('+1*V(:,1)');
topoplot(-1*V(:,1),chaninfo29,'nosedir','+X');title('-1*V(:,1)');

% stimTrials_synchronization_erp=mean(stimTrials_synchronization_r,3);
stimTrials_synchronization_erp_estimate_trial=zeros(size(stimTrials_synchronization_r,1),1);
for i=1:240
    stimTrials_synchronization_erp_estimate_trial(:,i)=squeeze(stimTrials_synchronization_r(:,4:32,i))*(-1*V(:,1));
end
    
figure;
for i=1:16
    subplot(4,4,i)
    plot(stimTrials_synchronization_erp_estimate_trial(:,i));
    title(['trial' num2str(i)]); xline(2000,'m',{'task'});
end
suptitle('single trial estimate of ERP (multiply with V(:,1), look for negative peaks)');