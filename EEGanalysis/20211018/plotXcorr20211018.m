run organize_Behavioral_step1.m

cd /home/zhibin/Documents/GitHub/Motor_cordination/EEGanalysis/20211018
run organize_EEG_filter_step2.m

cd /home/zhibin/Documents/Acquisition/bimanual_Reach_zhibin_20211018
load filtered_data.mat

% filtered_data=filtered_broadband;
% filtered_data=filtered_delta;
% filtered_data=filtered_theta;
% filtered_data=filtered_alpha;
% filtered_data=filtered_beta;
% filtered_data=filtered_gamma;


run integrate_EEG_into_data_trials_step3.m
run Remove_EEG_artifact_Rereference_step6.m

%%
UniCondi=unique(CondiData);

for u=1:length(UniCondi)% u=2
    
    indtemp=find(CondiData==UniCondi(u));
    
    % pick a condition from randomized set allPerm
    conditionSelected = UniCondi(u);

    %subplot(1,length(UniCondi),u);
    subplot(7,13,abs(7-conditionSelected)*13+(14-conditionSelected));
    
    rAll=[];rAllMaxLag=[];
    
    yline(0);hold on;

    for indt=1:length(indtemp)

        C3EEG=reRef_data(:,28,indtemp(indt)); % Left hemisphere
        C4EEG=reRef_data(:,32,indtemp(indt)); % Right hemisphere

        [r,lags]=xcorr(C3EEG, C4EEG); % if r peak is positive, then the right hemisphere is leading
        rMaxLag=lags(find(r==max(r)));
        xline(rMaxLag);
        rAll=[rAll r];
        rAllMaxLag=[rAllMaxLag rMaxLag];
    end
    
    rAllmean=mean(rAll,2);
    %rPositiveMax=max(rAll,[],1);rNegativeMin=min(rAll,[],1);
    plot(lags,rAllmean);xlabel('time [ms]');ylabel('xcorr'); xlim([-100 100]); %ylim([-10000 10000]);% ylim([-1*1e4 1*1e4]);
    xline(0,'r');
    title(['condition: ' num2str(conditionNames(UniCondi(u))) ' (radian)']);
    hold off;
    
end

suptitle(['mean xcorr']);