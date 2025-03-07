% reorganzie pow_all according to the 4 states
condi_pow_all=cell(1,4); % 4 state {51 freq x 32 chan x 72 trials}
condi_coh_all=cell(1,4); % 4 state {51 freq x 32 chan x 32 chan x 72 trials}
condi_ampcorr_all=cell(1,4); % 4 state {51 freq x 32 chan x 32 chan x 72 trials}
tic
for s=1:numSes
        for t=1:12 % 12 trials
            if condition_all(s,t)==1; % Uncouple
                condi_pow_all{1}=cat(3,condi_pow_all{1}, squeeze(pow_all(s,1,t,:,:)));
                condi_pow_all{1}=cat(3,condi_pow_all{1}, squeeze(pow_all(s,2,t,:,:)));
                condi_coh_all{1}=cat(4,condi_coh_all{1}, squeeze(coh_all(s,1,t,:,:,:)));
                condi_coh_all{1}=cat(4,condi_coh_all{1}, squeeze(coh_all(s,2,t,:,:,:)));
                condi_ampcorr_all{1}=cat(4,condi_ampcorr_all{1}, squeeze(ampcorr_all(s,1,t,:,:,:)));
                condi_ampcorr_all{1}=cat(4,condi_ampcorr_all{1}, squeeze(ampcorr_all(s,2,t,:,:,:)));
            elseif condition_all(s,t)==2; % L-lead
                condi_pow_all{2}=cat(3,condi_pow_all{2}, squeeze(pow_all(s,1,t,:,:)));% L leading
                condi_pow_all{3}=cat(3,condi_pow_all{3}, squeeze(pow_all(s,2,t,:,:)));% R following
                condi_coh_all{2}=cat(4,condi_coh_all{2}, squeeze(coh_all(s,1,t,:,:,:)));
                condi_coh_all{3}=cat(4,condi_coh_all{3}, squeeze(coh_all(s,2,t,:,:,:)));
                condi_ampcorr_all{2}=cat(4,condi_ampcorr_all{2}, squeeze(ampcorr_all(s,1,t,:,:,:)));
                condi_ampcorr_all{3}=cat(4,condi_ampcorr_all{3}, squeeze(ampcorr_all(s,2,t,:,:,:)));
            elseif condition_all(s,t)==3; % R-lead
                condi_pow_all{2}=cat(3,condi_pow_all{2}, squeeze(pow_all(s,2,t,:,:)));% R leading
                condi_pow_all{3}=cat(3,condi_pow_all{3}, squeeze(pow_all(s,1,t,:,:)));% L following
                condi_coh_all{2}=cat(4,condi_coh_all{2}, squeeze(coh_all(s,2,t,:,:,:)));
                condi_coh_all{3}=cat(4,condi_coh_all{3}, squeeze(coh_all(s,1,t,:,:,:)));
                condi_ampcorr_all{2}=cat(4,condi_ampcorr_all{2}, squeeze(ampcorr_all(s,2,t,:,:,:)));
                condi_ampcorr_all{3}=cat(4,condi_ampcorr_all{3}, squeeze(ampcorr_all(s,1,t,:,:,:)));
            else
                condi_pow_all{4}=cat(3,condi_pow_all{4}, squeeze(pow_all(s,1,t,:,:)));
                condi_pow_all{4}=cat(3,condi_pow_all{4}, squeeze(pow_all(s,2,t,:,:)));
                condi_coh_all{4}=cat(4,condi_coh_all{4}, squeeze(coh_all(s,1,t,:,:,:)));
                condi_coh_all{4}=cat(4,condi_coh_all{4}, squeeze(coh_all(s,2,t,:,:,:)));
                condi_ampcorr_all{4}=cat(4,condi_ampcorr_all{4}, squeeze(ampcorr_all(s,1,t,:,:,:)));
                condi_ampcorr_all{4}=cat(4,condi_ampcorr_all{4}, squeeze(ampcorr_all(s,2,t,:,:,:)));
            end
        end
end
toc % Elapsed time is 5.812487 seconds.

% % Z-score across 72 values of power separately in each of the 4 states
% for st=1:4
%     condi_pow_all{st}=zscore(condi_pow_all{st},0,3);
% end

% Z-score across 4x72 values of power in the 4 states
grand_mean=mean(cat(3,condi_pow_all{1:4}),3);
grand_std=std(cat(3,condi_pow_all{1:4}),0,3);
% same as cat(3,condi_pow_all{1},condi_pow_all{2},condi_pow_all{3},condi_pow_all{4});
for st=1:4
    condi_pow_all{st}=(condi_pow_all{st}-repmat(grand_mean,1,1,72))./repmat(grand_std,1,1,72);
end

% get average pow for each state
avg_state_pow=nan(4,51,32); % 4 state x 51 freq x 32 chan
% get average corr and ampcorr for each state
avg_state_coh=nan(4,51,32,32); % 4 state x 51 freq x 32 chan x 32 chan
avg_state_ampcorr=nan(4,51,32,32); % 4 state x 51 freq x 32 chan x 32 chan
for st=1:4
    for freq=1:51
    avg_state_pow(st,freq,:)=mean(condi_pow_all{st}(freq,:,:),3);
    avg_state_coh(st,freq,:,:)=mean(condi_coh_all{st}(freq,:,:,:),4);
    avg_state_ampcorr(st,freq,:,:)=mean(condi_ampcorr_all{st}(freq,:,:,:),4);
    end
end

%% check size of outputs
% size(avg_state_coh,1:4)
% size(avg_state_ampcorr,1:4)

%% tryout
% ss=4; freq_s=4; chan_s=5;
% oneseries=squeeze(condi_pow_all{ss}(freq_s,chan_s,:))
% mean_one=mean(oneseries)
% std_one=std(oneseries)
% figure;
% subplot(1,3,1); plot(oneseries); title('72 pow values')
% subplot(1,3,2); plot(zscore(oneseries)); title('after zscoring');
% hold on; yline(mean(zscore(oneseries)),'m-','mean');
% yline(std(zscore(oneseries)),'b-','std'); hold off;
% subplot(1,3,3); plot((oneseries-mean_one)/std_one);title('manual zscore')
% sgtitle(['condition ' states4names{ss} ';  freq ' num2str(freq_s) ' Hz;  chan ' num2str(chan_s)]);
