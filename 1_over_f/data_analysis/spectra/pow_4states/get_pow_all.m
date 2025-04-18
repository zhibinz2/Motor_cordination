sr=2000; maxfreq=50;
% get the EEG pow from L and R
% organize pow from all sessions into a matrix 
tic
pow_all=nan(12,2,12,51,32); % session x L/R X trials x pow [51 freq x 32 chan]
coh_all=nan(12,2,12,51,32,32); % session x L/R X trials x pow [51 freq x 32 chan x 32 chan]
ampcorr_all=nan(12,2,12,51,32,32); % session x L/R X trials x pow [51 freq x 32 chan x 32 chan]
for s=1:numSes
    clear data
    runid=num2str(seeds(s,:));
    data=load(['/ssd/zhibin/1overf/' runid '_2P/Cleaned_data/clean_' runid '.mat'],'dataL','dataR');
    for sub=1:2 % R/L
        for t=1:12 % 12 trials
        clear epoches;
        if sub==1; 
        % cut data into epoches of 1 sec
        [epoches] = intoepoches(data.dataL{t}(:,1:32),sr,1);
        else
        [epoches] = intoepoches(data.dataR{t}(:,1:32),sr,1);
        end
        % compute pow
        [pow_all(s,sub,t,:,:),~,~,coh_all(s,sub,t,:,:,:),~,ampcorr_all(s,sub,t,:,:,:),~] ...
            = allspectra1(epoches,sr,maxfreq);
        end
    end
end
toc % Elapsed time is 158.317634 seconds.

%% check output size
% size(pow_all,1:5)
% size(coh_all,1:6)
% size(ampcorr_all,1:6)