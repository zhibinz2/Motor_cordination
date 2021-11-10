%% Prefrontal
% Left prefrontal
Leftticks=[1 4 6 7 8 9 10 72 73 74];
% Right prefrontal
Rightticks=[3 5 12 13 14 15 16 75 76 77];

%% plot just one condtion and only the motor chans
% Left hemisphere
Leftticks=[18 19 20 27 28 29 93 94 95];
% Right hemisphere
Rightticks=[22 23 24 31 32 33 96 97 98];

%% plot just one condtion and only the parietal chans

% Left hemisphere
Leftticks=[36 37 38 46 47 48 49 109 110 111];
% Right hemisphere
Rightticks=[40 41 42 51 52 53 54 112 113 114];


%%
u=3; % u=7; % plot just one condtion
    
indtemp=find(CondiData==UniCondi(u));

%%
uall=[7 3];

for sp=1:2
    u=uall(sp);
    indtemp=find(CondiData==UniCondi(u));
    subplot(1,2,sp);


    [out,idx] = sort(TrialScores(indtemp)); % example: [out,idx] = sort([14 8 91 19])

    % Average every n elements
    n=10;
    % OutMean = arrayfun(@(i) mean(out(i:i+n-1)),1:n:length(out)-n+1);
    OutMean = arrayfun(@(i) mean(out(i:i+n-1)),1:1:length(out)-n+1);
    OutStd = arrayfun(@(i) std(out(i:i+n-1)),1:1:length(out)-n+1);


    CohMean=zeros(1,length(out)-n+1);

    for i=1:1:length(out)-n+1
        % First, plot the scalp coherence
        [pow,freqs,df,eppow,corr,cprod,fcoef] = allspectra(reRef_data,rate,maxfreq,idx(i:i+n-1),win);
        % calculate the coherence
        coh=abs(corr).^2;
        % Plot coh of specific frequency between the two hemispheres
        CohMean(i)=mean(squeeze(coh(46,Leftticks,Rightticks)),'all');
    end

    yyaxis left; 
    % errorbar(OutMean,OutStd);
    plot(1:length(out)-n+1,OutMean);
    ylabel('Average accuracy');hold on;
    yyaxis right; plot(1:length(out)-n+1,CohMean);ylabel('Average coherence');xlabel('per 10 trials');
    xlim([1 length(out)-n+1]);
    title({'45 Hz', 'coh - condition: ' num2str(conditions(u))});hold off;

end
    
