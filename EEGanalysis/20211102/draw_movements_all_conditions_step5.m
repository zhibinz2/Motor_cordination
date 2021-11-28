%% draw condition trace

% plot the mouse trace
UniCondi=unique(CondiData);

figure;
for u=1:length(UniCondi)% u=2
    
    indtemp=find(CondiData==UniCondi(u));
    
    % pick a condition from randomized set allPerm
    conditionSelected = UniCondi(u);
    % produced the position parameters
        rad_ang=conditions(conditionSelected);
        [x,y] = drawReach(radius,rad_ang, xCenter, yCenter);
            
    
    subplot(1,length(UniCondi),length(UniCondi)+1-u);
    % subplot(7,13,abs(7-conditionSelected)*13+(14-conditionSelected));
    
    plot(x,y,'-', 'color', grey, 'LineWidth',5);
    hold on;

%     for indt=1:length(indtemp)
%     TracexL=squeeze(data_trials(:,NumEEGChannels+1,indtemp(indt)));
%     TraceyL=squeeze(data_trials(:,NumEEGChannels+2,indtemp(indt)));
%     TracexR=squeeze(data_trials(:,NumEEGChannels+3,indtemp(indt)));
%     TraceyR=squeeze(data_trials(:,NumEEGChannels+4,indtemp(indt)));
%     TracexJ=squeeze(data_trials(:,NumEEGChannels+5,indtemp(indt)));
%     TraceyJ=squeeze(data_trials(:,NumEEGChannels+6,indtemp(indt)));
%     plot(TracexL,TraceyL,'-r',TracexR,TraceyR,'-b',TracexJ,TraceyJ,'-g');
%     end
    
    % Without using EEG data or resampling is better
    for indt=1:length(indtemp)
        dataxLRJyLRJtemp=BehavData{indtemp(indt)};
        TracexL=dataxLRJyLRJtemp(:,1);
        TraceyL=dataxLRJyLRJtemp(:,2);
        TracexR=dataxLRJyLRJtemp(:,3);
        TraceyR=dataxLRJyLRJtemp(:,4);
        TracexJ=dataxLRJyLRJtemp(:,5);
        TraceyJ=dataxLRJyLRJtemp(:,6);
    plot(TracexL,TraceyL,'-r',TracexR,TraceyR,'-b',TracexJ,TraceyJ,'-g');
    end
    
    
    set(gca, 'YDir', 'reverse');
    hold off;
    
    xlim([0 screenXpixels]);ylim([0 screenYpixels]);
    xlabel('Screen-x-pixel');ylabel('Screen-y-pixel');
    title(['condition ' conditionNames{UniCondi(u)}]);
    %title(['condition: ' num2str(conditionNames(UniCondi(u))) ' (radian)']);
  
end

