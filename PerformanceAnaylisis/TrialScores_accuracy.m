for conditionSelected=1:length(conditions)
    
    indtemp=find(CondiData==conditionSelected);
    
    subplot(1,5,6-conditionSelected);
    BlockScores=TrialScores(indtemp);
    
    tbl=table((1:length(BlockScores))',BlockScores');mdl = fitlm(tbl,'linear');
    plot(mdl);ylim([0 1]);
    
    xlabel('trial');ylabel('acurracy');
    % title(['conditions : ' num2str(conditions(conditionSelected)) ' (radian)']);
    title(['condition ' conditionNames{UniCondi(conditionSelected)}]);
end
