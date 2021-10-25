for conditionSelected=1:length(conditions)
    
    indtemp=find(CondiData==conditionSelected);
    
    subplot(4,4,conditionSelected);
    BlockScores=TrialScores(indtemp);
    
    tbl=table((1:30)',BlockScores');mdl = fitlm(tbl,'linear');
    plot(mdl);ylim([0 1]);
    
    xlabel('trial');ylabel('acurracy');
    title(['conditions : ' num2str(conditions(conditionSelected)) ' (radian)']);
end
