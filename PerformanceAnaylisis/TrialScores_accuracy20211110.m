figure;
for conditionSelected=1:length(conditions)
    
    indtemp=find(CondiData==conditionSelected);
    
    subplot(1,length(conditions),length(conditions)+1-conditionSelected);
    BlockScores=TrialScores(indtemp);
    
    tbl=table((1:length(BlockScores))',BlockScores');mdl = fitlm(tbl,'linear');
    plot(mdl);ylim([0 1]);
    xlabel('trial');ylabel('acurracy');
    % title(['conditions : ' num2str(conditions(conditionSelected)) ' (radian)']);
    Rsquared=mdl.Rsquared.Ordinary;
    title({'condition ' conditionNames{UniCondi(conditionSelected)} 'Rsquared ' num2str(Rsquared)});
%     subplot(2,length(conditions),2*length(conditions)+1-conditionSelected);
%     BlockScores_predict=mdl.predict((1:length(BlockScores))')';
%     plot(1:length(BlockScores),BlockScores,'or',1:length(BlockScores),BlockScores_predict,'+b-');
end


%% plot different sessions with different color in one plot