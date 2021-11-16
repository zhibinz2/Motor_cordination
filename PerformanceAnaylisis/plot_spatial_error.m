cd /home/zhibin/Documents/GitHub/Motor_cordination/behaviraldata
clear;
% load('matlab2021101802.mat');


cd /home/zhibin/Documents/GitHub/Motor_cordination/EEGanalysis/20211102
run organize_Behavioral_step1.m

% ans=BehavData{1};
% CondiData

NumTrialsRecorded=length(allPerm);
UniCondi=unique(allPerm);

% % examine behaviraldata
% ans=BehavData{1};
% ans=RSBehavData{51};


% allPerm
% 
% 
% b=11;t=5;
% dataxLyLtemp=behaviraldata.dataBlock(b).dataTrialL(t).xLyL;
% dataxRyRtemp=behaviraldata.dataBlock(b).dataTrialR(t).xRyR;
% dataxJyJtemp=behaviraldata.dataBlock(b).dataTrialJ(t).xJyJ;

% behaviraldata
% BehavData

for conditionSelected=1:length(UniCondi)
    
    indtemp=find(allPerm==conditionSelected);

    % produced the position parameters
    rad_ang=conditions(conditionSelected);
    [x,y] = drawReach(radius,rad_ang, xCenter, yCenter);
    
    % subplot(7,15,abs(7-conditionSelected)*15+(15-conditionSelected-1));
    subplot(2,length(UniCondi),length(UniCondi)+1-conditionSelected);
    
    plot(x,y,'k');
    hold on;
    for i=1:length(indtemp)
        % plot(dataxJyJtemp(:,1),dataxJyJtemp(:,2),'g');
        % dataxJyJtemp=cell2mat(RSBehavData{indtemp(i)});
        dataxJyJtemp=BehavData{indtemp(i)};
        % AllxJ=dataxJyJtemp(5,:); ALLyJ=dataxJyJtemp(6,:);
        ALLxJ=dataxJyJtemp(:,5); ALLyJ=dataxJyJtemp(:,6);
        plotx(ALLxJ,ALLyJ,'g');
        ylim([100 600]);xlim([500 1500]);title(conditionNames(length(UniCondi)));

        ALLd=[];
        for p=1:length(ALLxJ)
            v1=[x(1) y(1) 0];v2=[x(end) y(end) 0];pt=[ALLxJ(p) ALLyJ(p) 0];
            d = point_to_line(pt, v1, v2);
            ALLd=[ALLd d];
        end
        SpatialError(i)=sqrt((sumsqr(ALLd))/length(ALLd))/TrialScores(indtemp(i));
    end

    set(gca, 'YDir', 'reverse');
    hold off;
    % xlim([0 screenXpixels]);ylim([0 screenYpixels]);
    
    % subplot(7,15,abs(7-conditionSelected)*15+(15-conditionSelected));
    subplot(2,length(UniCondi),2*length(UniCondi)+1-conditionSelected);
    plot(1:length(indtemp),SpatialError,'go');hold on;
    tbl=table((1:length(indtemp))',SpatialError');mdl = fitlm(tbl,'linear');
    plot(mdl);hold off;ylim([0 120]);
    ylabel('spatial error (pixels/score)');xlabel('trial'); title(['r^2 = ' num2str(mdl.Rsquared.Ordinary)]);
    %---delete the legend-------------
    % s=findobj('type','legend');delete(s);
end

% suptitle('spatial error');
