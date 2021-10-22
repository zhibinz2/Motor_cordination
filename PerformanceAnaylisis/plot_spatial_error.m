cd('/home/zhibin/Documents/GitHub/Motor_cordination/behaviraldata')
load('matlab2021101802.mat');

% ans=BehavData{1};
% CondiData

NumTrialsRecorded=length(CondiData);

% % examine RSBehavData
% ans=RSBehavData{1};
% ans=RSBehavData{51};


% allPerm
% 
% 
% b=11;t=5;
% dataxLyLtemp=behaviraldata.dataBlock(b).dataTrialL(t).xLyL;
% dataxRyRtemp=behaviraldata.dataBlock(b).dataTrialR(t).xRyR;
% dataxJyJtemp=behaviraldata.dataBlock(b).dataTrialJ(t).xJyJ;

BehavData


for conditionSelected=1:length(conditions)
    
    indtemp=find(CondiData==conditionSelected);

    % produced the position parameters
    rad_ang=conditions(conditionSelected);
    [x,y] = drawReach(radius,rad_ang, xCenter, yCenter);
    
    subplot(7,15,abs(7-conditionSelected)*15+(15-conditionSelected-1));
    
    plot(x,y,'k');
    hold on;
    for i=1:length(indtemp)
        % plot(dataxJyJtemp(:,1),dataxJyJtemp(:,2),'g');
        % dataxJyJtemp=cell2mat(RSBehavData{indtemp(i)});
        dataxJyJtemp=BehavData{indtemp(i)};
        % AllxJ=dataxJyJtemp(5,:); ALLyJ=dataxJyJtemp(6,:);
        ALLxJ=dataxJyJtemp(:,5); ALLyJ=dataxJyJtemp(:,6);
        plotx(ALLxJ,ALLyJ,'g');
        ylim([100 600]);xlim([500 1500]);

        ALLd=[];
        for p=1:length(ALLxJ)
            v1=[x(1) y(1) 0];v2=[x(end) y(end) 0];pt=[ALLxJ(p) ALLyJ(p) 0];
            d = point_to_line(pt, v1, v2);
            ALLd=[ALLd d];
        end
        SpatialError(i)=sqrt(sumsqr(ALLd));
    end

    set(gca, 'YDir', 'reverse');
    hold off;
    % xlim([0 screenXpixels]);ylim([0 screenYpixels]);
    
    subplot(7,15,abs(7-conditionSelected)*15+(15-conditionSelected));

    plot(1:10,SpatialError,'go');hold on;
    tbl=table((1:10)',SpatialError');mdl = fitlm(tbl,'linear');
    plot(mdl);hold off;ylim([0 200]);
    ylabel('spatial error');xlabel('trial'); title('spatial error');
    %---delete the legend-------------
    s=findobj('type','legend');delete(s);
end

