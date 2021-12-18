% tryout the calculation
t = 1:10; %time
h = t.^2; %distance
dt = diff(t) 
dh = diff(h)
% v - velocity
v = dh./dt
dv = diff(v)
% a - acceleration
a = dv./dt(2:end)
plot(t,h,'-ob');
hold on;
plot(t(2:end),v,'-og');
plot(t(3:end),a,'-or');
legend({'h', 'v','a'})

%%
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
%     rad_ang=conditions(conditionSelected);
%     [x,y] = drawReach(radius,rad_ang, xCenter, yCenter);
    
%     plot(x,y,'k');
%     hold on;

    ALLhL=[];ALLhR=[];ALLvL=[];ALLvR=[];
    for i=1:length(indtemp) %loop through every trial
        % plot(dataxJyJtemp(:,1),dataxJyJtemp(:,2),'g');
        % dataxLRJyLRJtemp=cell2mat(RSBehavData{indtemp(i)});
        dataxLRJyLRJtemp=BehavData{indtemp(i)};
        
        ALLxL=dataxLRJyLRJtemp(:,1); ALLyL=dataxLRJyLRJtemp(:,2);
        ALLxR=dataxLRJyLRJtemp(:,3); ALLyR=dataxLRJyLRJtemp(:,4);
        % AllxJ=dataxLRJyLRJtemp(5,:); ALLyJ=dataxLRJyLRJtemp(6,:);
%         ALLxJ=dataxLRJyLRJtemp(:,5); ALLyJ=dataxLRJyLRJtemp(:,6);

%         plotx(ALLxL,ALLyL,'r');hold on
%         plotx(ALLxR,ALLyR,'b');
%         plotx(ALLxJ,ALLyJ,'g');hold off
%         ylim([100 600]);xlim([500 1500]);
%         xlim([0 screenXpixels]);ylim([0 screenYpixels]);
        
        t=1:length(ALLxL);
        hL=xCenter-ALLxL';hR=ALLxR'-xCenter;
        ALLhL=[ALLhL; hL]; ALLhR=[ALLhR; hR];
        dt=diff(t);
        dhL=diff(hL); dhR=diff(hR);
        
        % v - velocity
        vL = dhL./dt; vR = dhR./dt;
        ALLvL=[ALLvL; vL]; ALLvR=[ALLvR; vR];
        dvL = diff(vL); dvR = diff(vR);
        % a - acceleration
        aL = dvL./dt(2:end); aR = dvR./dt(2:end);
        
        
%         subplot(7,15,abs(7-conditionSelected)*15+(15-conditionSelected-1));
%         
%         % plot the distance
%         plot(t,hL,'.r');hold on; plot(t,hR,'.b');legend({'hL','hR'});
%         xlabel('t');ylabel('distance (pixels)');
%         %---delete the legend-------------
%         s=findobj('type','legend');delete(s);

%         % plot the velocity
%         plot(t(2:end),vL,'-r');hold on; plot(t(2:end),vR,'-b');hold off;
%         legend({'vL','vR');
        
%         % plot the acceleration
%         plot(t(3:end),aL,'or');hold on; plot(t(3:end),aR,'ob');hold off;
%         legend({'aL','aR'});
        
%         hold off;
%         legend({'hL','hR', 'vL','vR','aL','aR'});
        

%         ALLd=[];
%         for p=1:length(ALLxJ)
%             v1=[x(1) y(1) 0];v2=[x(end) y(end) 0];pt=[ALLxJ(p) ALLyJ(p) 0];
%             d = point_to_line(pt, v1, v2);
%             ALLd=[ALLd d];
%         end
%         SpatialError(i)=sqrt(sumsqr(ALLd));
%     end

%     set(gca, 'YDir', 'reverse');
%     hold off;
    % xlim([0 screenXpixels]);ylim([0 screenYpixels]);
    
%     subplot(7,15,abs(7-conditionSelected)*15+(15-conditionSelected));
%     
%     % plot the velocity
%     plot(t(2:end),vL,'-r');hold on; plot(t(2:end),vR,'-b');
%     legend({'vL','vR'});xlabel('t');ylabel('velocity (pixels/frame)');

%     plot(1:10,SpatialError,'go');hold on;
%     tbl=table((1:10)',SpatialError');mdl = fitlm(tbl,'linear');
%     plot(mdl);hold off;ylim([0 200]);
%     ylabel('spatial error');xlabel('trial'); title('spatial error');
%     %---delete the legend-------------
%     s=findobj('type','legend');delete(s);
    end
    
%     subplot(7,15,abs(7-conditionSelected)*15+(15-conditionSelected-1));
    subplot(2,7,7+1-conditionSelected);

    % plot the distance
    averagehL=mean(ALLhL,1);averagehR=mean(ALLhR,1);
    plot(t,averagehL,'.r');hold on; plot(t,averagehR,'.b');hold off;
    xlabel('t');ylabel('average distance (pixels)');
    title(['condition ' conditionNames{UniCondi(conditionSelected)}]);
    legend({'Left','Right'});

    %subplot(7,15,abs(7-conditionSelected)*15+(15-conditionSelected));
    subplot(2,7,14+1-conditionSelected);
    
    % plot the velocity
    averagevL=mean(ALLvL,1);averagevR=mean(ALLvR,1);
    plot(t(2:end),averagevL,'-r');hold on; plot(t(2:end),averagevR,'-b');hold off;
    xlabel('t');ylabel('average velocity (pixels/frame)');legend({'Left','Right'});
    
    %---delete the legend-------------


%     hold off;
end
