% tryout the calculation
%     t = 1:10; %time
%     h = t.^2; %distance
%     dt = diff(t) 
%     dh = diff(h)
%     % v - velocity
%     v = dh./dt
%     dv = diff(v)
%     % a - acceleration
%     a = dv./dt(2:end)
%     plot(t,h,'-ob');
%     hold on;
%     plot(t(2:end),v,'-og');
%     plot(t(3:end),a,'-or');
%     legend({'h', 'v','a'})

%% Using the original behaviral data (averaged)
% cd('/home/zhibin/Documents/GitHub/Motor_cordination/behaviraldata')
% load('matlab2021101802.mat');

% original behaviral data
% ans=BehavData{1};
% CondiData

% NumTrialsRecorded=length(CondiData);

% % examine RSBehavData
% ans=RSBehavData{1};
% ans=RSBehavData{51};


% allPerm

% b=11;t=5;
% dataxLyLtemp=behaviraldata.dataBlock(b).dataTrialL(t).xLyL;
% dataxRyRtemp=behaviraldata.dataBlock(b).dataTrialR(t).xRyR;
% dataxJyJtemp=behaviraldata.dataBlock(b).dataTrialJ(t).xJyJ;

BehavData;


figure('units','normalized','outerposition',[0 0 1 1]);

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
    legend({'Left','Right'},'Location','northwest');
    ylim([0 1000]);
    
    %subplot(7,15,abs(7-conditionSelected)*15+(15-conditionSelected));
    subplot(2,7,14+1-conditionSelected);
    
    % plot the velocity
    averagevL=mean(ALLvL,1);averagevR=mean(ALLvR,1);
    plot(t(2:end),averagevL,'-r');hold on; plot(t(2:end),averagevR,'-b');hold off;
    xlabel('t');ylabel('average velocity (pixels/frame)');legend({'Left','Right'},'Location','northwest');
    ylim([-5 45]);
    
    %---delete the legend-------------


%     hold off;
end


%% Using the resampled behaviral data (all trials)
% cd('/home/zhibin/Documents/GitHub/Motor_cordination/behaviraldata')
% load('matlab2021101802.mat');

% Original behaviral data
% ans=BehavData{1};
% CondiData

% NumTrialsRecorded=length(CondiData);

% % examine RSBehavData (the resampled one)
% data_trials;
% ans=cell2mat(RSBehavData{1});

close;
figure('units','normalized','outerposition',[0 0 1 1]);

for conditionSelected=1:length(conditions)
    
    indtemp=find(CondiData==conditionSelected);

    ALLhL=[];ALLhR=[];ALLvL=[];ALLvR=[];ALLaL=[];ALLaR=[];
    for i=1:length(indtemp) %loop through every trial
        % plot(dataxJyJtemp(:,1),dataxJyJtemp(:,2),'g');
        dataxLRJyLRJtemp=cell2mat(RSBehavData{indtemp(i)});
        dataxLRJyLRJtemp=dataxLRJyLRJtemp';
        
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
        
        % distance
        hL=xCenter-ALLxL';hR=ALLxR'-xCenter;
        % all trials
        ALLhL=[ALLhL; hL]; ALLhR=[ALLhR; hR];
        
        % change in time
        dt=diff(t);
        % change in distance
        dhL=diff(hL); dhR=diff(hR);
        
        % v - velocity
        vL = dhL./dt; vR = dhR./dt;
        % all trials
        ALLvL=[ALLvL; vL]; ALLvR=[ALLvR; vR];
        
        % change in velocity
        dvL = diff(vL); dvR = diff(vR); 
        % a - acceleration
        aL = dvL./dt(2:end); aR = dvR./dt(2:end);
        % all trials
        ALLaL =[ALLaL; aL]; ALLaR=[ALLaR; aR];
        
        
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
    

    subplot(3,7,7+1-conditionSelected);
    % plot the distance
    averagehL=mean(ALLhL,1);averagehR=mean(ALLhR,1);
    plot(t,averagehL,'.r');hold on; plot(t,averagehR,'.b');hold off;
    xlabel('t');ylabel('average distance (pixels)');
    title(['condition ' conditionNames{UniCondi(conditionSelected)}]);
    legend({'Left','Right'},'Location','northwest');
    ylim([0 1000]);
    
    subplot(3,7,14+1-conditionSelected);
    % plot the velocity
    averagevL=mean(ALLvL,1);averagevR=mean(ALLvR,1);
    plot(t(2:end),averagevL,'-r');hold on; plot(t(2:end),averagevR,'-b');hold off;
    xlabel('t');ylabel('average velocity (pixels/1000ms)');legend({'Left','Right'},'Location','northwest');
    title('velocity');
    ylim([-0.2 2]);
    
    
    subplot(3,7,21+1-conditionSelected);
    % plot the acceleration
    averageaL=mean(ALLaL,1);averageaR=mean(ALLaR,1);
    if conditionSelected > 4
    plot(t(3:end),averageaL,'-r');hold on; 
    plot(t(3:end),averageaR,'-b');hold off;
    elseif conditionSelected < 5
    plot(t(3:end),averageaR,'-b'); hold on; 
    plot(t(3:end),averageaL,'-r');hold off;
    end
    xlabel('t');ylabel('average acceleration (pixels/1000ms^2)');legend({'Left','Right'},'Location','southwest');
    title('acceleration');
    % ylim([-0.2 2]);
    %---delete the legend-------------


%     hold off;
end


%% Using the resampled behaviral data and apply Savitzky-Golay Filter (average)        
close;
figure('units','normalized','outerposition',[0 0 1 1]);

for conditionSelected=1:length(conditions)
    
    indtemp=find(CondiData==conditionSelected);

    ALLhL=[];ALLhR=[];ALLvL=[];ALLvR=[];ALLaL=[];ALLaR=[];
    for i=1:length(indtemp) %loop through every trial
        % plot(dataxJyJtemp(:,1),dataxJyJtemp(:,2),'g');
        dataxLRJyLRJtemp=cell2mat(RSBehavData{indtemp(i)});
        dataxLRJyLRJtemp=dataxLRJyLRJtemp';
        
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
        
        % distance
        hL=xCenter-ALLxL';hR=ALLxR'-xCenter;
        % all trials
        ALLhL=[ALLhL; hL]; ALLhR=[ALLhR; hR];
        
        % change in time
        dt=diff(t);
        % change in distance
        dhL=diff(hL); dhR=diff(hR);
        
        % v - velocity
        vL = dhL./dt; vR = dhR./dt;
        % Savitzky-Golay Filter 
        vL = sgolayfilt(vL,1,251); vR = sgolayfilt(vR,1,251); 
        % all trials
        ALLvL=[ALLvL; vL]; ALLvR=[ALLvR; vR];
        
        % change in velocity
        dvL = diff(vL); dvR = diff(vR); 
        % a - acceleration
        aL = dvL./dt(2:end); aR = dvR./dt(2:end);
        % Savitzky-Golay Filter 
        aL = sgolayfilt(aL,1,3); aR = sgolayfilt(aR,1,3); 
        % all trials
        ALLaL =[ALLaL; aL]; ALLaR=[ALLaR; aR];
        
        
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
    

    subplot(3,7,7+1-conditionSelected);
    % plot the distance
    averagehL=mean(ALLhL,1);averagehR=mean(ALLhR,1);
    plot(t,averagehL,'.r');hold on; plot(t,averagehR,'.b');hold off;
    xlabel('t');ylabel('average distance (pixels)');
    title(['condition ' conditionNames{UniCondi(conditionSelected)}]);
    legend({'Left','Right'},'Location','northwest');
    ylim([0 1000]);
    
    subplot(3,7,14+1-conditionSelected);
    % plot the velocity
    averagevL=mean(ALLvL,1);averagevR=mean(ALLvR,1);
    plot(t(2:end),averagevL,'-r');hold on; plot(t(2:end),averagevR,'-b');hold off;
    xlabel('t');ylabel('average velocity (pixels/1000ms)');legend({'Left','Right'},'Location','northwest');
    title('velocity');
    ylim([-0.2 2]);
    
    
    subplot(3,7,21+1-conditionSelected);
    % plot the acceleration
    averageaL=mean(ALLaL,1);averageaR=mean(ALLaR,1);
    if conditionSelected > 4
    plot(t(3:end),averageaL,'-r');hold on; 
    plot(t(3:end),averageaR,'-b');hold off;
    elseif conditionSelected < 5
    plot(t(3:end),averageaR,'-b'); hold on; 
    plot(t(3:end),averageaL,'-r');hold off;
    end
    xlabel('t');ylabel('average acceleration (pixels/1000ms^2)');legend({'Left','Right'},'Location','southwest');
    title('acceleration');
    ylim([-7e-3 7e-3]);
    %---delete the legend-------------


%     hold off;
end

%% Using the resampled behaviral data and apply Savitzky-Golay Filter (trial level) (for correlation)
% To Add to data_trials as another 4 channels: vL vR aL aR (maybe not, unequal length)
data_trials; 
% To produce the following array
ALLhL=[];ALLhR=[];% travel distance
ALLvL=[];ALLvR=[];% velocity
ALLaL=[];ALLaR=[];% acceleration

for i=1:length(CondiData) %loop through every trial
    % plot(dataxJyJtemp(:,1),dataxJyJtemp(:,2),'g');
    dataxLRJyLRJtemp=cell2mat(RSBehavData{i});
    dataxLRJyLRJtemp=dataxLRJyLRJtemp';

    ALLxL=dataxLRJyLRJtemp(:,1); ALLyL=dataxLRJyLRJtemp(:,2);
    ALLxR=dataxLRJyLRJtemp(:,3); ALLyR=dataxLRJyLRJtemp(:,4);

    t=1:length(ALLxL);

    % distance
    hL=xCenter-ALLxL';hR=ALLxR'-xCenter;
    % all trials
    ALLhL=[ALLhL; hL]; ALLhR=[ALLhR; hR];

    % change in time
    dt=diff(t);
    % change in distance
    dhL=diff(hL); dhR=diff(hR);

    % v - velocity
    vL = dhL./dt; vR = dhR./dt;
    % Savitzky-Golay Filter 
    vL = sgolayfilt(vL,1,251); vR = sgolayfilt(vR,1,251); 
    % all trials
    ALLvL=[ALLvL; vL]; ALLvR=[ALLvR; vR];

    % change in velocity
    dvL = diff(vL); dvR = diff(vR); 
    % a - acceleration
    aL = dvL./dt(2:end); aR = dvR./dt(2:end);
    % Savitzky-Golay Filter 
    aL = sgolayfilt(aL,1,3); aR = sgolayfilt(aR,1,3); 
    % all trials
    ALLaL =[ALLaL; aL]; ALLaR=[ALLaR; aR];
end
    
