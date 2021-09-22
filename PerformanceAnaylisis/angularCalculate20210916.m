%% obtain x y
[x,y] = drawA1(steplength);

%% find the center of the circle
CenterXL = steplength + steplength/2; 
CenterYL = steplength + (yCenter - max(y) + (max(y) - min(y))/2); 
CenterXR = steplength + screenXpixels/2 + steplength/6; 
CenterYR = steplength + (yCenter - max(y) + (max(y) - min(y))/2); 

% plot the routes to check
plot(PosL(1,:),PosL(2,:),'-g', 'LineWidth',10);
hold on;
plot(PosR(1,:),PosR(2,:),'-g','LineWidth',10); 
hold on;
set(gca, 'YDir', 'reverse');
% plot the two ceter point
plot(CenterXL, CenterYL,'go');
plot(CenterXR, CenterYR,'go');

%% Plot one last trace to check
plot(TracexL,TraceyL,'-r',TracexR,TraceyR,'-b');

%% find the route start from 0.75*pi to 2.25*pi (condition A2&A2)
StartPointxL=PosL(1,1);StartPointyL=PosL(2,1); 
plot(StartPointxL, StartPointyL, 'ro');
StartPointxR=PosR(1,1);StartPointyR=PosR(2,1);
plot(StartPointxR, StartPointyR, 'Bo');


%% Vector of the start side 0.75*pi (condition A2)
StartSideL=[(StartPointxL-CenterXL) (StartPointyL-CenterYL)];
StartSideR=[(StartPointxR-CenterXR) (StartPointyR-CenterYR)];

figure;
hold on;
plot(StartSideL(1),StartSideL(2),'ro');
plot(StartSideR(1),StartSideR(2),'bo');
xlim([-2*steplength 2*steplength]);ylim([-2*steplength 2*steplength]);
set(gca, 'YDir', 'reverse');

for ang=[1 500 1000 1500 2000 3000 4000 5000 6000 7000 8000 9000 10000]
    MovingSideL=[(TracexL(ang)-CenterXL) (TraceyL(ang)-CenterYL)];
    MovingSideR=[(TracexR(ang)-CenterXR) (TraceyR(ang)-CenterYR)];
    hold on;
    plot(MovingSideL(1),MovingSideL(2),'ro');
    plot(MovingSideR(1),MovingSideR(2),'bo');
    pause(0.5);
end

%% https://www.mathworks.com/matlabcentral/answers/364763-how-to-calculate-angle-from-coordinate
% https://stackoverflow.com/questions/1211212/how-to-calculate-an-angle-from-three-points
% https://www.mathworks.com/matlabcentral/answers/101590-how-can-i-determine-the-angle-between-two-vectors-in-matlab


figure; hold on;
plot(StartSideL(1),StartSideL(2),'ro');
plot(StartSideR(1),StartSideR(2),'bo');
xlim([-2*steplength 2*steplength]);ylim([-2*steplength 2*steplength]);
set(gca, 'YDir', 'reverse');

% % ****************Play in For loop take too much time
% for ang=1:length(TracexL)
%     MovingSideL=[(TracexL(ang)-CenterXL) (TraceyL(ang)-CenterYL)];
%     MovingSideR=[(TracexR(ang)-CenterXR) (TraceyR(ang)-CenterYR)];
%     plot(MovingSideL(1),MovingSideL(2),'ro');
%     plot(MovingSideR(1),MovingSideR(2),'bo');
% %     pause(0.0001);
%     AngL(ang) = atan2(abs(det([MovingSideL;StartSideL])),dot(MovingSideL,StartSideL)); 
%     AngR(ang) = atan2(abs(det([MovingSideR;StartSideR])),dot(MovingSideR,StartSideR)); 
%     DiffLR(ang) = AngR(ang)-AngL(ang);
% end
%     
% hold off;***************

% calculate all at once
MovingSideL=[(TracexL-CenterXL) (TraceyL-CenterYL)];
MovingSideR=[(TracexR-CenterXR) (TraceyR-CenterYR)];

plot(MovingSideL(:,1),MovingSideL(:,2),'ro');
plot(MovingSideR(:,1),MovingSideR(:,2),'bo');
    
    for ang=1:length(TracexL)
        % temperal
        AngL(ang)= atan2(abs(det([MovingSideL(ang,:);StartSideL])),dot(MovingSideL(ang,:),StartSideL)); 
        AngR(ang)= atan2(abs(det([MovingSideR(ang,:);StartSideR])),dot(MovingSideR(ang,:),StartSideR)); 
%         AngL(ang)= atan2(norm(cross(MovingSideL(ang,:),StartSideL)),dot(MovingSideL(ang,:),StartSideL)); 
%         AngR(ang)= atan2(norm(cross(MovingSideR(ang,:),StartSideR)),dot(MovingSideR(ang,:),StartSideR)); 
        % spatial
        RadiL(ang)=rssq(MovingSideL(ang,:));
        RadiR(ang)=rssq(MovingSideR(ang,:));
    end

% Temperal
plot(AngL,'r');hold on; plot(AngR,'b');    
DiffLR= AngR-AngL;
plot(DiffLR);
legend({'L','R','R-F'});title('temperal accuracy');



% rssq(DiffLR);
% Root Mean Squared Error
rmseAng=sqrt(mean(DiffLR).^2); % in radiaus

% Spatial
DiffRadiL=RadiL-steplength;
DiffRadiR=RadiR-steplength;
plot(DiffRadiL,'r');hold on;plot(DiffRadiR,'b');
legend({'L','R'});title('spatial accuracy')
% Root Mean Squared Error
rmseSpatialL=sqrt(mean(DiffRadiL).^2); % in radiaus
rmseSpatialR=sqrt(mean(DiffRadiR).^2); % in radiaus

%% trial angle difference

allDiffLR=[]; allAngL=[];allAngR=[];


for u=1:length(UniCondi)
    
    indtemp=find(CondiData==UniCondi(u));
    
    [xL,yL,PosL,PosL0,PosL1,PosL2,PosL3,ThicknessL,ColorL,...
    xR,yR,PosR,PosR0,PosR1,PosR2,PosR3,ThicknessR,ColorR] ...
    = conditionfuctions{UniCondi(u)}(steplength,yCenter,screenXpixels);
    
    % find starting point
    StartPointxL=PosL(1,1);StartPointyL=PosL(2,1); 
    StartPointxR=PosR(1,1);StartPointyR=PosR(2,1);
    StartSideL=[(StartPointxL-CenterXL) (StartPointyL-CenterYL)];
    StartSideR=[(StartPointxR-CenterXR) (StartPointyR-CenterYR)];
    
    for indt=1:length(indtemp)
        TracexL=squeeze(data_trials(:,NumEEGChannels+1,indtemp(indt)));
        TraceyL=squeeze(data_trials(:,NumEEGChannels+2,indtemp(indt)));
        TracexR=squeeze(data_trials(:,NumEEGChannels+3,indtemp(indt)));
        TraceyR=squeeze(data_trials(:,NumEEGChannels+4,indtemp(indt)));
    %   plot(TracexL,TraceyL,'-r',TracexR,TraceyR,'-b');

        % calculate all at once
        MovingSideL=[(TracexL-CenterXL) (TraceyL-CenterYL)];
        MovingSideR=[(TracexR-CenterXR) (TraceyR-CenterYR)];

            for ang=1:length(TracexL)
                % temperal
                AngL(ang)= atan2(abs(det([MovingSideL(ang,:);StartSideL])),dot(MovingSideL(ang,:),StartSideL)); 
                AngR(ang)= atan2(abs(det([MovingSideR(ang,:);StartSideR])),dot(MovingSideR(ang,:),StartSideR)); 

            end
        % Temperal
        allAngL=[allAngL; AngL];
        allAngR=[allAngR; AngR];
        DiffLR= AngR-AngL;
        allDiffLR=[allDiffLR;DiffLR];
        
        

    end
    
    meanAngL=mean(allAngL,1);
    meanAngR=mean(allAngR,1);
    MeanAllDiffLR=mean(allDiffLR,1);
    timestamp=[1:length(AngL)]/2;
    
    subplot(1,4,u);
    h(1:50)=plot(timestamp,allAngL,'Color',[1, 0, 0, 0.2]);hold on; h(51:100)=plot(timestamp,meanAngL,'r','LineWidth',5);
    h(101:150)=plot((timestamp),allAngR,'Color',[0, 0, 1, 0.2]);   h(151:200)=plot((timestamp),meanAngR,'b','LineWidth',5);   
    h(201:250)=plot(timestamp,allDiffLR,'Color',[0, 1, 0, 0.2]);    h(251:300)=plot((timestamp),MeanAllDiffLR,'g','LineWidth',5);   
    xlabel('time');ylabel('angle (radian)');
    functionName=func2str(conditionfuctions{UniCondi(u)});
    functionName=strrep(functionName, '_', '-');
    title(['condition ' functionName 'realtime angles']);
    
    hleg=legend(h(1:50:300),{'L','L mean','R','R mean','R-L', 'R-L mean'});

    hold off;
    
    
    allDiffLR=[]; allAngL=[];allAngR=[];
end







%% Plot temperal and spatial error over 50 trials of 4 blocks

rmseAngTrend=[]; rmseSpatialLTrend=[];rmseSpatialRTrend=[];

for u=1:length(UniCondi)
    
    indtemp=find(CondiData==UniCondi(u));
    
    [xL,yL,PosL,PosL0,PosL1,PosL2,PosL3,ThicknessL,ColorL,...
    xR,yR,PosR,PosR0,PosR1,PosR2,PosR3,ThicknessR,ColorR] ...
    = conditionfuctions{UniCondi(u)}(steplength,yCenter,screenXpixels);
    
    % find starting point
    StartPointxL=PosL(1,1);StartPointyL=PosL(2,1); 
    StartPointxR=PosR(1,1);StartPointyR=PosR(2,1);
    StartSideL=[(StartPointxL-CenterXL) (StartPointyL-CenterYL)];
    StartSideR=[(StartPointxR-CenterXR) (StartPointyR-CenterYR)];
    
    for indt=1:length(indtemp)
        TracexL=squeeze(data_trials(:,NumEEGChannels+1,indtemp(indt)));
        TraceyL=squeeze(data_trials(:,NumEEGChannels+2,indtemp(indt)));
        TracexR=squeeze(data_trials(:,NumEEGChannels+3,indtemp(indt)));
        TraceyR=squeeze(data_trials(:,NumEEGChannels+4,indtemp(indt)));
    %   plot(TracexL,TraceyL,'-r',TracexR,TraceyR,'-b');

        % calculate all at once
        MovingSideL=[(TracexL-CenterXL) (TraceyL-CenterYL)];
        MovingSideR=[(TracexR-CenterXR) (TraceyR-CenterYR)];

            for ang=1:length(TracexL)
                % temperal
                AngL(ang)= atan2(abs(det([MovingSideL(ang,:);StartSideL])),dot(MovingSideL(ang,:),StartSideL)); 
                AngR(ang)= atan2(abs(det([MovingSideR(ang,:);StartSideR])),dot(MovingSideR(ang,:),StartSideR)); 
                % spatial
                RadiL(ang)=rssq(MovingSideL(ang,:));
                RadiR(ang)=rssq(MovingSideR(ang,:));
            end
        % Temperal
        DiffLR= AngR-AngL;
        rmseAng=sqrt(mean(DiffLR).^2); % in radiaus
        
        % Spatial
        DiffRadiL=RadiL-steplength;
        DiffRadiR=RadiR-steplength;
        rmseSpatialL=sqrt(mean(DiffRadiL).^2); % in radiaus
        rmseSpatialR=sqrt(mean(DiffRadiR).^2); % in radiaus
        
        %Cumulate the rmse
        rmseAngTrend=[rmseAngTrend rmseAng];
        rmseSpatialLTrend=[rmseSpatialLTrend rmseSpatialL];
        rmseSpatialRTrend=[rmseSpatialRTrend rmseSpatialR];
    end
    
    subplot(4,2,2*u-1);
    plot(rmseAngTrend);xlabel('trials');ylabel('time acurracy (radian)');
    functionName=func2str(conditionfuctions{UniCondi(u)});
    functionName=strrep(functionName, '_', '-');
    title(['condition ' functionName]);
    subplot(4,2,2*u);
    plot(rmseSpatialLTrend,'r');xlabel('trials');ylabel('spatial acurracy (pixels)');
    hold on;
    plot(rmseSpatialRTrend,'b');legend({'L','R'});
    title(['condition ' functionName]);
    hold off;
    
    rmseAngTrend=[]; rmseSpatialLTrend=[];rmseSpatialRTrend=[];
end


