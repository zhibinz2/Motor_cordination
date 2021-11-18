% conditions = [0:1:12]*(pi/2/6);% 13 conditions from [0-180] degrees


conditions = [0 2*atan(1/4) 2*atan(1/2) 2*atan(1) 2*atan(2/1) 2*atan(4/1) pi];
conditionNames={'0:4' '1:4' '1:2' '1:1' '2:1' '4:1' '4:0'};  % conditionNames{1}

% pick a condition from randomized set allPerm
% conditionSelected = allPerm(numTrials*(block-1)+t);

%% plot out all condition targets
figure;
for conditionSelected=1:length(conditions)
    % produced the position parameters
        rad_ang=conditions(conditionSelected);
        [x,y] = drawReach(radius,rad_ang, xCenter, yCenter);
    plot(x,y,'linewidth',3,'color',[rand(1,3)]);
    hold on;
end
hold off;
ylim([0 screenYpixels]);xlim([0 screenXpixels]);set(gca, 'YDir', 'reverse');
xlabel('Screen-x-pixel');ylabel('Screen-y-pixel');
legend(conditionNames);


%% plot out all correct movements of standard
figure; es=2.25;

%     % another alternative calculation (tangent calculation)
%     XL=xCenter-xML; XR=xMR-xCenter;
%     distance=sqrt(XL.^2+XR.^2); %travel distance
%     theta=2*atan(-XL/XR)
%     XJ=distance*cos(theta);
%     YJ=distance*sin(theta);
%     xJ=xCenter+XJ/es; yJ=yCenter-(-YJ)/es;
    
    
for conditionSelected=1:length(conditions)
    % produced the position parameters
        rad_ang=conditions(conditionSelected);
        [x,y] = drawReach(radius,rad_ang, xCenter, yCenter);   
    % theta
    theta=rad_ang/2;
    
    % left hand movement
%     Ltotal=radius*sin(theta)*sqrt(sqrt(2))*es;
    Ltotal=radius*sin(theta)*es;
%     Ltotal=(yCenter-y(end)/es);

    % right hand movement
%     Rtotal=radius*cos(theta)*sqrt(sqrt(2))*es;
    Rtotal=radius*cos(theta)*es;
%     Rtotal=(x(end)/es-xCenter);

%     subplot(7,13,abs(7-conditionSelected)*13+(14-conditionSelected));
%     subplot(4,7,abs(4-conditionSelected)*7+(7+1-conditionSelected));
    subplot(1,7,7+1-conditionSelected);
    plot(x-xCenter,yCenter-y,'g','linewidth',2);hold on;
    xlabel('Screen-x-pixel');ylabel('Screen-y-pixel');
    title([conditionNames{conditionSelected}]);
    if conditionSelected==1:6;
        plot([-Ltotal 0],[0 0],'r','linewidth',2);
    else
        plot([-Ltotal 0],[0 0],'r','linewidth',2);
    end
    if conditionSelected==1:6;
        plot([-Rtotal 0], [0 0],'b','linewidth',2);
    else
        plot([Rtotal 0], [0 0],'b','linewidth',2);
    end
    plot(0,0,'go');hold off;
    ylim([-screenYpixels/2 screenYpixels/2]);xlim([-screenXpixels/2 screenXpixels/2]); 
    
    
%     plot(x-xCenter,yCenter-y,'k');hold on;
%     plot([-Ltotal 0],[0 0],'r');
%     plot([Rtotal 0], [0 0],'b');
%     plot(0,0,'go');hold off;
%     ylim([0 screenYpixels]);xlim([0 screenXpixels]);set(gca, 'YDir', 'reverse');
end

%% wrong calculation
figure; es=1;
for conditionSelected=1:length(conditions)
    % produced the position parameters
        rad_ang=conditions(conditionSelected);
        [x,y] = drawReach(radius,rad_ang, xCenter, yCenter);
        
    % left hand movement
    Ltotal=(yCenter-y(end)/es);
    % right hand movement
    Rtotal=(x(end)/es-xCenter);
    
    subplot(1,length(conditions),conditionSelected);
    plot(x-xCenter,yCenter-y,'k');hold on;
    if conditionSelected==1:6;
        plot([-Ltotal 0],[0 0],'r');
    else
        plot([-Ltotal 0],[0 0],'r');
    end
    if conditionSelected==1:6;
        plot([Rtotal 0], [0 0],'b');
    else
        plot([-Rtotal 0], [0 0],'b');
    end
    plot(0,0,'go');hold off;
    ylim([-screenYpixels/2 screenYpixels/2]);xlim([-screenXpixels/2 screenXpixels/2]); 
end

