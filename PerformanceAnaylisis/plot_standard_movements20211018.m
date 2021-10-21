conditions = [0:1:12]*(pi/2/6);% 13 conditions from [0-180] degrees

conditions


% pick a condition from randomized set allPerm
% conditionSelected = allPerm(numTrials*(block-1)+t);

%% plot out all condition targets
figure;
for conditionSelected=1:length(conditions)
    % produced the position parameters
        rad_ang=conditions(conditionSelected);
        [x,y] = drawReach(radius,rad_ang, xCenter, yCenter);
    plot(x,y,'k');
    hold on;
end
hold off;
ylim([0 screenYpixels]);xlim([0 screenXpixels]);set(gca, 'YDir', 'reverse');


%% plot out all correct movements of standard
figure; es=2.5;

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

    subplot(7,13,abs(7-conditionSelected)*13+(14-conditionSelected));
    plot(x-xCenter,yCenter-y,'k');hold on;
    if conditionSelected==1:6;
        plot([-Ltotal 0],[0 0],'r');
    else
        plot([-Ltotal 0],[0 0],'r');
    end
    if conditionSelected==1:6;
        plot([-Rtotal 0], [0 0],'b');
    else
        plot([Rtotal 0], [0 0],'b');
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

