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
    plot([-Ltotal 0],[0 0],'r');
    plot([Rtotal 0], [0 0],'b');
    plot(0,0,'go');hold off;
    ylim([0 screenYpixels]);xlim([0 screenXpixels]);set(gca, 'YDir', 'reverse');
end


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

