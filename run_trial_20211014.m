% This edition change the 3/4 circle into 1/2 circle. And remove the connecting dots

while n <  numFrames %~KbCheck 
    % If esc is press, break out of the while loop and close the screen
    [keyIsDown, keysecs, keyCode] = KbCheck;
    if keyCode(KbName('escape'))
        Screen('CloseAll');
        break;
    end

    % Draw the route
    Screen('DrawDots', windowPtr, [x;y], Thickness, white, [0 0], 2);
    
    % Show the movement ratio
    textRatio=conditionNames{conditionSelected};
    DrawFormattedText2(textRatio,'win',windowPtr,...
            'sx',x(end),'sy', y(end)-Thickness, 'xalign','center','yalign','top','baseColor',white);
    
    % Draw guiding arrows
    % create a line and a triangle on the left
    Screen('DrawLine', windowPtr, white, xCenter-screenXpixels/6, yCenter, xCenter-screenXpixels/4, yCenter, 1);
    head   = [ xCenter-screenXpixels/4, yCenter ]; % coordinates of head
    width  = 10;           % width of arrow head
    points = [ head-[0, width]         % left corner
               head+[0, width]         % right corner
               head-[width, 0] ];      % vertex
    Screen('FillPoly', windowPtr, white , points);
    % create a line and a triangle on the right
    Screen('DrawLine', windowPtr, white, xCenter+screenXpixels/6, yCenter, xCenter+screenXpixels/4, yCenter, 1);
    head   = [ xCenter+screenXpixels/4, yCenter ]; % coordinates of head
    width  = 10;           % width of arrow head
    points = [ head-[0, width]         % left corner
               head+[0, width]         % right corner
               head+[width, 0] ];      % vertex
    Screen('FillPoly', windowPtr, white , points);
    
    
    % *******************************************************************
    % Get mouse location
    [xML0, yML0] = GetMouse(windowPtr,mice(2));
    [xMR0, yMR0] = GetMouse(windowPtr,mice(1));
    
    % Shift the mouse location to map the tablets
%     xML=xML0/2; % the upper left quatrand of the scrren
%     yML=yML0/2; % the upper left quatrand of the scrren
%     xMR=xCenter+xMR0/2; % the upper right quatrand of the screen
%     yMR=yMR0/2;% the upper right quatrand of the screen
    % xML=xML0/2; % the left side of the scrren
    xML=xML0*(11/20); % a bit more than half of the left side of the scrren
    yML=yML0; % the left side of the scrren
    % xMR=xCenter+xMR0/2; % the right side of the screen
    xMR=(xCenter-screenXpixels/20)+xMR0*(11/20); % a bit more than half of the right side of the scrren
    yMR=yMR0;% the right side of the screen
    
    
    %********************************************************************
    % limit the mice movement only going outward, cannot return
%     if xML > xMLbefore
%         xML=xMLbefore; 
%     end
%     if xMR < xMRbefore
%         xMR=xMRbefore; 
%     end
    
    %********************************************************************
    % Hide the cursor
    HideCursor(windowPtr,mice(2));
    HideCursor(windowPtr,mice(1));
     
    % Display the cursor as a dot
%     Screen('DrawDots', windowPtr, [xML yML], Thickness, red, [], 2);
%     Screen('DrawDots', windowPtr, [xMR yMR], Thickness, blue, [], 2);


    % *******************************************************************
    % fash the upper left corner once at the start of planing phase (and bottom right)
    if n == 1
    Screen('FillRect', windowPtr, white, LeftUpperSquare);
    Screen('FillRect', windowPtr, white, RightBottomSquare);
    end
    
    % Planning phase (hold the mice)
    if n< numFramesPlan 
%         SetMouse(x(1),y(1), windowPtr, mice(2));
%         SetMouse(x(1),y(1), windowPtr, mice(1)); 
        % Hide the cursor
%         HideCursor(windowPtr,mice(2));
%         HideCursor(windowPtr,mice(1)); 
        textPlan=['Wait...'];
        DrawFormattedText2(textPlan,'win',windowPtr,...
            'sx','center','sy', yCenter+screenYpixels/20,'xalign','center','yalign','top','baseColor',white);
        %Screen('DrawText', windowPtr, ['Place mice at the starting points'], xCenter-steplength/1.5, yCenter, white);
        %DrawFormattedText2(['Place mice at starting point and start drawing in 1 second!'],'win',windowPtr,...
        %'sx','center','sy','center','xalign','center','yalign','center','baseColor',white);
        
    end
    
    %flash the upper right corner once at the end of planning phase (and bottom right)
    if n == numFramesPlan 
    Screen('FillRect', windowPtr, white, RightUpperSquare);
    Screen('FillRect', windowPtr, white, RightBottomSquare);
    end
    
    
    % ******
    % Flash the bottom left corner every other second during the whole trial
    if (~isempty(find([1:2:numFrames]==n))) % every two frames
    Screen('FillRect', windowPtr, white, LeftBottomSquare);
    end
    
    % ******
    % Flash the upper left corner once at the end of the trial (and bottom right)
    if n == (numFrames-1)
    Screen('FillRect', windowPtr, white, LeftUpperSquare);
    Screen('FillRect', windowPtr, white, RightBottomSquare);
    end
    
    
     
    %********************************************************************
    % Enlarge the scale of movement by es times
%     es=screenXpixels/screenYpixels; % it has a limit of 0 ~ (screenXpixels/screenYpixels)
    es=2.25; % 2.5 % for biggest movement
    
    %********************************************************************
%     % Display the joint position as a dot
%     if conditionSelected < median(1:length(conditions));
%         xJ=xCenter+(xMR-xCenter)/es; yJ=yCenter-(xCenter-xML)/es;
% %         if abs(xJ) < abs(xJbefore) | abs(yJ) < abs(yJbefore)
% %             xJ=xJbefore; yJ=yJbefore;
% %         end % not moving at all
%         Screen('DrawDots', windowPtr, [xJ yJ], Thickness, Jcolor, [], 2);
%     end
%     if conditionSelected == median(1:length(conditions));
%         xJ=xCenter+((xMR-xCenter)-(xCenter-xML))/es; yJ=yCenter-(sqrt((xMR-xCenter).^2+(xCenter-xML).^2))/es;
% %         if abs(xJ) < abs(xJbefore) | abs(yJ) < abs(yJbefore)
% %             xJ=xJbefore; yJ=yJbefore;
% %         end % not moving at all
%         Screen('DrawDots', windowPtr, [xJ yJ], Thickness, Jcolor, [], 2);
%     end
%     if conditionSelected > median(1:length(conditions));
%         xJ=xCenter-(xCenter-xML)/es; yJ=yCenter-(xMR-xCenter)/es;
% %         if abs(xJ) < abs(xJbefore) | abs(yJ) < abs(yJbefore)
% %             xJ=xJbefore; yJ=yJbefore;
% %         end % not moving at all
%         Screen('DrawDots', windowPtr, [xJ yJ], Thickness, Jcolor, [], 2);
%     end

    %********************************************************************
    % alternative calculation (imperfection: at 45 and 135 degree, one hand don't have to move at all)
%     xJ=xCenter+((xMR-xCenter)-(xCenter-xML))/es; yJ=yCenter-(sqrt((xMR-xCenter).^2+(xCenter-xML).^2))/es;
%     Screen('DrawDots', windowPtr, [xJ yJ], Thickness, Jcolor, [], 2);

    %********************************************************************
    % another alternative calculation (tangent calculation)
    XL=xCenter-xML; XR=xMR-xCenter;
    distance=sqrt(XL.^2+XR.^2); %travel distance
    theta=2*atan(-XL/XR)
    XJ=distance*cos(theta);
    YJ=distance*sin(theta);
    xJ=xCenter+XJ/es; yJ=yCenter-(-YJ)/es;
    
    %********************************************************************
    % When both hands falls within each other near the curve, both cursors turn green
    % When Joint dot near the line 
    pt=[xJ yJ]; v1=dotCenter; v2=[x(end) y(end)];
    distance=point_to_line_distance(pt, v1, v2);
    withinRadius=norm(pt-dotCenter);
    withinRadiusReverse=norm(pt-v2);
    inside = (distance<Thickness) & (withinRadius < radius+Thickness) & (withinRadiusReverse < radius+Thickness);
%     inside = (distance<Thickness) & (withinRadius < radius+Thickness) & (yJ+Thickness)>0 & ...
%     (xML < xMLbefore) & (xMR > xMRbefore);
    
    if n > numFramesPlan % start counting scores after planning phase
        % MLcolor=red;MRcolor=blue;
        textMove=['Go!'];
        DrawFormattedText2(textMove,'win',windowPtr,...
            'sx','center','sy', yCenter+screenYpixels/20,'xalign','center','yalign','top','baseColor',green);
        % display the joint movement
        Screen('DrawDots', windowPtr, [xJ yJ], Thickness, Jcolor, [], 2);
        if inside == 1
            Jcolor = green;
            
            % Count green frames (not optimal)
%             if ((xML~=xMLbefore)|(yML~=yMLbefore)|(xMR~=xMRbefore)|(yMR~=yMRbefore))==1 % at least one coortinate moved
%                 ScoreLR=ScoreLR+round(fullBonusPerTrial/numFrames,4); % add a bonus score
%             end

            % Count covered green pixels (not optimal)
%             GreenPixelsCoveredstr=strcat(num2str(round(xJ)),num2str(round(xJ)));
%             GreenPixelsCoverednum=str2num(GreenPixelsCoveredstr);
%             GreenPixelsCovered=[GreenPixelsCovered; GreenPixelsCoverednum]; % The green pixels covered in this trial 
%             ScoreRealtime=length(unique(GreenPixelsCovered)); % real time score, get updated in every frame
%             ScoreLR=ScoreRealtime/(radius*Thickness);

            % Count added union of polyshape greenpixel areas
            %PolyshapGreenPixel=polyshape([xJ-Thickness/2 xJ+Thickness/2 xJ+Thickness/2 xJ-Thickness/2],...
            %    [yJ-Thickness/2 yJ-Thickness/2 yJ+Thickness/2 yJ+Thickness/2]);
            
            % the circle version
            PolyshapGreenPixel=polycircle(xJ,yJ,Thickness/2+2); % add 2 pixels of freedom so that can cover 100%
            
            % hold on; plot(PolyshapGreenPixel);ylim([0 screenYpixels]);xlim([0 screenXpixels]);set(gca, 'YDir', 'reverse');
            PolyshapeUnion=union(PolyshapeUnion,intersect(PolyshapGreenPixel,polyTrajatory));
            % figure;plot(PolyshapeUnion);ylim([0 screenYpixels]);xlim([0 screenXpixels]);set(gca, 'YDir', 'reverse');
            
        else % insideLR == 0
            Jcolor = white;
        end
    end
     
    % Flip to the screen
    vbl  = Screen('Flip', windowPtr, vbl + (waitframes -0.5) * ifi);
    
    % Keep track of the nummers of flames when dot is inside the rect
    NumInside(n)=inside;    
        
    % keep a record of the mouse position to be use as previous position in
    % the next loop
    xMLbefore=xML;
    yMLbefore=yML;
    xMRbefore=xMR;
    yMRbefore=yMR;
    xJbefore=xJ;
    yJbefore=yJ;
    
    % Keep track of the mouse trace
    xLyL(n,1) = xML;
    xLyL(n,2) = yML;
    xRyR(n,1) = xMR;
    xRyR(n,2) = yMR;
    xJyJ(n,1)=xJ;
    xJyJ(n,2)=yJ;
    % update n
    n = n+1;


end

    % Reset
    Jcolor = green;