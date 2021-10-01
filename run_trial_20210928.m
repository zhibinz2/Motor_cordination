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

    % ******
    % fash the upper left corner once at the start of planing phase (and bottom right)
    if n == 1
    Screen('FillRect', windowPtr, white, LeftUpperSquare);
    Screen('FillRect', windowPtr, white, RightBottomSquare);
    end
    
    % Planning phase (hold the mice)
    if n< numFramesPlan 
        SetMouse(x(1),y(1), windowPtr, mice(2));
        SetMouse(x(1),y(1), windowPtr, mice(1)); 
        textPlan=['Place mice at the starting points'];
        DrawFormattedText2(textPlan,'win',windowPtr,...
            'sx','center','sy','center','xalign','center','yalign','center','baseColor',white);
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
    
    
    % Get mouse location
    [xML, yML] = GetMouse(windowPtr,mice(2));
    [xMR, yMR] = GetMouse(windowPtr,mice(1));
    
    % Display the cursor as a dot
    Screen('DrawDots', windowPtr, [xML yML], 30, red, [], 2);
    Screen('DrawDots', windowPtr, [xMR yMR], 30, blue, [], 2);

    % Display the joint position as a dot
    if (conditionSelected < length(conditions)/2) | (conditionSelected == 2);
        xJ=xMR; yJ=yCenter-(xCenter-xML);
        Screen('DrawDots', windowPtr, [xJ yJ], 30, Jcolor, [], 2);
    end
    if conditionSelected > length(conditions)/2;
        xJ=xML; yJ=yCenter-(xMR-xCenter);
        Screen('DrawDots', windowPtr, [xJ yJ], 30, Jcolor, [], 2);
    end
    
    % When both hands falls within each other near the curve, both cursors turn green
    % When Joint dot near the line
    pt=[xJ yJ]; v1=dotCenter; v2=[x(end) y(end)];
    distance=point_to_line_distance(pt, v1, v2);
    withinRadius=norm(pt-dotCenter);
    inside = (distance<ConnectDotSize) & (withinRadius < radius) & (yJ+ConnectDotSize)>0;
    
    if n > numFramesPlan % start counting scores after planning phase
        % MLcolor=red;MRcolor=blue;
        textMove=['Move now!'];
        DrawFormattedText2(textMove,'win',windowPtr,...
            'sx','center','sy','center','xalign','center','yalign','center','baseColor',white);
        %Screen('DrawText', windowPtr, ['Move now!'], xCenter, yCenter, white);
        if inside == 1
            Jcolor = green;
            if ((xML~=xMLbefore)|(yML~=yMLbefore)|(xMR~=xMRbefore)|(yMR~=yMRbefore))==1 % at least one coortinate moved
                ScoreLR=ScoreLR+round(fullBonusPerTrial/numFrames,4); % add a bonus score
            end
        else % insideLR == 0
            Jcolor = white;
        end
    end
     
    % Flip to the screen
    vbl  = Screen('Flip', windowPtr, vbl + (waitframes -0.5) * ifi);
    
    % Keep track of the nummers of flames when dot is inside the rect
    % NumInside(n)=inside;    
        
    % keep a record of the mouse position to be use as previous position in
    % the next loop
    xMLbefore=xML;
    yMLbefore=yML;
    xMRbefore=xMR;
    yMRbefore=yMR;
    
    % Keep track of the mouse trace
    xLyL(n,1) = xML;
    xLyL(n,2) = yML;
    xRyR(n,1) = xMR;
    xRyR(n,2) = yMR;
    % update n
    n = n+1;


end

    % Reset
    Jcolor = green;