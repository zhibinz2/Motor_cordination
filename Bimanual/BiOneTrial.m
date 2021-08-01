while n <  numFrames %~KbCheck

    % If esc is press, break out of the loop and close the screen
    [keyIsDown, keysecs, keyCode] = KbCheck;
    if keyCode(KbName('escape'))
        Screen('CloseAll');
        break;
    end
    

    % Draw the route
    Screen('DrawDots', windowPtr, PosL, ThicknessL, ColorL, [0 0], 2);
    Screen('DrawDots', windowPtr, PosR, ThicknessR, ColorR, [0 0], 2);
    % Draw the connecting points
    Screen('DrawDots', windowPtr, PosL0, ConnectDotSize, ConnectDotColorL0, [0 0], 2);
    Screen('DrawDots', windowPtr, PosL1, ConnectDotSize, ConnectDotColorL1, [0 0], 2);
    Screen('DrawDots', windowPtr, PosL2, ConnectDotSize, ConnectDotColorL2, [0 0], 2);
    Screen('DrawDots', windowPtr, PosL3, ConnectDotSize, ConnectDotColorL3, [0 0], 2);
    Screen('DrawDots', windowPtr, PosR0, ConnectDotSize, ConnectDotColorR0, [0 0], 2);
    Screen('DrawDots', windowPtr, PosR1, ConnectDotSize, ConnectDotColorR1, [0 0], 2);
    Screen('DrawDots', windowPtr, PosR2, ConnectDotSize, ConnectDotColorR2, [0 0], 2);
    Screen('DrawDots', windowPtr, PosR3, ConnectDotSize, ConnectDotColorR3, [0 0], 2);
 
    % position of the right-bottom to draw the dot
%     RightBottom = [screenXpixels-PhotosensorSize screenYpixels-PhotosensorSize]; 
%     if n == 1
%         % Draw a white dot on the right-bottom corner of the screen for the photosensor
%         Screen('DrawDots', windowPtr, RightBottom, PhotosensorSize, ones(3,1), [0 0], 2);
%     elseif n == length(steplength)-1
%         Screen('DrawDots', windowPtr, RightBottom, PhotosensorSize, ones(3,1), [0 0], 2);
%     elseif n == 2*length(steplength)-1
%         Screen('DrawDots', windowPtr, RightBottom, PhotosensorSize, ones(3,1), [0 0], 2);
%     elseif n == length(PosL)-1
%         Screen('DrawDots', windowPtr, RightBottom, PhotosensorSize, ones(3,1), [0 0], 2);
%     end 
%   
    % Flash the four corners every other 5 s
    if (~isempty(find([1:5:numFrames]==n))) % every five frames
    Screen('FillRect', windowPtr, white, RightBottomSquare);
    Screen('FillRect', windowPtr, white, RightUpperSquare);
    Screen('FillRect', windowPtr, white, LeftBottomSquare);
    Screen('FillRect', windowPtr, white, LeftUpperSquare);
    end
    
    % Planning phase
    if n<60
        Screen('DrawText', windowPtr, ['You have 3 sec to draw, start now!'], screenXpixels/3, screenYpixels/2, white); 
    end
    
%     % Locate the moving rectangle
%     centeredRect = CenterRectOnPointd(baseRect, PosL(1,n), PosL(2,n));
%     % Draw the rect to the screen
%     Screen('FillRect', windowPtr, rectColor, centeredRect);

    % Setting default mouse Position
    SetMouse(PosL(1,1),PosL(2,1), window, mice(2));
    SetMouse(PosR(1,1),PosR(2,1), window, mice(1));
    
    % Get mouse location
    [xML, yML] = GetMouse(windowPtr,mice(2));
    [xMR, yMR] = GetMouse(windowPtr,mice(1));
    
    % Display the cursor as a dot
    Screen('DrawDots', windowPtr, [xML yML], 16, red, [], 2);
    Screen('DrawDots', windowPtr, [xMR yMR], 16, blue, [], 2);
    
    % See if the mouse cursor is inside the connecting dots
    % Left side
    insideL0 = IsInDot(xML, yML, PosL0, ConnectDotSize);
    if insideL0 == 1
        ConnectDotColorL0 = red;
    elseif insideL0 == 0
        ConnectDotColorL0 = white;
    end
    insideL1 = IsInDot(xML, yML, PosL1, ConnectDotSize);
    if insideL1 == 1
        ConnectDotColorL1 = red;
    elseif insideL1 == 0
        ConnectDotColorL1 = white;
    end
    inside2L = IsInDot(xML, yML, PosL2, ConnectDotSize);
    if inside2L == 1
        ConnectDotColorL2 = red;
    elseif inside2L == 0
        ConnectDotColorL2 = white;
    end
    inside3L = IsInDot(xML, yML, PosL3, ConnectDotSize);
    if inside3L == 1
        ConnectDotColorL3 = red;
    elseif inside3L == 0
        ConnectDotColorL3 = white;
    end
    % Right side
    insideR0 = IsInDot(xMR, yMR, PosR0, ConnectDotSize);
    if insideR0 == 1
        ConnectDotColorR0 = blue;
    elseif insideR0 == 0
        ConnectDotColorR0 = white;
    end
    insideR1 = IsInDot(xMR, yMR, PosR1, ConnectDotSize);
    if insideR1 == 1
        ConnectDotColorR1 = blue;
    elseif insideR1 == 0
        ConnectDotColorR1 = white;
    end
    inside2R = IsInDot(xMR, yMR, PosR2, ConnectDotSize);
    if inside2R == 1
        ConnectDotColorR2 = blue;
    elseif inside2R == 0
        ConnectDotColorR2 = white;
    end
    inside3R = IsInDot(xMR, yMR, PosR3, ConnectDotSize);
    if inside3R == 1
        ConnectDotColorR3 = blue;
    elseif inside3R == 0
        ConnectDotColor3R = white;
    end
    % Both side
    insideBi0 = IsInDot(xML, yML, PosL0, ConnectDotSize) & IsInDot(xMR, yMR, PosR0, ConnectDotSize);
    if insideBi0 == 1
        ConnectDotColorL0 = green;
        ConnectDotColorR0 = green;
    elseif insideBi0 == 0
        ConnectDotColorL0 = white;
        ConnectDotColorR0 = white;
    end
    insideBi1 = IsInDot(xML, yML, PosL1, ConnectDotSize) & IsInDot(xMR, yMR, PosR1, ConnectDotSize);
    if insideBi1 == 1
        ConnectDotColorL1 = green;
        ConnectDotColorR1 = green;
    elseif insideBi1 == 0
        ConnectDotColorL1 = white;
        ConnectDotColorR1 = white;
    end
    inside2R = IsInDot(xML, yML, PosL2, ConnectDotSize) & IsInDot(xMR, yMR, PosR2, ConnectDotSize);
    if inside2R == 1
        ConnectDotColorL2 = green;
        ConnectDotColorR2 = green;
    elseif inside2R == 0
        ConnectDotColorL2 = white;
        ConnectDotColorR2 = white;
    end
    inside3R = IsInDot(xML, yML, PosL3, ConnectDotSize) & IsInDot(xMR, yMR, PosR3, ConnectDotSize);
    if inside3R == 1
        ConnectDotColor3L = green;
        ConnectDotColor3R = green;
    elseif inside3R == 0
        ConnectDotColor3L = white;
        ConnectDotColor3R = white;
    end
    
    
    % Flip to the screen
    vbl  = Screen('Flip', windowPtr, vbl + (waitframes -0.5) * ifi);
    
    % Keep track of the nummers of flames when dot is inside the rect
%     NumInside(n)=inside;    
    
    % Keep track of the mouse trace
    xLyL(n,1) = xML;
    xLyL(n,2) = yML;
    xRyR(n,1) = xMR;
    xRyR(n,2) = yMR;
    % update n
    n = n+1;

end

