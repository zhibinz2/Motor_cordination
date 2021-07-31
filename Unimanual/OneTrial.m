while n <  numFrames %~KbCheck

    % If esc is press, break out of the loop and close the screen
    [keyIsDown, keysecs, keyCode] = KbCheck;
    if keyCode(KbName('escape'))
        Screen('CloseAll');
        break;
    end
    

    % Draw the route
    Screen('DrawDots', windowPtr, Pos, Thickness, Color, [0 0], 2);
    % Draw the connecting points
    Screen('DrawDots', windowPtr, Pos0, ConnectDotSize, ConnectDotColor0, [0 0], 2);
    Screen('DrawDots', windowPtr, Pos1, ConnectDotSize, ConnectDotColor1, [0 0], 2);
    Screen('DrawDots', windowPtr, Pos2, ConnectDotSize, ConnectDotColor2, [0 0], 2);
    Screen('DrawDots', windowPtr, Pos3, ConnectDotSize, ConnectDotColor3, [0 0], 2);
    
    % position of the right-bottom to draw the dot
%     RightBottom = [screenXpixels-PhotosensorSize screenYpixels-PhotosensorSize]; 
%     if n == 1
%         % Draw a white dot on the right-bottom corner of the screen for the photosensor
%         Screen('DrawDots', windowPtr, RightBottom, PhotosensorSize, ones(3,1), [0 0], 2);
%     elseif n == length(steplength)-1
%         Screen('DrawDots', windowPtr, RightBottom, PhotosensorSize, ones(3,1), [0 0], 2);
%     elseif n == 2*length(steplength)-1
%         Screen('DrawDots', windowPtr, RightBottom, PhotosensorSize, ones(3,1), [0 0], 2);
%     elseif n == length(Pos)-1
%         Screen('DrawDots', windowPtr, RightBottom, PhotosensorSize, ones(3,1), [0 0], 2);
%     end 
%   
    % Position the right-bottom to draw a square 
    if (~isempty(find([1:5:numFrames]==n))) % every five frames
    Screen('FillRect', windowPtr, white, RightBottomSquare);
    Screen('FillRect', windowPtr, white, RightUpperSquare);
    Screen('FillRect', windowPtr, white, LeftBottomSquare);
    Screen('FillRect', windowPtr, white, LeftUpperSquare);
    end
    
%     % Locate the moving rectangle
%     centeredRect = CenterRectOnPointd(baseRect, Pos(1,n), Pos(2,n));
%     % Draw the rect to the screen
%     Screen('FillRect', windowPtr, rectColor, centeredRect);

    % Get mouse location
    [xM, yM, buttons] = GetMouse(windowPtr);

    % See if the mouse cursor is inside the connecting dots
    inside0 = IsInDot(xM, yM, Pos0, ConnectDotSize);
    if inside0 == 1
        ConnectDotColor0 = red;
    elseif inside0 == 0
        ConnectDotColor0 = blue;
    end
    inside1 = IsInDot(xM, yM, Pos1, ConnectDotSize);
    if inside1 == 1
        ConnectDotColor1 = red;
    elseif inside1 == 0
        ConnectDotColor1 = blue;
    end
    inside2 = IsInDot(xM, yM, Pos2, ConnectDotSize);
    if inside2 == 1
        ConnectDotColor2 = red;
    elseif inside2 == 0
        ConnectDotColor2 = blue;
    end
    inside3 = IsInDot(xM, yM, Pos3, ConnectDotSize);
    if inside3 == 1
        ConnectDotColor3 = red;
    elseif inside3 == 0
        ConnectDotColor3 = blue;
    end
    
    % Display the cursor as a dot
    Screen('DrawDots', windowPtr, [xM yM], 16, grey, [], 2);
    
    % Flip to the screen
    vbl  = Screen('Flip', windowPtr, vbl + (waitframes -0.5) * ifi);
    
    % Keep track of the nummers of flames when dot is inside the rect
%     NumInside(n)=inside;    
    
    % Keep track of the mouse trace
    xy(n,1) = xM;
    xy(n,2) = yM;
    n = n+1;
%     Screen('Flip', windowPtr); 
end

