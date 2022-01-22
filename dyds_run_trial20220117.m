% This edition change the 3/4 circle into 1/2 circle. And remove the connecting dots

% while n <  numFrames %~KbCheck 
while numberOfSecondsRemaining > 0

    % If esc is press, break out of the while loop and close the screen
    [keyIsDown, keysecs, keyCode] = KbCheck;
    if keyCode(KbName('escape'))
        Screen('CloseAll');
        break;
    end
    
    % Update the while loop with time
    numberOfSecondsElapsed = round((now - startTime) * 10 ^ 5);
    numberOfSecondsRemaining = durationInSeconds - numberOfSecondsElapsed;
    % Display time remaining
    Screen('DrawText', windowPtr, 'Press the left or right arrow key to move, or the escape key to quit.', 20, 20, white);
    Screen('DrawText', windowPtr, sprintf('%i seconds remaining...', numberOfSecondsRemaining), 20, 50, white);


    xOffsetL = 0; xOffsetR = 0;
    yOffsetL = elevationL; yOffsetR = elevationR;
    offsetspotRectL = OffsetRect(startspotRectL, xOffsetL, -yOffsetL);
    offsetspotRectR = OffsetRect(startspotRectR, xOffsetR, -yOffsetR);
    offsetbarRect   = startbarRect + [xOffsetL -yOffsetL xOffsetR -yOffsetR];
    Screen('FillOval', windowPtr, red, offsetspotRectL');
    Screen('FillOval', windowPtr, blue, offsetspotRectR');
    if abs(elevationR-elevationL)<3*spotDiameter
        Screen('DrawLine', windowPtr, green, screenXpixels*1/3-spotRadius , screenYpixels-3*spotRadius-yOffsetL, screenXpixels*2/3+spotRadius, screenYpixels-3*spotRadius-yOffsetR ,10);
        % Flip to the screen
        vbl  = Screen('Flip', windowPtr, vbl + (waitframes -0.5) * ifi);
    else
        Screen('DrawLine', windowPtr, red, screenXpixels*1/3+4*spotDiameter , screenYpixels-3*spotRadius-yOffsetL, screenXpixels*2/3-4*spotDiameter, screenYpixels-3*spotRadius-yOffsetR ,10);
        % Flip to the screen
        vbl  = Screen('Flip', windowPtr, vbl + (waitframes -0.5) * ifi);
        break;
    end
    

     


    Screen('FillOval', windowPtr0, red, offsetspotRectL');
    Screen('DrawDots', windowPtr0, [FixCrX;FixCrY-yOffsetL], Thickness/10, white, [0 0], 2);
    Screen('FillOval', windowPtr0, blue, offsetspotRectR');
    Screen('DrawLine', windowPtr0, green, screenXpixels*1/3-spotRadius , screenYpixels-3*spotRadius-yOffsetL, screenXpixels*2/3+spotRadius, screenYpixels-3*spotRadius-yOffsetR ,10);

    vbl  = Screen('Flip', windowPtr0, vbl + (waitframes -0.5) * ifi);

    % Update the while loop
    % [ keyIsDown, seconds, keyCode ] = KbCheck;
    [keyboardIndices, productNames, allInfos] = GetKeyboardIndices
%     deviceNumberL=10;deviceNumberR=15;
    deviceNumberL=6;deviceNumberR=10;
    [keyIsDownL, secsL, keyCodeL, deltaSecsL] = KbCheck([deviceNumberL]);
    [keyIsDownR, secsR, keyCodeR, deltaSecsR] = KbCheck([deviceNumberR]);

    escapeKey = KbName('ESCAPE');

    if keyIsDownL
        if keyCodeL (leftKey)
            if elevationL < screenYpixels - spotDiameter
                elevationL = elevationL + spotDiameter;
            else
                elevationL = 0;
            end
        end
    end

    if keyIsDownR
        if keyCodeR (rightKey)
            if elevationR < screenYpixels - spotDiameter
                elevationR = elevationR + spotDiameter;
            else
                elevationR = 0;
            end
        end
    end
    if keyCodeL (escapeKey) | keyCodeR (escapeKey) 
            break;
    end

    % If the user holds down a key, KbCheck will report multiple events.
    % To condense multiple 'keyDown' events into a single event, we wait until all
    % keys have been released.
    KbReleaseWait;

end

    % Reset
%     Jcolor = green;