% This version uses two Cedrus RB

% run SuppressWarning.m 

%initialize variables
pressedL1=[]; RBkeyL1=[]; %portL1=[]; stampL1=[];
pressedL2=[]; RBkeyL2=[]; %portL2=[]; stampL2=[];
pressedR1=[]; RBkeyR1=[]; %portR1=[]; stampR1=[];
pressedR2=[]; RBkeyR2=[]; %portR2=[]; stampR2=[];
elevationL = 0; elevationR = 0;

startTime = now;
numberOfSecondsRemaining=durationInSeconds;

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

    offsetspotRectL = OffsetRect(startspotRectL, 0, -elevationL);
    offsetspotRectR = OffsetRect(startspotRectR, 0, -elevationR);
    offsetbarRect   = startbarRect + [0 -elevationL 0 -elevationR];
    Screen('FillOval', windowPtr, red, offsetspotRectL');
    Screen('FillOval', windowPtr, blue, offsetspotRectR');
    if abs(elevationR-elevationL)<4*spotDiameter
        Screen('DrawLine', windowPtr, green, screenXpixels*1/3-spotRadius , screenYpixels-3*spotRadius-elevationL, screenXpixels*2/3+spotRadius, screenYpixels-3*spotRadius-elevationR ,10);
        % Flip to the screen
        vbl  = Screen('Flip', windowPtr, vbl + (waitframes -0.5) * ifi);
    else
        Screen('DrawLine', windowPtr, red, screenXpixels*1/3+4*spotDiameter , screenYpixels-3*spotRadius-elevationL, screenXpixels*2/3-4*spotDiameter, screenYpixels-3*spotRadius-elevationR ,10);
        textMissionFailed=['Mission Failed!'];
        DrawFormattedText2(textMissionFailed,'win',windowPtr,...
        'sx','center','sy', 'center','xalign','center','yalign','top','baseColor',white);
        % Flip to the screen
        vbl  = Screen('Flip', windowPtr, vbl + (waitframes -0.5) * ifi);
        WaitSecs(1)
        break;
    end

    
    % left player
    try
        tic
        % Left player
        [pressedL1, RBkeyL1]=readCedrusRB(deviceL, keymapL);
        [pressedL2, RBkeyL2]=readCedrusRB(deviceL, keymapL);
        toc
    catch
        % do nothing; % if no nothing is pressed or released continue to next loop without throwing error 
    end
    
    % Left player
    if pressedL1 ==1 | pressedL2 == 1 % key press detected
        if RBkeyL1 == 4 | RBkeyL2 == 4
            if elevationL < screenYpixels - 3*spotDiameter
                elevationL = elevationL + spotDiameter;
            else
%                 elevationL = elevationL;
            end
        end
    end
    % Right player
    try
        tic
        % Right player
        [pressedR1, RBkeyR1]=readCedrusRB(deviceR, keymapR);
        [pressedR2, RBkeyR2]=readCedrusRB(deviceR, keymapR);
        toc
    catch
        % do nothing; 
    end
    % Right player
    if pressedR1 ==1 | pressedR2 == 1 % key press detected
        if RBkeyR1 == 2 | RBkeyR2 == 2
            if elevationR < screenYpixels - 3*spotDiameter
                elevationR = elevationR + spotDiameter;
            else
%                 elevationR = elevationR;
            end
        end
    end
    % reset
    pressedL1=[]; RBkeyL1=[]; %portL1=[]; stampL1=[];
    pressedL2=[]; RBkeyL2=[]; %portL2=[]; stampL2=[];
    pressedR1=[]; RBkeyR1=[]; %portR1=[]; stampR1=[];
    pressedR2=[]; RBkeyR2=[]; %portR2=[]; stampR2=[];
    

end

        