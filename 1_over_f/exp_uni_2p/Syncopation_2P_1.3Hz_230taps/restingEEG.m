% time for resting EEG (EO=eye open; EC= eye close)
TimeRestingEEG=150; % 2.5 min = 150 seconds @@@@@@@@@@@@@@@@@@@@@@@@@2
numFramesRestEye=round (TimeRestingEEG/(ifi*waitframes)); 

%*************************************************************************
% Start taking eye open and eye close resting stage EEG
instructionStart=['Rest with eye open for ' num2str(TimeRestingEEG/60) ' min'];
DrawFormattedText2(instructionStart,'win',windowPtr,...
    'sx', 'center','sy', 'center','xalign','center','yalign','top','baseColor',white);
DrawFormattedText2(instructionStart,'win',windowPtr,...
    'sx', xCenterL,'sy', 'center','xalign','center','yalign','top','baseColor',white);
DrawFormattedText2(instructionStart,'win',windowPtr,...
    'sx', xCenterR,'sy', 'center','xalign','center','yalign','top','baseColor',white);
Screen('TextSize',windowPtr, 10);% reset font size to default
instructionStart=['Experimenter hit a key to start.'];
DrawFormattedText2(instructionStart,'win',windowPtr,...
    'sx', 'center','sy', screenYpixels*0.8,'xalign','center','yalign','top','baseColor',white);
DrawFormattedText2(instructionStart,'win',windowPtr,...
    'sx', xCenterL,'sy', screenYpixels*0.8,'xalign','center','yalign','top','baseColor',white);
DrawFormattedText2(instructionStart,'win',windowPtr,...
    'sx', xCenterR,'sy', screenYpixels*0.8,'xalign','center','yalign','top','baseColor',white);
Screen('Flip',windowPtr);
Screen('TextSize',windowPtr, 18);% reset font size to default
% hit a key to continue
KbStrokeWait;
Screen('Flip',windowPtr);% Reponse with a black screen
pause(1); % given time to return the position of arm after key press 

% Flash photocells once at the start and end of eye open period 
% Flash photocells once at the start of eye open period ###### (upper left and bottom right)
Screen('FillRect', windowPtr, white, RightUpperSquare);  % event type = 1200001
Screen('FillRect', windowPtr, white, LeftUpperSquare);
Screen('Flip',windowPtr);
pause(0.02);

% get a timestamp and begin taking resting EEG ****************************
vbl = Screen('Flip', windowPtr);
i=1;
while i<numFramesRestEye
% If esc is press, break out of the while loop and close the screen
    [keyIsDown, keysecs, keyCode] = KbCheck;
    if keyCode(KbName('escape'))
        Screen('CloseAll');
        break;
    end
    % If F1 is pressed, exit and continue to the next condition
    if keyCode(KbName('F1'));
        i=numFramesRestEye;
    end
    % Show the fixation cross
    Screen('DrawDots', windowPtr, [FixCrX;FixCrY], screenYpixels/400, white, [0 0], 2);
    Screen('DrawDots', windowPtr, [FixCrX-screenXpixels/3;FixCrY], screenYpixels/400, white, [0 0], 2);
    Screen('DrawDots', windowPtr, [FixCrX+screenXpixels/3;FixCrY], screenYpixels/400, white, [0 0], 2);
    % Show the central dot
    % Screen('DrawDots', windowPtr, [xCenter;yCenter], Thickness, white, [0 0], 2);
    % Show reminder
    Screen('TextSize',windowPtr, 18);% reset font size to default
    instructionStart=['Look at fixation cross and relax'];
    DrawFormattedText2(instructionStart,'win',windowPtr,...
        'sx', 'center','sy', screenYpixels*0.8,'xalign','center','yalign','top','baseColor',white);
    DrawFormattedText2(instructionStart,'win',windowPtr,...
        'sx', xCenterL,'sy', screenYpixels*0.8,'xalign','center','yalign','top','baseColor',white);
    DrawFormattedText2(instructionStart,'win',windowPtr,...
        'sx', xCenterR,'sy', screenYpixels*0.8,'xalign','center','yalign','top','baseColor',white);
    % Flip the black screen
    vbl  = Screen('Flip', windowPtr, vbl + (waitframes -0.5) * ifi);
    % update the while loop
    i=i+1;
end

% Flash photocells once to mark the end of open eye resting EEG collection ###### (upper left and bottom right)
Screen('FillRect', windowPtr, white, RightUpperSquare);  % event type = 1200001
Screen('FillRect', windowPtr, white, LeftUpperSquare);
Screen('Flip',windowPtr);
pause(0.02); % enable left monitor to detect it.

%****************************************************************************
instructionStart=['Eyes closed for ' num2str(TimeRestingEEG/60) 'min'...
    '\nOpen your eyes when being tapped on the shoulder'...
    '\nor when you hear a beep sound.'];
DrawFormattedText2(instructionStart,'win',windowPtr,...
    'sx', 'center','sy', 'center','xalign','center','yalign','top','baseColor',white);
DrawFormattedText2(instructionStart,'win',windowPtr,...
    'sx', xCenterL,'sy', 'center','xalign','center','yalign','top','baseColor',white);
DrawFormattedText2(instructionStart,'win',windowPtr,...
    'sx', xCenterR,'sy', 'center','xalign','center','yalign','top','baseColor',white);
Screen('TextSize',windowPtr, 10); % reset font size
instructionStart=['Experimenter hit a key to continue.'];
DrawFormattedText2(instructionStart,'win',windowPtr,...
    'sx', 'center','sy', screenYpixels*0.8,'xalign','center','yalign','top','baseColor',white);
DrawFormattedText2(instructionStart,'win',windowPtr,...
    'sx', xCenterL,'sy', screenYpixels*0.8,'xalign','center','yalign','top','baseColor',white);
DrawFormattedText2(instructionStart,'win',windowPtr,...
    'sx', xCenterR,'sy', screenYpixels*0.8,'xalign','center','yalign','top','baseColor',white);
Screen('Flip',windowPtr);
Screen('TextSize',windowPtr, 18); % reset font size to default
% hit a key to continue
KbStrokeWait;
Screen('Flip',windowPtr);% Response with a black screen
pause(1); % given time to return the position of arm after key press 

% Flash photocells to mark the start and end of eye closed period 
% Flash once to start ###### (upper left and bottom right)
Screen('FillRect', windowPtr, white, RightUpperSquare);  % event type = 1200001
Screen('FillRect', windowPtr, white, LeftUpperSquare);
Screen('Flip',windowPtr);
pause(0.02); % enable left monitor to detect it.

% get a timestamp and begin taking resting EEG
vbl = Screen('Flip', windowPtr);
i=1;
while i<numFramesRestEye
    % If esc is press, break out of the while loop and close the screen
    [keyIsDown, keysecs, keyCode] = KbCheck;
    if keyCode(KbName('escape'))
        Screen('CloseAll');
        break;
    end
    % If F1 is pressed, exit and continue to the next condition
    if keyCode(KbName('F1'));
        i=numFramesRestEye;
    end
    % Show reminder
    Screen('TextSize',windowPtr, 18);% reset font size to default
    instructionStart=['No peeking. Eyes closed please!'];
    DrawFormattedText2(instructionStart,'win',windowPtr,...
        'sx', 'center','sy', screenYpixels*0.8,'xalign','center','yalign','top','baseColor',white);
    DrawFormattedText2(instructionStart,'win',windowPtr,...
        'sx', xCenterL,'sy', screenYpixels*0.8,'xalign','center','yalign','top','baseColor',white);
    DrawFormattedText2(instructionStart,'win',windowPtr,...
        'sx', xCenterR,'sy', screenYpixels*0.8,'xalign','center','yalign','top','baseColor',white);
    % Flip the black screen
    vbl  = Screen('Flip', windowPtr, vbl + (waitframes -0.5) * ifi);
    % update the while loop
    i=i+1;
end
    

% Flash photocells again to end ###### (upper left and bottom right)
i=1;
while i < 2;
    Screen('FillRect', windowPtr, white, RightUpperSquare);  % event type = 1200001
    Screen('FillRect', windowPtr, white, LeftUpperSquare);
    Screen('Flip',windowPtr);
    pause(0.02); % enable left monitor to detect it.
    i=i+1;
end


% ************************************************************************
% show the end
instructionStart=['We will now start the experiment!'] % Tell subject to open eye and start
DrawFormattedText2(instructionStart,'win',windowPtr,...
    'sx', 'center','sy', 'center','xalign','center','yalign','top','baseColor',white);
DrawFormattedText2(instructionStart,'win',windowPtr,...
    'sx', xCenterL,'sy', 'center','xalign','center','yalign','top','baseColor',white);
DrawFormattedText2(instructionStart,'win',windowPtr,...
    'sx', xCenterR,'sy', 'center','xalign','center','yalign','top','baseColor',white);
Screen('TextSize',windowPtr, 10); % font size
instructionStart=['Experimenter hit a key to continue.'];
DrawFormattedText2(instructionStart,'win',windowPtr,...
    'sx', 'center','sy', screenYpixels*0.8,'xalign','center','yalign','top','baseColor',white);
DrawFormattedText2(instructionStart,'win',windowPtr,...
    'sx', xCenterL,'sy', screenYpixels*0.8,'xalign','center','yalign','top','baseColor',white);
DrawFormattedText2(instructionStart,'win',windowPtr,...
    'sx', xCenterR,'sy', screenYpixels*0.8,'xalign','center','yalign','top','baseColor',white);
Screen('Flip',windowPtr);
Screen('TextSize',windowPtr, 18); % reset font size to default
% play a beep sound to signal subject to open eyes
Beeper(400,0.9,1.5);

% hit a key to continue
KbStrokeWait;
Screen('Flip',windowPtr);% Response with a black screen
pause(1); % given time to return the position of arm after key press 

%*************************************************************************
