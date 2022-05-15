% time for resting EEG (EO=eye open; EC= eye close)
TimeRestingEEG=30; % 30 seconds @@@@@@@@@@@@@@@@@@@@@@@@@@@@2 
numFramesRestEye=round (TimeRestingEEG/(ifi*waitframes)); 

%*************************************************************************
% Start taking eye open and eye close resting stage EEG
instructionStart=['Eye open, rest for ' num2str(TimeRestingEEG/60) ' min' '\n Experimenter hit a key to start'];
DrawFormattedText2(instructionStart,'win',windowPtr,...
    'sx', 'center','sy', 'center','xalign','center','yalign','top','baseColor',white);
DrawFormattedText2(instructionStart,'win',windowPtr,...
    'sx', xCenterL,'sy', 'center','xalign','center','yalign','top','baseColor',white);
DrawFormattedText2(instructionStart,'win',windowPtr,...
    'sx', xCenterR,'sy', 'center','xalign','center','yalign','top','baseColor',white);
Screen('Flip',windowPtr);
% hit a key to continue
KbStrokeWait;
Screen('Flip',windowPtr);% Reponse with a black screen
pause(1); % given time to return the position of arm after key press 

% Flash photocells once at the start and end of eye open period 
% Flash photocells once at the start of eye open period ###### (upper left and bottom right)
Screen('FillRect', windowPtr, white, RightUpperSquare);  % event type = 1200001
Screen('FillRect', windowPtr, white, LeftUpperSquare);
Screen('Flip',windowPtr);



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
    % Show the fixation cross
    Screen('DrawDots', windowPtr, [FixCrX;FixCrY], screenYpixels/400, white, [0 0], 2);
    Screen('DrawDots', windowPtr, [FixCrX-screenXpixels/3;FixCrY], screenYpixels/400, white, [0 0], 2);
    Screen('DrawDots', windowPtr, [FixCrX+screenXpixels/3;FixCrY], screenYpixels/400, white, [0 0], 2);
    % Show the central dot
    % Screen('DrawDots', windowPtr, [xCenter;yCenter], Thickness, white, [0 0], 2);
    % Flip the black screen
    vbl  = Screen('Flip', windowPtr, vbl + (waitframes -0.5) * ifi);
    % update the while loop
    i=i+1;
end

% Flash photocells once to mark the end of open eye resting EEG collection ###### (upper left and bottom right)
Screen('FillRect', windowPtr, white, RightUpperSquare);  % event type = 1200001
Screen('FillRect', windowPtr, white, LeftUpperSquare);
Screen('Flip',windowPtr);

%****************************************************************************
instructionStart=['        Ready to continue?' '\n Experimenter hit a key to resume.']
DrawFormattedText2(instructionStart,'win',windowPtr,...
    'sx', 'center','sy', 'center','xalign','center','yalign','top','baseColor',white);
DrawFormattedText2(instructionStart,'win',windowPtr,...
    'sx', xCenterL,'sy', 'center','xalign','center','yalign','top','baseColor',white);
DrawFormattedText2(instructionStart,'win',windowPtr,...
    'sx', xCenterR,'sy', 'center','xalign','center','yalign','top','baseColor',white);
Screen('Flip',windowPtr);
% hit a key to continue
KbStrokeWait;
Screen('Flip',windowPtr);% Response with a black screen
pause(1); % given time to return the position of arm after key press 

%*************************************************************************
