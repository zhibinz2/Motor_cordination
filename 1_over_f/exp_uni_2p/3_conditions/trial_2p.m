% run SuppressWarning.m 

%initialize variables
pressedL1=[]; RBkeyL1=[]; %portL1=[]; stampL1=[];
pressedL2=[]; RBkeyL2=[]; %portL2=[]; stampL2=[];
pressedR1=[]; RBkeyR1=[]; %portR1=[]; stampR1=[];
pressedR2=[]; RBkeyR2=[]; %portR2=[]; stampR2=[];
% k1=[];k2=[];
fliptimes=[];

% startTime = now;
% numberOfSecondsElapsed = 0;
% tic 
% while numberOfSecondsElapsed < 10;
% while n <  numFrames 
numFrames = max(Showframeselected)+NumFramesInterval; 

while n < numFrames
    tic
    % If esc is press, break out of the while loop and close the screen
    [keyIsDown, keysecs, keyCode] = KbCheck;
    if keyCode(KbName('escape'));
        Screen('CloseAll');
        break;
    end
        
    %% draw instruction for each trial / condition with different color
    DrawFormattedText2(conditionNames{conditionSelected},'win',windowPtr,...
    'sx','center','sy', screenYpixels*0.8,'xalign','center','yalign','top','baseColor',color);
    DrawFormattedText2(conditionNames{conditionSelected},'win',windowPtr,...
    'sx',xCenterL,'sy', screenYpixels*0.8,'xalign','center','yalign','top','baseColor',color);
    DrawFormattedText2(conditionNames{conditionSelected},'win',windowPtr,...
    'sx',xCenterR,'sy', screenYpixels*0.8,'xalign','center','yalign','top','baseColor',color);

    %% draw fixation crosses on 3 monitors (different color in each condition)
    Screen('DrawDots', windowPtr, [FixCrX;FixCrY], screenYpixels/400, color, [0 0], 2); % center monitor
    Screen('DrawDots', windowPtr, [FixCrX-screenXpixels/3;FixCrY], screenYpixels/400, color, [0 0], 2); % Left monitor
    Screen('DrawDots', windowPtr, [FixCrX+screenXpixels/3;FixCrY], screenYpixels/400, color, [0 0], 2); % Right monitor

    % Update the while loop with time
%     numberOfSecondsElapsed = (now - startTime) * 10 ^ 5;
     
%     if conditionSelected ==1;
%         if mod(n,2) == 1; % if n is odd number
%             Screen('DrawTexture', windowPtr, textureIndex1, [],posL);
%         end
%         if mod(n,2) == 0; % if n is even number
%         Screen('DrawTexture', windowPtr, textureIndex2, [],posL);
%         end
%     end
% 
%     if conditionSelected==2;
%         if mod(n,2) == 1; % if n is odd number
%             Screen('DrawTexture', windowPtr, textureIndex1, [],posR);
%         end
%         if mod(n,2) == 0; % if n is even number
%             Screen('DrawTexture', windowPtr, textureIndex2, [],posR);
%         end
%     end

    
    % *****OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    %% Flash photocell every other second during the whole trial (bottom left)
%     if (~isempty(find(([1:2:numFrames]==n)))) % every two frames
%     Screen('FillRect', windowPtr, white, LeftBottomSquare);
%     Screen('FillRect', windowPtr, white, RightBottomSquare); % event type = 1200002
%     end
    
    % ******#################################
    % Flash photocells at start&end of each trial/condition (upper left)
%     if n == 1 | n == numFrames
%     Screen('FillRect', windowPtr, white, LeftUpperSquare); % event type = 1200002
%     end

    %% show stimulus and photocells on the top at the same time
    if any(Showframeselected(:) == n)
        % Show the pacer stimulus
        Screen('DrawDots', windowPtr, [xCenter;yCenter], screenYpixels/30, color, [0 0], 2); % center monitor
        Screen('DrawDots', windowPtr, [xCenter-screenXpixels/3;yCenter], screenYpixels/30, color, [0 0], 2); % Left monitor
        Screen('DrawDots', windowPtr, [xCenter+screenXpixels/3;yCenter], screenYpixels/30, color, [0 0], 2); % Right monitor
        % Flash the photocell
        Screen('FillRect', windowPtr, white, LeftUpperSquare);
        Screen('FillRect', windowPtr, white, RightUpperSquare);  % event type = 1200001 for neuroScan
        
    end

    % show noise stimulus
    if conditionSelected ==3
        if any(Noiseframes3(:) == n)
            % Show the stimuli
            Screen('DrawDots', windowPtr, [xCenter;yCenter], screenYpixels/30, white, [0 0], 2); % center moniter
            Screen('DrawDots', windowPtr, [xCenter-screenXpixels/3;yCenter], screenYpixels/30, white, [0 0], 2); % Left
            Screen('DrawDots', windowPtr, [xCenter+screenXpixels/3;yCenter], screenYpixels/30, white, [0 0], 2); % Right
        end
    end

    %% Left Player
    % left player (read button press)
%     tic
%     k1 = read(deviceL,6,"char") % if no key pressed or release, k will be empty
%     k2 = read(deviceL,6,"char") % if k is not empty, then a least a key was pressed or relesed
%     toc
%     % read function takes about 30 ms if nothing in the buffer for the
%     % device, takes about 2 ms if a key pressed or release is in the
%     % buffer.
%     if  ~isempty(k1) | ~isempty(k2) % only when the key press and release both detected in the same loop
%         Screen('DrawDots', windowPtr, [xCenter;yCenter], screenXpixels/30, megenta, [0 0], 2); % center monitor
%         Screen('DrawDots', windowPtr, [xCenter+screenXpixels/3;yCenter], screenXpixels/30, megenta, [0 0], 2); % Right monitor
%         % flash photocell on the bottom at the buttom press
%         Screen('FillRect', windowPtr, white, LeftBottomSquare);
%         Screen('FillRect', windowPtr, white, RightBottomSquare); % event type = 1200002
%     end

    try
        % Left player
        [pressedL1, RBkeyL1]=readCedrusRB(deviceL, keymapL); % extract first key press
        % [pressedL2, RBkeyL2]=readCedrusRB(deviceL, keymapL); % extract first key release
    catch
        % do nothing; % if no nothing is pressed or released continue to next loop without throwing error 
    end
    % read empty buffer and return function is time consuming
    
    % Left player (provide feedback)
    if pressedL1 ==1 %| pressedL2 == 1 % at least one key press detected in the frist two events of the previous buffer
        if RBkeyL1 == 3 %| RBkeyL2 == 3 % confirm it is the middle key on RB740
            % show feedback for other paticipants
            Screen('DrawDots', windowPtr, [xCenter-screenYpixels/30;yCenter], screenYpixels/30, megenta, [0 0], 2); % center monitor
            Screen('DrawDots', windowPtr, [xCenter+screenXpixels/3-screenYpixels/30;yCenter], screenYpixels/30, megenta, [0 0], 2); % Right monitor
            % flash photocell on the bottom at the buttom press
            % Screen('FillRect', windowPtr, white, LeftBottomSquare);
            Screen('FillRect', windowPtr, white, RightBottomSquare); % event type = 1200002
        end
    end
    
    %% Right Player
    % Right player  (read button press)
    try
        % Right player
        [pressedR1, RBkeyR1]=readCedrusRB(deviceR, keymapR);
        % [pressedR2, RBkeyR2]=readCedrusRB(deviceR, keymapR);
    catch
        % do nothing; 
    end

    % Right player (provide feedback)
    if pressedR1 ==1 %| pressedR2 == 1  % at least one key press detected in the frist two events of the previous buffer
        if RBkeyR1 == 3 %| RBkeyR2 == 3  % confirm it is the middle key on RB740
            % show feedback for other paticipants
            Screen('DrawDots', windowPtr, [xCenter+screenYpixels/30;yCenter], screenYpixels/30, cyan, [0 0], 2); % center monitor
            Screen('DrawDots', windowPtr, [xCenter-screenXpixels/3+screenYpixels/30;yCenter], screenYpixels/30, cyan, [0 0], 2); % Right monitor
            % flash photocell on the bottom at the buttom press
            Screen('FillRect', windowPtr, white, LeftBottomSquare);
            % Screen('FillRect', windowPtr, white, RightBottomSquare); % event type = 1200002
        end
    end

    %% Flip to the screen
    vbl  = Screen('Flip', windowPtr, vbl + (waitframes -0.5) * ifi);

     % LSL marker to check screen flip frequency
    % send markers into the outlet
%     mrk = markers{5};
%     outlet2.push_sample({mrk});   % note that the string is wrapped into a cell-array

    %% reset
    pressedL1=[]; RBkeyL1=[]; %portL1=[]; stampL1=[];
    pressedL2=[]; RBkeyL2=[]; %portL2=[]; stampL2=[];
    pressedR1=[]; RBkeyR1=[]; %portR1=[]; stampR1=[];
    pressedR2=[]; RBkeyR2=[]; %portR2=[]; stampR2=[];
%     k1=[];k2=[];

    % update n
    n = n+1;
    fliptime=toc;
    fliptimes=[fliptimes fliptime];
end

% TrialDuration=toc
% TrialDurations=[TrialDurations;TrialDuration];