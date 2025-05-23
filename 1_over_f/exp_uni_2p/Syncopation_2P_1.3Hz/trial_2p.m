% run SuppressWarning.m 
% Empty extra taps in the memory before the trial begin
i=1;
while i<100
    try
     [pressedL1, RBkeyL1]=readCedrusRB(deviceL, keymapL);% extract first key press  
     [pressedR1, RBkeyR1]=readCedrusRB(deviceR, keymapR);% extract first key press  
    catch
        % do nothing
    end
    i=i+1;
end

%initialize variables
pressedL1=[]; RBkeyL1=[]; 
pressedL2=[]; RBkeyL2=[]; 
pressedR1=[]; RBkeyR1=[];
pressedR2=[]; RBkeyR2=[]; 

% fliptimes=[];

% initial or reset trial frame number
n=1;

% initialized number of taps recorded
tapsRecordedL=0;tapsRecordedR=0;

numFrames = NumFramesInterval13Hz*600*3; % add numFrames to allow more time

MaxTaps; numTaps;

while (n < numFrames) & (tapsRecordedL < MaxTaps ) & (tapsRecordedR < MaxTaps ) % either one reach 600 taps  @@@@@@@@
    tic
    % If esc is pressed, break out of the while loop and close the screen
    [keyIsDown, keysecs, keyCode] = KbCheck;
    if keyCode(KbName('escape'));
        Screen('CloseAll');
        break;
    end
    % If F1 is pressed, exit and continue to the next condition 
    if keyCode(KbName('F1'));
        n=numFrames;
    end
        
    %% draw instruction for each trial / condition with different color
    % for experimentor
    DrawFormattedText2(ConditionInstructions{conditionSelected},'win',windowPtr,...
    'sx','center','sy', screenYpixels*0.8,'xalign','center','yalign','top','baseColor',color);
    % for player L
    if conditionSelected == 1;
        DrawFormattedText2(ConditionInstructionsL{conditionSelected},'win',windowPtr,...
    'sx',xCenterL,'sy', screenYpixels*0.8,'xalign','center','yalign','top','baseColor',red);
    elseif conditionSelected == 4;
        DrawFormattedText2(ConditionInstructionsL{conditionSelected},'win',windowPtr,...
    'sx',xCenterL,'sy', screenYpixels*0.8,'xalign','center','yalign','top','baseColor',blue);
    else
        DrawFormattedText2(ConditionInstructionsL{conditionSelected},'win',windowPtr,...
    'sx',xCenterL,'sy', screenYpixels*0.8,'xalign','center','yalign','top','baseColor',color);
    end
    % for player R
    if conditionSelected == 1;
        DrawFormattedText2(ConditionInstructionsR{conditionSelected},'win',windowPtr,...
    'sx',xCenterR,'sy', screenYpixels*0.8,'xalign','center','yalign','top','baseColor',blue);
    elseif conditionSelected == 4;
        DrawFormattedText2(ConditionInstructionsR{conditionSelected},'win',windowPtr,...
    'sx',xCenterR,'sy', screenYpixels*0.8,'xalign','center','yalign','top','baseColor',red);
    else
        DrawFormattedText2(ConditionInstructionsR{conditionSelected},'win',windowPtr,...
    'sx',xCenterR,'sy', screenYpixels*0.8,'xalign','center','yalign','top','baseColor',color);
    end

    %% draw fixation crosses on 3 monitors (different color in each condition)
    Screen('DrawDots', windowPtr, [FixCrX;FixCrY], screenYpixels/400, color, [0 0], 2); % center monitor
    Screen('DrawDots', windowPtr, [FixCrX-screenXpixels/3;FixCrY], screenYpixels/400, color, [0 0], 2); % Left monitor
    Screen('DrawDots', windowPtr, [FixCrX+screenXpixels/3;FixCrY], screenYpixels/400, color, [0 0], 2); % Right monitor

    %% draw labels for feedback and stimulus only in the center moniter
    % stimulus label
%     DrawFormattedText2('pacer','win',windowPtr,...
%     'sx',xCenter,'sy',yCenter-screenYpixels/15,'xalign','center','yalign','top','baseColor',color);
%     DrawFormattedText2('L','win',windowPtr,...
%     'sx',xCenter-screenYpixels/4,'sy',yCenter-screenYpixels/15,'xalign','center','yalign','top','baseColor',color);
%     DrawFormattedText2('R','win',windowPtr,...
%     'sx',xCenter+screenYpixels/4,'sy',yCenter-screenYpixels/15,'xalign','center','yalign','top','baseColor',color);

    %% Flash photocell every other second during the whole trial (bottom right)
    if (~isempty(find(([1:2:numFrames]==n)))) % every two frames
    Screen('FillRect', windowPtr, white, RightBottomSquare+[-screenXpixels/3*2 0 -screenXpixels/3*2 0]); % Left monitor (Right Bottom)
    Screen('FillRect', windowPtr, white, RightBottomSquare+[-screenXpixels/3 0 -screenXpixels/3 0]); % Middle monitor (Right Bottom)
    Screen('FillRect', windowPtr, white, RightBottomSquare+[-screenXpixels/3+PhotosensorSize*2 0 -screenXpixels/3+PhotosensorSize*2 0]); % Right monitor (Left Bottom)
    end
    
    %% show resampled recording as pacers and photocells on the top at the same time
    
    % show one's own resampled recording
    if (conditionSelected == 1) ;
        if any(ShowframesL(:) == n); 
            Screen('DrawDots', windowPtr, [xCenter-screenXpixels/3;yCenter], screenYpixels/20-2, red, [0 0], 2); % Left monitor
            Screen('FillRect', windowPtr, white, LeftUpperSquare); % photocell
        end
        if any(ShowframesR(:) == n); 
            Screen('DrawDots', windowPtr, [xCenter+screenXpixels/3;yCenter], screenYpixels/20-2, blue, [0 0], 2); % Right monitor
            Screen('FillRect', windowPtr, white, RightUpperSquare); % photocell
        end
    end     
    if (conditionSelected == 2) ;
        if any(ShowframesL(:) == n); 
            Screen('DrawDots', windowPtr, [xCenter-screenXpixels/3;yCenter], screenYpixels/20-2, red, [0 0], 2); % Left monitor
            Screen('FillRect', windowPtr, white, LeftUpperSquare); % photocell
            Screen('FillRect', windowPtr, white, RightUpperSquare); % photocell
        end
    end
    if (conditionSelected == 3) ;
        if any(ShowframesR(:) == n); 
            Screen('DrawDots', windowPtr, [xCenter+screenXpixels/3;yCenter], screenYpixels/20-2, blue, [0 0], 2); % Right monitor
            Screen('FillRect', windowPtr, white, LeftUpperSquare); % photocell
            Screen('FillRect', windowPtr, white, RightUpperSquare); % photocell
        end
    end
    if (conditionSelected == 4) ;
        if any(ShowframesL(1:30) == n); % Pace right subject to get started
            Screen('DrawDots', windowPtr, [xCenter+screenXpixels/3;yCenter], screenYpixels/20-2, red, [0 0], 2); % Right monitor
            Screen('FillRect', windowPtr, white, LeftUpperSquare); % photocell
            Screen('FillRect', windowPtr, white, RightUpperSquare); % photocell
        end
    end
    
    %% Left Player
    % Left player  (read button press)
    try
        % Left player
        [pressedL1, RBkeyL1]=readCedrusRB(deviceL, keymapL); % extract first key press  
    catch
        % do nothing; % if no nothing is pressed or released continue to next loop without throwing error 
    end
    % read empty buffer and return function is time consuming
    % Left player (leader, providing feedback)
    if pressedL1 ==1 %| pressedL2 == 1 % at least one key press detected in the frist two events of the previous buffer
        % if RBkeyL1 == 3 %| RBkeyL2 == 3 % confirm it is the middle key on RB740
            % For the experimentor (feedback and its photocell)
            Screen('DrawDots', windowPtr, [xCenter-screenXpixels/24;yCenter], screenYpixels/20-2, red, [0 0], 2); % center monitor (left shift -screenYpixels/24)
            % show bottomright photocell on the other side undercovered
            Screen('FillRect', windowPtr, white, RightBottomSquare); 
            % show upperleft photocell on the other side undercovered - right monitor
            Screen('FillRect', windowPtr, white, RightUpperSquare+[-screenXpixels/3+PhotosensorSize*2 0 -screenXpixels/3+PhotosensorSize*2 0]); 

            % feedback to the other paticipants
            if (conditionSelected == 2) | (conditionSelected == 4) % show feedback to other
                Screen('DrawDots', windowPtr, [xCenter+screenXpixels/3;yCenter], screenYpixels/20-2, red, [0 0], 2); % Right monitor
            end

            % update number of taps recorded
            tapsRecordedL=tapsRecordedL+1;
        % end
    end
    
    %% Right Player
    % Right player  (read button press)
    try
        % Right player
        [pressedR1, RBkeyR1]=readCedrusRB(deviceR, keymapR);
    catch
        % do nothing; 
    end

    % Right player (leader, providing feedback)
    if pressedR1 ==1 %| pressedR2 == 1  % at least one key press detected in the frist two events of the previous buffer
        % if RBkeyR1 == 3 %| RBkeyR2 == 3  % confirm it is the middle key on RB740
            % show feedback for the experimentor 
            % Screen('DrawDots', windowPtr, [xCenter;yCenter], screenYpixels/20-2, blue, [0 0], 2); % center monitor 
            Screen('DrawDots', windowPtr, [xCenter+screenXpixels/24;yCenter], screenYpixels/20-2, blue, [0 0], 2); % center monitor (right shift +screenYpixels/24)
            % show bottomleft photocell on the other side undercovered
            Screen('FillRect', windowPtr, white, LeftBottomSquare); 
            % show upperleft photocell on the other side undercovered - left monitor
            Screen('FillRect', windowPtr, white, LeftUpperSquare+[+screenXpixels/3-PhotosensorSize*2 0 +screenXpixels/3-PhotosensorSize*2 0]); 
          
            % provide feedback stim to the other player 
            if (conditionSelected == 3) | (conditionSelected == 4)
                Screen('DrawDots', windowPtr, [xCenter-screenXpixels/3;yCenter], screenYpixels/20-2, blue, [0 0], 2); % Left monitor
            end
           
            % update number of taps recorded
            tapsRecordedR=tapsRecordedR+1;
        % end
    end

    %% Flip to the screen
%     fliptime=toc;
%     fliptimes=[fliptimes fliptime];
    vbl  = Screen('Flip', windowPtr, vbl + (waitframes -0.5) * ifi);
%     Screen('Flip', windowPtr);

    %% reset
    pressedL1=[]; RBkeyL1=[]; %portL1=[]; stampL1=[];
    pressedL2=[]; RBkeyL2=[]; %portL2=[]; stampL2=[];
    pressedR1=[]; RBkeyR1=[]; %portR1=[]; stampR1=[];
    pressedR2=[]; RBkeyR2=[]; %portR2=[]; stampR2=[];

    % update n
    n = n+1;
    % fliptime=toc;
    % fliptimes=[fliptimes fliptime];
end