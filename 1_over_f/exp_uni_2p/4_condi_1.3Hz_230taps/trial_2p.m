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

while (n < numFrames) & (tapsRecordedL < numTaps ) & (tapsRecordedR < numTaps ) % either one reach 600 taps  @@@@@@@@
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
    DrawFormattedText2([ConditionInstructions{conditionSelected} ' -- block ' num2str(block) ' trial ' num2str(t)],'win',windowPtr,...
    'sx','center','sy', screenYpixels*0.8,'xalign','center','yalign','top','baseColor',color);
    % for player L
    DrawFormattedText2(ConditionInstructionsL{conditionSelected},'win',windowPtr,...
    'sx',xCenterL,'sy', screenYpixels*0.8,'xalign','center','yalign','top','baseColor',color);
    % for player R
    DrawFormattedText2(ConditionInstructionsR{conditionSelected},'win',windowPtr,...
    'sx',xCenterR,'sy', screenYpixels*0.8,'xalign','center','yalign','top','baseColor',color);

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
    
    %% show stimulus(pacer) and photocells on the top at the same time
    if any(Showframeselected(:) == n)
        % Show the pacer on center monitor
        Screen('FrameOval', windowPtr,green, [xCenter-screenXpixels/24-screenYpixels/35 yCenter-screenYpixels/35 xCenter-screenXpixels/24+screenYpixels/35 yCenter+screenYpixels/35],1,1); % green pacer center monitor left
        Screen('FrameOval', windowPtr,green, [xCenter+screenXpixels/24-screenYpixels/35 yCenter-screenYpixels/35 xCenter+screenXpixels/24+screenYpixels/35 yCenter+screenYpixels/35],1,1); % green pacer center monitor right
        % Show pacer on the left monitor
         if (conditionSelected == 1) | (conditionSelected == 2) | (conditionSelected == 4)| (conditionSelected == 5) % shown on left monitor
            % Green outer circle
            Screen('FrameOval', windowPtr,green, [xCenter-screenXpixels/3-screenYpixels/35 yCenter-screenYpixels/35 xCenter-screenXpixels/3+screenYpixels/35 yCenter+screenYpixels/35],1,1);
            % Red dot to pace oneself on the left monitor (this might be confusing)
            % Screen('DrawDots', windowPtr, [xCenter-screenXpixels/3;yCenter], screenYpixels/20-2, red, [0 0], 2); % Left monitor
         end
        % Show pacer on the right monitor 
         if (conditionSelected == 1) | (conditionSelected == 3) | (conditionSelected == 4) | (conditionSelected == 5) % shown on right monitor
            % Green outer circle
            Screen('FrameOval', windowPtr,green, [xCenter+screenXpixels/3-screenYpixels/35 yCenter-screenYpixels/35 xCenter+screenXpixels/3+screenYpixels/35 yCenter+screenYpixels/35],1,1);
            % Blue dot to pace oneself on the left monitor (this might be confusing)
            % Screen('DrawDots', windowPtr, [xCenter+screenXpixels/3;yCenter], screenYpixels/20-2, blue, [0 0], 2); % Left monitor
         end
        % Flash the photocell to indicate the pacers
        Screen('FillRect', windowPtr, white, LeftUpperSquare);
        Screen('FillRect', windowPtr, white, RightUpperSquare);    
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
            if RBkeyL1==6;
                % show warning dot if last key was pressed
                Screen('DrawDots', windowPtr, [xCenter-screenXpixels/24;yCenter], screenYpixels/20-2, yellow, [0 0], 2); % center monitor (left shift -screenYpixels/24)
            else
                Screen('DrawDots', windowPtr, [xCenter-screenXpixels/24;yCenter], screenYpixels/20-2, red, [0 0], 2); % center monitor (left shift -screenYpixels/24)
            end
            % show bottomright photocell on the other side undercovered - right monitor
            Screen('FillRect', windowPtr, white, RightBottomSquare); 
            % show upperleft photocell on the other side undercovered - right monitor
            Screen('FillRect', windowPtr, white, RightUpperSquare+[-screenXpixels/3+PhotosensorSize*2 0 -screenXpixels/3+PhotosensorSize*2 0]); 

            % show leader player's own buttom press
            if (conditionSelected == 1) | (conditionSelected == 2) 
                Screen('DrawDots', windowPtr, [xCenter-screenXpixels/3;yCenter], screenYpixels/20-2, red, [0 0], 2); % Left monitor
            end

            % feedback to the other paticipants
            if (conditionSelected == 2) % show feedback
                Screen('DrawDots', windowPtr, [xCenter+screenXpixels/3;yCenter], screenYpixels/20-2, red, [0 0], 2); % Right monitor
            end

            if (conditionSelected == 4)
                if (n < NumFramesInterval13Hz*30) | (n == NumFramesInterval13Hz*30) % feed to oneself
                    Screen('DrawDots', windowPtr, [xCenter-screenXpixels/3;yCenter], screenYpixels/20-2, red, [0 0], 2); % Left monitor
                end
                if n > NumFramesInterval13Hz*30 % feed to the other
                    Screen('DrawDots', windowPtr, [xCenter+screenXpixels/3;yCenter], screenYpixels/20-2, red, [0 0], 2); % Right monitor
                end
            end

            if (conditionSelected == 5)
                if (n < NumFramesInterval13Hz*30) | (n == NumFramesInterval13Hz*30) % feed to oneself
                    Screen('DrawDots', windowPtr, [xCenter-screenXpixels/3;yCenter], screenYpixels/20-2, red, [0 0], 2); % Left monitor
                end
                if n > NumFramesInterval13Hz*30 % feed to the other
                    Screen('DrawDots', windowPtr, [xCenter+screenXpixels/3;yCenter], screenYpixels/20-2, red, [0 0], 2); % Right monitor
                end
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
        % [pressedR2, RBkeyR2]=readCedrusRB(deviceR, keymapR);
    catch
        % do nothing; 
    end

    % Right player (leader, providing feedback)
    if pressedR1 ==1 %| pressedR2 == 1  % at least one key press detected in the frist two events of the previous buffer
        % if RBkeyR1 == 3 %| RBkeyR2 == 3  % confirm it is the middle key on RB740
            % show feedback for the experimentor 
            if RBkeyR1==6;
                % Screen('DrawDots', windowPtr, [xCenter;yCenter], screenYpixels/20-2, blue, [0 0], 2); % center monitor 
                % show warning dot if last key was pressed
                Screen('DrawDots', windowPtr, [xCenter+screenXpixels/24;yCenter], screenYpixels/20-2, yellow, [0 0], 2); % center monitor (right shift +screenYpixels/24)
            else
                Screen('DrawDots', windowPtr, [xCenter+screenXpixels/24;yCenter], screenYpixels/20-2, blue, [0 0], 2); % center monitor (right shift +screenYpixels/24)
            end
            % show bottomleft photocell on the other side undercovered - left monitor
            Screen('FillRect', windowPtr, white, LeftBottomSquare);
            % show upperleft photocell on the other side undercovered - left monitor
            Screen('FillRect', windowPtr, white, LeftUpperSquare+[+screenXpixels/3-PhotosensorSize*2 0 +screenXpixels/3-PhotosensorSize*2 0]); 
            
            % show one's own buttom press
            if (conditionSelected == 1) | (conditionSelected == 3) 
                Screen('DrawDots', windowPtr, [xCenter+screenXpixels/3;yCenter], screenYpixels/20-2, blue, [0 0], 2); % Right monitor
            end

            % provide feedback stim to the other player 
            if (conditionSelected == 3)
                Screen('DrawDots', windowPtr, [xCenter-screenXpixels/3;yCenter], screenYpixels/20-2, blue, [0 0], 2); % Left monitor
            end
            
            if (conditionSelected == 4) 
                if (n < NumFramesInterval13Hz*30) | (n == NumFramesInterval13Hz*30) % feed to oneself
                    Screen('DrawDots', windowPtr, [xCenter+screenXpixels/3;yCenter], screenYpixels/20-2, blue, [0 0], 2); % Left monitor
                end
                if n > NumFramesInterval13Hz*30 % feed to the other
                    Screen('DrawDots', windowPtr, [xCenter-screenXpixels/3;yCenter], screenYpixels/20-2, blue, [0 0], 2); % Right monitor
                end
            end
            
            if (conditionSelected == 5) 
                if (n < NumFramesInterval13Hz*30) | (n == NumFramesInterval13Hz*30) % feed to oneself
                    Screen('DrawDots', windowPtr, [xCenter+screenXpixels/3;yCenter], screenYpixels/20-2, blue, [0 0], 2); % Left monitor
                end
                if n > NumFramesInterval13Hz*30 % feed to the other
                    Screen('DrawDots', windowPtr, [xCenter-screenXpixels/3;yCenter], screenYpixels/20-2, blue, [0 0], 2); % Right monitor
                end
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

% % Flash photocells once to mark the end of the trial ###### (upper left and bottom right)
% Screen('FillRect', windowPtr, white, RightUpperSquare);  % event type = 1200001
% Screen('FillRect', windowPtr, white, LeftUpperSquare);
% Screen('Flip',windowPtr);
% pause(0.02); % enable left monitor to detect it.