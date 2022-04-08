% initial or reset trial frame number
n=1;
% get a timestamp and begin the trial
vbl = Screen('Flip', windowPtr);
            
% startTime = now;
% numberOfSecondsElapsed = 0;
tic 
% while numberOfSecondsElapsed < 10;
while n <  numFrames 
    % If esc is press, break out of the while loop and close the screen
    [keyIsDown, keysecs, keyCode] = KbCheck;
    if keyCode(KbName('escape'));
        Screen('CloseAll');
        break;
    end

    % Update the while loop with time
%     numberOfSecondsElapsed = (now - startTime) * 10 ^ 5;
    
    % photocell at beginning and end
    if n < 2 | n > numFrames-3; % flash for the frist 2 video frames (at least 75/2 Hz), because fNIRS sampled at 75Hz
        Screen('FillRect', windowPtr, white, RightUpperSquare);  % event type = 1200001
        Screen('FillRect', windowPtr, white, RightBottomSquare);
    end
     
    
    if conditionSelected ==1;
        if mod(n,2) == 1; % if n is odd number
            Screen('DrawTexture', windowPtr, textureIndex1, [],posL);
        end
        if mod(n,2) == 0; % if n is even number
            Screen('DrawTexture', windowPtr, textureIndex2, [],posL);
        end
    end

    if conditionSelected==2;
        if mod(n,2) == 1; % if n is odd number
            Screen('DrawTexture', windowPtr, textureIndex1, [],posR);
        end
        if mod(n,2) == 0; % if n is even number
            Screen('DrawTexture', windowPtr, textureIndex2, [],posR);
        end
    end
    
    % Show the fixation cross
    Screen('DrawDots', windowPtr, [FixCrX;FixCrY], screenXpixels/300, white, [0 0], 2);

    % Flip to the screen
    vbl  = Screen('Flip', windowPtr, vbl + (waitframes -0.5) * ifi);

%      % LSL marker to check screen flip frequency
%     % send markers into the outlet
%     mrk = markers{5};
%     outlet2.push_sample({mrk});   % note that the string is wrapped into a cell-array

    % update n
    n = n+1;
end
TrialDuration=toc
TrialDurations=[TrialDurations;TrialDuration];
