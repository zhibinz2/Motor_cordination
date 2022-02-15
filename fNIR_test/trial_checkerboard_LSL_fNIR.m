startTime = now;
numberOfSecondsElapsed = 0;

while numberOfSecondsElapsed < 10
    % If esc is press, break out of the while loop and close the screen
    [keyIsDown, keysecs, keyCode] = KbCheck;
    if keyCode(KbName('escape'))
        Screen('CloseAll');
        break;
    end

    % Update the while loop with time
    numberOfSecondsElapsed = round((now - startTime) * 10 ^ 5);
    
    if conditionSelected ==1
        if mod(n,2) == 1 % if n is odd number
            Screen('DrawTexture', windowPtr, textureIndex1, [],posL);
        end
        if mod(n,2) == 0 % if n is even number
        Screen('DrawTexture', windowPtr, textureIndex2, [],posL);
        end
    end

    if conditionSelected==2
        if mod(n,2) == 1 % if n is odd number
            Screen('DrawTexture', windowPtr, textureIndex1, [],posR);
        end
        if mod(n,2) == 0 % if n is even number
            Screen('DrawTexture', windowPtr, textureIndex2, [],posR);
        end
    end

    % Flip to the screen
    vbl  = Screen('Flip', windowPtr, vbl + (waitframes -0.5) * ifi);

    % LSL markers to the local network
    % send data into the outlet, sample by sample
%     outlet1.push_sample(1); % send data number 1
    % send markers into the outlet
    mrk = markers{1};
    outlet2.push_sample({mrk});   % note that the string is wrapped into a cell-array

    % update n
    n = n+1;
end
