 
while n < numFrames
    % If esc is press, break out of the while loop and close the screen
    [keyIsDown, keysecs, keyCode] = KbCheck;
    if keyCode(KbName('escape'))
        Screen('CloseAll');
        break;
    end

    Screen('DrawTexture', windowPtr, textureIndex, [],posL);
    Screen('DrawTexture', windowPtr, textureIndex, [],posR);

    % Flip to the screen
    vbl  = Screen('Flip', windowPtr, vbl + (waitframes -0.5) * ifi);

    % update n
    n = n+1;
end
