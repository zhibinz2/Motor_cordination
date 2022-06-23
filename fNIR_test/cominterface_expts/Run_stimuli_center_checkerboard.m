% This version shows checkerboard on center
% Zhibin 05/05/2022

% Clear the screen display and all variables in the workspace
sca; clc; close all; clear all; clearvars; 

% set the random number seed as the date of today in formate such as 20220505
seed=input('enter the date in format YYYYMMDD:');
rng(seed);

% Break and issue an error message if the installed Psychtoolbox is not
% based on OpenGL or Screen() is not working properly.
AssertOpenGL;


% Set up cominterface
oxy = actxserver('oxysoft.oxyapplication')

%% ########################################################################
try      % if anything went wrong, exit the display and show the error on the Command window
    
    % Here we call some default settings for setting up Psychtoolbox 
    PsychDefaultSetup(2);

    % Start with black screen
    % Removes the blue screen flash and minimize extraneous warnings.
    Screen('Preference', 'VisualDebugLevel', 1); 
    Screen('Preference', 'SkipSyncTests', 1);
    Screen('Preference', 'SuppressAllWarnings', 1);

    % Get the screen numbers. This gives us a number for each of the screens
    % attached to our computer. For help see: Screen Screens?
    screens = Screen('Screens');
    % Draw we select the maximum of these numbers. So in a situation where we
    % have two screens attached to our monitor we will draw to the external
    % screen. When only one screen is attached to the monitor we will draw to
    % this. For help see: help max
    screenNumber = max(screens);

    % Define black and white (white will be 1 and black 0). This is because 
    % luminace values are (in general) defined between 0 and 1.0
    % For help see: help WhiteIndex and help BlackIndex
    white = WhiteIndex(screenNumber);
    black = BlackIndex(screenNumber);
    white = [white white white];
    black = [black black black];
    % Initialize some other colors
    green = [0 1 0];

    % Open an on screen window and color it black
    % For help see: Screen Openwindow?
    % This will draw on a black backgroud with a size of [0 0 800 600] or
    % full screen, and then return a window pointer windowPtr
    % [windowPtr, windowRect] = PsychImaging('Openwindow', screenNumber, black, [0 0 800 600]); 
    [windowPtr, windowRect] = PsychImaging('Openwindow', screenNumber, black); % show on a full screen

    % Get the size of the on screen windowPtr in pixels
    % For help see: Screen windowSize?
    [screenXpixels, screenYpixels] = Screen('windowSize', windowPtr);

    % Get the centre coordinate of the window in pixels
    % For help see: help RectCenter
    [xCenter, yCenter] = RectCenter(windowRect); 

    % Enable alpha blending for anti-aliasing
    % For help see: Screen BlendFunction?
    % Also see: Chapter 6 of the OpenGL programming guide
    Screen('BlendFunction', windowPtr, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    Screen('ColorRange', windowPtr, 1, -1,1);

    % Set text display options for operating system other than Linux.
    % For help see: Screen TextFont?
    if ~IsLinux
        Screen('TextFont', windowPtr, 'Arial');
        Screen('TextSize', windowPtr, 18);
    end

    % Retreive the maximum priority number
    topPriorityLevel = MaxPriority(windowPtr) ; 
    % topPriorityLevel0 = MaxPriority(windowPtr0); 
    % set Priority once at the start of a script after setting up onscreen window.
    Priority(topPriorityLevel);

    % Measure the vertical refresh rate of the monitor
    ifi = Screen('GetFlipInterval', windowPtr);
    % Check if ifi=0.0167 second (60Hz)
    if round(1/ifi)~=60
      error('Error: Screen flash frequency is not set at 60Hz.');    
    end

    %% Checkerboard Design*************************************************
    % Green and white squares in 2x2 tiles fit into 5x5 arrangement
    greenSquare=repmat(cat(3,0,1,0),10,10);
    whiteSquare=repmat(cat(3,1,1,1),10,10);
    QuadraSquare1=[greenSquare whiteSquare;whiteSquare greenSquare];
    QuadraSquare2=[whiteSquare greenSquare;greenSquare whiteSquare];
    K1=repmat(QuadraSquare1,5,5);
    K2=repmat(QuadraSquare2,5,5);
    % MakeTexture (two reversing checkerboard arrangements)
    textureIndex1=Screen('MakeTexture',windowPtr,K1);
    textureIndex2=Screen('MakeTexture',windowPtr,K2);
    % Position of the checkerboard
    pos = [xCenter-screenXpixels/4,yCenter-screenXpixels/4,xCenter+screenXpixels/4,yCenter+screenXpixels/4];

    %% Set size of the squares for photocell ##############################
    PhotosensorSize=60; 
    % Positions of the photocell at he bottom rihgt corner of the screen
    RightBottomSquare= [screenXpixels-PhotosensorSize screenYpixels-PhotosensorSize screenXpixels screenYpixels];

    %% Timing of each video frame to flip the checkboard or the black screen
    % Numer of frames to wait when specifying good timing. Note: the use of
    % wait frames is to show a generalisable coding. For example, by using
    % waitframes = 2 one would flip on every other frame. See the PTB
    % documentation for details. In what follows we flip every frame. 
    % In order to reverse the checkerboard at 8Hz
    % We have to reverse the checkerboard every "(1/8 seconds)/ifi" video frames
    waitframes = (1/8)/ifi;

    %% Hide the cursor of the mouse in the display*****************
    % Get handles for all virtual pointing devices, aka cursors:
    typeOnly='masterPointer'; 
    mice = GetMouseIndices(typeOnly);  
    HideCursor(windowPtr,mice);

    %% Setting time variables for each stimulus**********************
    % total number of video frames to flip in each stimlus
    StimulusDuration=10; % 10 seconds 
    numFramesStim=round(StimulusDuration/ifi/waitframes);

    %% Tell subject to start
    instructionStart=['Each stimulus will show up after a resting period of 20-30 seconds with a black screen.'...
        '\nJust stay still and look at the screen for about 24 minutes.' ...
        ''...
        '\nNow. Press any key to start!']; 
    DrawFormattedText2(instructionStart,'win',windowPtr,...
        'sx','center','sy','center','xalign','center','yalign','center','baseColor',white);
    Screen('Flip',windowPtr);
    % hit a key to continue
    KbStrokeWait;
    % Flip Black Screen
    Screen('Flip',windowPtr);

    %% initialize some variables before entering the loops of blocks
    numBlock=40; % Do 40 stimuli
    RestDurations=[]; % Record the durations of resting intervals
    StimDurations=[]; % Record the durations of stimuli
    
    %% Entering the loop for each block
    for block=1:numBlock
        % Breakout the experiment by hitting esc
        [keyIsDown, keysecs, keyCode] = KbCheck;  
        if keyCode(KbName('escape'))
        Screen('CloseAll');
        break;
        end

        %% Resting with black screen before each stimulus ******************************* 
        % Set the time for the resting interval before/after each stimulus with randomized jittering  
        restLength = 25 + 10.*(rand(1)-0.5); % uniform distribution in seconds, with a mean of 25 s (20-30s)
        numFramesRest = round(restLength / ifi / waitframes);

        % flash photocell once at the Bottom right cornor on the screen to mark the beginning of resting
        % Screen('FillRect', windowPtr, white, RightBottomSquare);
        % Screen('Flip', windowPtr);
    
        % get a timestamp by flip a black screen at the start of resting 
        vbl = Screen('Flip', windowPtr);
		
        tic
        % initalize the while loop for displaying black screen
        n=1;
        while n < numFramesRest
            % If esc is press, break out of the while loop and exit experiment
            [keyIsDown, keysecs, keyCode] = KbCheck;
            if keyCode(KbName('escape'))
                Screen('CloseAll');
                break;
            end
            % Flip the black screen
            vbl  = Screen('Flip', windowPtr, vbl + (waitframes -0.5) * ifi);
            % Update the loop index n
            n=n+1;
        end
        % Display the duration in command window
        RestDuration=toc % check the time length in command window after experiment

        % flash the photocell again to mark the end of resting
        % Screen('FillRect', windowPtr, white, RightBottomSquare);
        % Screen('Flip', windowPtr);
        
        %% Show the stimulus
        % flash photocell once at the Bottom right cornor on the screen to
        % mark the beginning of stimulus
        % Screen('FillRect', windowPtr, white, RightBottomSquare);
        % Screen('Flip', windowPtr);
        % pause(0.3); 

          
        % get a timestamp by flip a black screen at the start of stimulus 
        vbl = Screen('Flip', windowPtr);
		oxy.WriteEvent('A','stim')
        tic
        % initalize the while loop for displaying black screen
        n=1;
        while n < numFramesStim
            % If esc is press, break out of the while loop and exit experiment
            [keyIsDown, keysecs, keyCode] = KbCheck;
            if keyCode(KbName('escape'))
                Screen('CloseAll');
                break;
            end

            % Draw two reversing checkerboard pattern alternatively
            if mod(n,2) == 1; % if n is odd number
                Screen('DrawTexture', windowPtr, textureIndex1, [],pos);
            end
            if mod(n,2) == 0; % if n is even number
                Screen('DrawTexture', windowPtr, textureIndex2, [],pos);
            end

            % Flip to the screen
            vbl  = Screen('Flip', windowPtr, vbl + (waitframes -0.5) * ifi);

            % Update the loop index n
            n=n+1;
        end
        
        oxy.WriteEvent('B','rest')
        % Display the duration in command window
        StimDuration=toc % check the time length in command window after experiment

        % flash the photocell again to mark the end of resting
        % Screen('FillRect', windowPtr, white, RightBottomSquare);
        % Screen('Flip', windowPtr);

        %% Save the time durations for checking after the experiment
        RestDurations=[RestDurations;RestDuration];
        StimDurations=[StimDurations;StimDuration];
        
    end

%% One last resting period after the last stimulus
% Set the time for the resting interval before/after each stimulus with randomized jittering  
restLength = 25 + 10.*(rand(1)-0.5); % uniform distribution in seconds, with a mean of 25 s (20-30s)
numFramesRest = round(restLength / ifi / waitframes);

% flash photocell once at the Bottom right cornor on the screen to mark the beginning of resting
% Screen('FillRect', windowPtr, white, RightBottomSquare);
% Screen('Flip', windowPtr);

% get a timestamp by flip a black screen at the start of resting 
vbl = Screen('Flip', windowPtr);
tic
% initalize the while loop for displaying black screen
n=1;
while n < numFramesRest
    % If esc is press, break out of the while loop and exit experiment
    [keyIsDown, keysecs, keyCode] = KbCheck;
    if keyCode(KbName('escape'))
        Screen('CloseAll');
        break;
    end
    % Flip the black screen
    vbl  = Screen('Flip', windowPtr, vbl + (waitframes -0.5) * ifi);
    % Update the loop index n
    n=n+1;
end
% Display the duration in command window
RestDuration=toc % check the time length in command window after experiment

% flash the photocell again to mark the end of resting
% Screen('FillRect', windowPtr, white, RightBottomSquare);
% 
%   
%   Screen('Flip', windowPtr);

% Save the time duration for checking after the experiment
RestDurations=[RestDurations;RestDuration];

%% Show The End
TheEnd = ['The End.' ...
    '\nHit a key to exit.'];
DrawFormattedText2(TheEnd,'win',windowPtr,...
    'sx','center','sy', 'center','xalign','center','yalign','top','baseColor',white);
vbl=Screen('Flip',windowPtr);

% hit a key to exit
KbStrokeWait;Screen('Flip',windowPtr);pause(1);
%*************************************
Priority(0);   
sca;

catch
    sca;
    psychrethrow(psychlasterror);
end  

%% Save stimuli data (all the variables created by the experiment in the Workspace)
filename=[num2str(seed) '.mat'];
save(filename);
 

 

















