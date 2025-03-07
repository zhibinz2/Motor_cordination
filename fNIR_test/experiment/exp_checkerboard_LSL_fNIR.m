% This version shows checkerboard on two hemi-fields and sends LSL markers
% to the fNIRS laptop
sca; clc; close all; clear all; clearvars; 

% LSL #####################################################################
% Addpath
addpath(genpath('/home/hnl/Documents/GitHub/labstreaminglayer/LSL/liblsl-Matlab')); % for LSL
% instantiate the library for LSL
lib = lsl_loadlib();
% make new stream outlets
% info1 = lsl_streaminfo(lib,'BioSemi','EEG',8,100,'cf_float32','sdfwerr32432');
% the name (here MyMarkerStream) is visible to the experimenter and should be chosen so that 
% it is clearly recognizable as your MATLAB software's marker stream
% The content-type should be Markers by convention, and the next three arguments indicate the 
% data format (1 channel, irregular rate, string-formatted).
% The so-called source id is an optional string that allows for uniquely identifying your 
% marker stream across re-starts (or crashes) of your script (i.e., after a crash of your script 
% other programs could continue to record from the stream with only a minor interruption).
info2 = lsl_streaminfo(lib,'MyMarkerStream','Markers',1,0,'cf_string','myuniquesourceid23443');
% open the outlets
% outlet1 = lsl_outlet(info1);
outlet2 = lsl_outlet(info2);
% Set markers
% markers = {'Baseline-Start', 'Baseline-End', 'Trial-Start', 'Trial-End','Screen_Flip'};
markers = {'C', 'C', 'A', 'B','Screen_Flip'};

% ##################################################################### LSL


% Break and issue an error message if the installed Psychtoolbox is not
% based on OpenGL or Screen() is not working properly.
AssertOpenGL;

if ~IsLinux
  error('Sorry, this demo currently only works on a Linux.');    
end

% Set trial conditions ****************************************************
conditions = [1 2];
conditionNames={'L' 'R'};  % conditionNames{1}


% Block & Trial number of the experiment **********************************
% % number of trials per block
numTrials=4;
% number of blocks
numBlock=3; % 3 blocks per session,6 subject
% total trial number
numtotal=numTrials*numBlock; 
% num of conditions in the experiment
numconditions=length(conditions);
% **********************************Block & Trial number of the experiment

% Randomization of the experiment #########################################
% set the random number seed as the date of today in formate such as 20210809
seed=input('enter the date in format YYYYMMDD:');
rng(seed);
% how many semirandom permutation set in the experiment 
numPerm=numtotal/numconditions;
% create the whole set of random conditions for the experiment
allPerm=[];
for p=1:numPerm
    allPerm=[allPerm randperm(numconditions)];
end
% ######################################### Randomization of the experiment 
% ***************************************** Randomization of the experiment 



% ########################################################################
try      
    
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
    red   = [1 0 0];
    blue  = [0 0 1];
    green = [0 1 0];
    grey  = [0.5 0.5 0.5];
    yellow  = [1 1 0];
    megenta = [1 0 1];

    % Open an on screen window and color it black
    % For help see: Screen Openwindow?
    % This will draw on a black backgroud with a size of [0 0 500 1000] and
    % return a window pointer windowPtr
%     [windowPtr, windowRect] = PsychImaging('Openwindow', screenNumber, black, [0 0 600 400]); 
    [windowPtr, windowRect] = PsychImaging('Openwindow', screenNumber, black); 

    % Get the size of the on screen windowPtr in pixels
    % For help see: Screen windowSize?
    [screenXpixels, screenYpixels] = Screen('windowSize', windowPtr);

    % Get the centre coordinate of the window in pixels
    % For help see: help RectCenter
    [xCenter, yCenter] = RectCenter(windowRect); 

    % Checkerboard Design****************************************
    % Black and white
%     K1 = (checkerboard > 0.5);
%     K2 = ~K1;
    % Green and white
    greenSquare=repmat(cat(3,0,1,0),10,10);
    whiteSquare=repmat(cat(3,1,1,1),10,10);
    QuadraSquare1=[greenSquare whiteSquare;whiteSquare greenSquare];
    QuadraSquare2=[whiteSquare greenSquare;greenSquare whiteSquare];
    K1=repmat(QuadraSquare1,5,5);
    K2=repmat(QuadraSquare2,5,5);
    % MakeTexture
%     textureIndex1=Screen('MakeTexture',windowPtr,double(K1));
%     textureIndex2=Screen('MakeTexture',windowPtr,double(K2));
    textureIndex1=Screen('MakeTexture',windowPtr,K1);
    textureIndex2=Screen('MakeTexture',windowPtr,K2);
    % Position of the checkerboard
    xLcenter=xCenter/2;xRcenter=xCenter/2*3;Shift=screenXpixels/12;
    posL = [xLcenter-screenXpixels/20,yCenter-screenXpixels/20,xLcenter+screenXpixels/20,yCenter+screenXpixels/20];
    posL = posL + [-Shift +Shift*2 -Shift +Shift*2]; % shift to the lower left corner
    posR = [xRcenter-screenXpixels/20,yCenter-screenXpixels/20,xRcenter+screenXpixels/20,yCenter+screenXpixels/20];
    posR = posR + [+Shift +Shift*2 +Shift +Shift*2]; % shift to the lower right corner
    % ****************************************Checkerboard Design

    % Set size of the squares for photocell ###############################
    PhotosensorSize=30; % 30;
    % Positions of the four corners
    RightBottomSquare= [screenXpixels-PhotosensorSize*2 screenYpixels-PhotosensorSize*2 screenXpixels screenYpixels];
    RightUpperSquare= [screenXpixels-PhotosensorSize*2 0 screenXpixels PhotosensorSize*2];
    LeftBottomSquare= [0 screenYpixels-PhotosensorSize*2 PhotosensorSize*2 screenYpixels];
    LeftUpperSquare= [0 0 PhotosensorSize*2 PhotosensorSize*2];
%     RightBottomSquare= [screenXpixels-PhotosensorSize*2 screenYpixels/2+230-PhotosensorSize screenXpixels screenYpixels/2+230+PhotosensorSize];
%     RightUpperSquare= [screenXpixels-PhotosensorSize*2 screenYpixels/2+110-PhotosensorSize screenXpixels screenYpixels/2+110+PhotosensorSize];
%     LeftBottomSquare= [0 screenYpixels/2+230-PhotosensorSize PhotosensorSize*2 screenYpixels/2+230+PhotosensorSize];
%     LeftUpperSquare= [0 screenYpixels/2+110-PhotosensorSize PhotosensorSize*2 screenYpixels/2+110+PhotosensorSize];
    % ############################### Set size of the squares for photocell 

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
    % Check if ifi=0.0167
    if round(1/ifi)~=60
      error('Error: Screen flash frequency is not set at 60Hz.');    
    end

    % Numer of frames to wait when specifying good timing. Note: the use of
    % wait frames is to show a generalisable coding. For example, by using
    % waitframes = 2 one would flip on every other frame. See the PTB
    % documentation for details. In what follows we flip every frame. 
    % In order to reverse the checkerboard at 8Hz
    % We have to reverse the checkerboard every (1/8)/ifi frames
    waitframes = (1/8)/ifi;
    
    % Setting time variables**********************************************
    % Length of one minute baseline
    BaselineLength =60; % in seconds 4; % 4 s to move
    numFramesBaseline = round(BaselineLength / ifi / waitframes);

    % total number of frames per trial
    StimulusDuration=2; % in seconds  % 10 seconds 
    numFrames=round(StimulusDuration/ifi/waitframes);
    % **********************************************Setting time variables
    
    % Hide Mice:****************************************************************
    % Get handles for all virtual pointing devices, aka cursors:
    typeOnly='masterPointer'; 
    mice = GetMouseIndices(typeOnly);  
    % DuoMice: ***************************************************************
    % Hide the system-generated cursors. We do this, because only the
    % first mouse cursor is hardware-accelerated, HideCursorie., a GPU created
    % hardware cursor. All other cursors are software-cursors, created
    % by the Windowing system. These tend to flicker badly in our use
    % case. Therefore we disable all system cursor images and draw our
    % cursors ourselves for a more beautiful look:
    % Hide the cursor
    HideCursor(windowPtr,mice);
    %HideCursor(windowPtr,mice(1));

    % Baseline taking 60s  *******************************  % Checkerboard Design****************************************
    % Black and white
%     K1 = (checkerboard > 0.5);
%     K2 = ~K1;
    % Green and white
    greenSquare=repmat(cat(3,0,1,0),10,10);
    whiteSquare=repmat(cat(3,1,1,1),10,10);
    QuadraSquare1=[greenSquare whiteSquare;whiteSquare greenSquare];
    QuadraSquare2=[whiteSquare greenSquare;greenSquare whiteSquare];
    K1=repmat(QuadraSquare1,5,5);
    K2=repmat(QuadraSquare2,5,5);
    % MakeTexture
%     textureIndex1=Screen('MakeTexture',windowPtr,double(K1));
%     textureIndex2=Screen('MakeTexture',windowPtr,double(K2));
    textureIndex1=Screen('MakeTexture',windowPtr,K1);
    textureIndex2=Screen('MakeTexture',windowPtr,K2);
    % Position of the checkerboard
    xLcenter=xCenter/2;xRcenter=xCenter/2*3;Shift=screenXpixels/12;
    posL = [xLcenter-screenXpixels/12,yCenter-screenXpixels/12,xLcenter+screenXpixels/12,yCenter+screenXpixels/12];
    posL = posL + [-Shift +Shift*1.5 -Shift +Shift*1.5]; % shift to the lower left corner
    posR = [xRcenter-screenXpixels/12,yCenter-screenXpixels/12,xRcenter+screenXpixels/12,yCenter+screenXpixels/12];
    posR = posR + [+Shift +Shift*1.5 +Shift +Shift*1.5]; % shift to the lower right corner
    % ****************************************Checkerboard Design

    instructionStart=['Hit any key and then look at the center of the screen for ' num2str(BaselineLength) ' seconds'];
    DrawFormattedText2(instructionStart,'win',windowPtr,...
        'sx','center','sy','center','xalign','center','yalign','center','baseColor',white);
    Screen('Flip',windowPtr);
    % hit a key to continue
%     pause(0.5);
    KbStrokeWait;Screen('Flip',windowPtr);pause(2);
    % Create a fixation cross
    FixCrX=[xCenter-round(screenXpixels/500):xCenter+round(screenXpixels/500) repmat(xCenter,1,round(screenXpixels/250)+1)];
    FixCrY=[repmat(yCenter,1,round(screenXpixels/250)+1) yCenter-round(screenXpixels/500):yCenter+round(screenXpixels/500)];
    
    % LSL markers to the local network
    % send data into the outlet, sample by sample
%         outlet1.push_sample(1); % send data number 1
    % send markers into the outlet
    mrk = markers{1};
    outlet2.push_sample({mrk});   % note that the string is wrapped into a cell-array
    
    tic
%     startTime = now;
%     numberOfSecondsElapsed = 0;
%     while numberOfSecondsElapsed < 60

    % flash a photocell to mark the beginning of baseline taking
    Screen('FillRect', windowPtr, white, RightUpperSquare);  % event type = 1200001
    Screen('FillRect', windowPtr, white, RightBottomSquare);
    Screen('Flip', windowPtr);

    % get a timestamp and begin the baseline 
    vbl = Screen('Flip', windowPtr);
    n=1;
    while n < numFramesBaseline
        % If esc is press, break out of the while loop and close the screen
        [keyIsDown, keysecs, keyCode] = KbCheck;
        if keyCode(KbName('escape'))
            Screen('CloseAll');
            break;
        end

        % Update the while loop with time
%         numberOfSecondsElapsed = (now - startTime) * 10 ^ 5;
        
        % Show the fixation cross
        Screen('DrawDots', windowPtr, [FixCrX;FixCrY], screenXpixels/700, white, [0 0], 2);
        vbl  = Screen('Flip', windowPtr, vbl + (waitframes -0.5) * ifi);

        % LSL marker to check screen flip frequency
        % send markers into the outlet
%         mrk = markers{5};
%         outlet2.push_sample({mrk});   % note that the string is wrapped into a cell-array

        n=n+1;
    end
    % send markers into the outlet
    mrk = markers{2};
    outlet2.push_sample({mrk});   % note that the string is wrapped into a cell-array
    BaselineDuration=toc
    % flash a photocell to mark the end of baseline taking
    Screen('FillRect', windowPtr, white, RightUpperSquare);  % event type = 1200001
    Screen('FillRect', windowPtr, white, RightBottomSquare);
    Screen('Flip',windowPtr);
    Screen('Flip',windowPtr);

    
    instructionStart=['OK. Press a key to start!']; % Tell subject to open eye and start
    DrawFormattedText2(instructionStart,'win',windowPtr,...
        'sx','center','sy','center','xalign','center','yalign','center','baseColor',white);
    Screen('Flip',windowPtr);
    % hit a key to continue
    %pause(1); 
    KbStrokeWait;Screen('Flip',windowPtr);pause(5);% to separate baseline markers and trial markers
    % ************************************************ Baseline taking 60s 

    % initialize some variables
    TrialDurations=[]

    %############################### Loop through block
    %*******************************Loop through block
    for block=1:numBlock

        
        %******************** Within block loop
        % how many trials 
        for t=1:numTrials 

            % Breakout by hitting esc
            [keyIsDown, keysecs, keyCode] = KbCheck;  
            if keyCode(KbName('escape'))
            Screen('CloseAll');
            break;
            end

            % Show trial and block number at the bottom
%             Showtrial=['Beginning trial ' num2str(t) ' / ' num2str(numTrials) ', in block ' num2str(block) ' / ' num2str(numBlock)];
%             DrawFormattedText2(Showtrial,'win',windowPtr,...
%             'sx','center','sy', 'center','xalign','center','yalign','top','baseColor',white);
% 
%             % flip to screen and pause for 1 sec
%             vbl=Screen('Flip', windowPtr);
%             pause(2); % or WaitSecs(1);

            % pick a condition from randomized set allPerm
            conditionSelected = allPerm(numTrials*(block-1)+t);
            
            % Run the trial #################
            run trial_checkerboard_LSL_fNIR.m
            
%             if block ~= numBlock | t ~= numTrials % only bypass the last trial
                % Show Resting
%                 Resting = ['Rest for 20s'];
%                 DrawFormattedText2(Resting,'win',windowPtr,...
%                     'sx','center','sy', 'center','xalign','center','yalign','top','baseColor',white);
%                 vbl=Screen('Flip',windowPtr);
                
%             end

        end
        
    end
% initials trial frame number
% Show The End
TheEnd = ['The End'];
DrawFormattedText2(TheEnd,'win',windowPtr,...
    'sx','center','sy', 'center','xalign','center','yalign','top','baseColor',white);
vbl=Screen('Flip',windowPtr);
WaitSecs(3)
% hit a key to continue
% pause(1);
KbStrokeWait;Screen('Flip',windowPtr);pause(2);

%*************************************
Priority(0);   
sca;

catch
    sca;
    psychrethrow(psychlasterror);
end  
