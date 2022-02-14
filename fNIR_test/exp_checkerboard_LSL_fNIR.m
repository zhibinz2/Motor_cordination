% This version shows checkerboard on two hemi-fields and sends LSL markers
% to the fNIRS laptop
sca; clc; close all; clear all; clearvars; 

% LSL #####################################################################
% Addpath
addpath(genpath('/home/hnl/Documents/GitHub/labstreaminglayer/LSL/liblsl-Matlab')); % for LSL
% instantiate the library for LSL
lib = lsl_loadlib();
% make new stream outlets
info1 = lsl_streaminfo(lib,'BioSemi','EEG',8,100,'cf_float32','sdfwerr32432');
% the name (here MyMarkerStream) is visible to the experimenter and should be chosen so that 
% it is clearly recognizable as your MATLAB software's marker stream
% The content-type should be Markers by convention, and the next three arguments indicate the 
% data format (1 channel, irregular rate, string-formatted).
% The so-called source id is an optional string that allows for uniquely identifying your 
% marker stream across re-starts (or crashes) of your script (i.e., after a crash of your script 
% other programs could continue to record from the stream with only a minor interruption).
info2 = lsl_streaminfo(lib,'MyMarkerStream','Markers',1,0,'cf_string','myuniquesourceid23443');
% open the outlets
outlet1 = lsl_outlet(info1);
outlet2 = lsl_outlet(info2);
% Set markers
markers = {'Test-Marker', 'Trial-Start', 'Trial-End'};

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
numBlock=3;
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
    [windowPtr, windowRect] = PsychImaging('Openwindow', screenNumber, black, [0 0 600 400]); 

    % Get the size of the on screen windowPtr in pixels
    % For help see: Screen windowSize?
    [screenXpixels, screenYpixels] = Screen('windowSize', windowPtr);

    % Get the centre coordinate of the window in pixels
    % For help see: help RectCenter
    [xCenter, yCenter] = RectCenter(windowRect); 

    % Checkerboard Design****************************************
    % Black and white
    K1 = (checkerboard > 0.5);
    K2 = ~K1;
    % Green and white

    % MakeTexture
    textureIndex1=Screen('MakeTexture',windowPtr,double(K1));
    textureIndex2=Screen('MakeTexture',windowPtr,double(K2));
    % Position of the checkerboard
    posL = [xCenter-screenXpixels/12*5,yCenter-screenYpixels/8*3,xCenter-screenXpixels/12*1,yCenter+screenYpixels/8*3];
    posR = [xCenter+screenXpixels/12*1,yCenter-screenYpixels/8*3,xCenter+screenXpixels/12*5,yCenter+screenYpixels/8*3];
    % ****************************************Checkerboard Design
    
    % Enable alpha blending for anti-aliasing
    % For help see: Screen BlendFunction?
    % Also see: Chapter 6 of the OpenGL programming guide
    Screen('BlendFunction', windowPtr, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    % Set text display options for operating system other than Linux.
    % For help see: Screen TextFont?
    if ~IsLinux
        Screen('TextFont', windowPtr, 'Arial');
        Screen('TextSize', windowPtr, 18);
    end

    % Retreive the maximum priority number
    topPriorityLevel = MaxPriority(windowPtr); 
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
    
    %% Setting time variables**********************************************
    % total number of frames per trial
    numFrames=round(1/ifi)*10; % 10 seconds 
    %% **********************************************Setting time variables

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

            % Show trial and block number
            Showtrial=['Beginning trial ' num2str(t) ' / ' num2str(numTrials) ', in block ' num2str(block) ' / ' num2str(numBlock)];
            DrawFormattedText2(Showtrial,'win',windowPtr,...
            'sx','center','sy', 'center','xalign','center','yalign','top','baseColor',white);

            % flip to screen and pause for 1 sec
            Screen('Flip', windowPtr);
            pause(1); % or WaitSecs(1);

            % pick a condition from randomized set allPerm
            conditionSelected = allPerm(numTrials*(block-1)+t);

            % initial or reset trial frame number
            n=1;
            % get a timestamp and begin the trial
            vbl = Screen('Flip', windowPtr);

            run trial_checkerboard_LSL_fNIR.m
            
            if block ~= numBlock & t ~=numTrials
                % Show Resting
                Resting = ['Rest for 20s'];
                DrawFormattedText2(Resting,'win',windowPtr,...
                    'sx','center','sy', 'center','xalign','center','yalign','top','baseColor',white);
                Screen('Flip',windowPtr);
                % Rest 20 sec
                pause(20)
            end

        end
        
        
        % Show The End
        TheEnd = ['The End'];
        DrawFormattedText2(TheEnd,'win',windowPtr,...
            'sx','center','sy', 'center','xalign','center','yalign','top','baseColor',white);
        Screen('Flip',windowPtr);
        WaitSecs(3)
        % hit a key to continue
        KbStrokeWait;

        %*************************************
        Priority(0);   
        sca;

    end
% initials trial frame number

catch
    sca
    psychrethrow(psychlasterror);
end  
