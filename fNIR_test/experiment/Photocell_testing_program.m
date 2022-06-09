.
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
numBlock=9; % 3 blocks per session,6 subject
% total trial number
numtotal=numTrials*numBlock; 
% **********************************Block & Trial number of the experiment


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

  
    % Set size of the squares for photocell ###############################
    PhotosensorSize=30; % 30;
    % Positions of the four corners
    RightBottomSquare= [screenXpixels-PhotosensorSize*2 screenYpixels-PhotosensorSize*2 screenXpixels screenYpixels];
    RightUpperSquare= [screenXpixels-PhotosensorSize*2 0 screenXpixels PhotosensorSize*2];
    LeftBottomSquare= [0 screenYpixels-PhotosensorSize*2 PhotosensorSize*2 screenYpixels];
    LeftUpperSquare= [0 0 PhotosensorSize*2 PhotosensorSize*2];

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
    % total number of frames per trial
    numFrames=round(10/ifi/waitframes); % 10 seconds 
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

    instructionStart=['OK. Press a key to start!']; % Tell subject to open eye and start
    DrawFormattedText2(instructionStart,'win',windowPtr,...
        'sx','center','sy','center','xalign','center','yalign','center','baseColor',white);
    Screen('Flip',windowPtr);
    % hit a key to continue
    %pause(1); 
    KbStrokeWait;Screen('Flip',windowPtr);pause(2);
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
            
            % LSL markers to the local network
            % send data into the outlet, sample by sample
        %     outlet1.push_sample(1); % send data number 1
            % send markers into the outlet
            mrk = markers{3};
            outlet2.push_sample({mrk});   % note that the string is wrapped into a cell-array
            
            % Run the trial #################
            % initial or reset trial frame number
            n=1;
            % get a timestamp and begin the trial
            vbl = Screen('Flip', windowPtr);
                        
            
            % while numberOfSecondsElapsed < 10;
            while n <  numFrames 
                % If esc is press, break out of the while loop and close the screen
                [keyIsDown, keysecs, keyCode] = KbCheck;
                if keyCode(KbName('escape'));
                    Screen('CloseAll');
                    break;
                end

                % Flash photocell every other second during the whole trial (bottom left)
                if (~isempty(find([1:2:numFrames]==n))) % every two frames
                    Screen('FillRect', windowPtr, white, RightUpperSquare);  % event type = 1200001
                    Screen('FillRect', windowPtr, white, RightBottomSquare);
                end
        
                % Flip to the screen
                vbl  = Screen('Flip', windowPtr, vbl + (waitframes -0.5) * ifi);
            
                n = n+1;
            end

            %% send markers into the outlet
            mrk = markers{4};
            outlet2.push_sample({mrk});   % note that the string is wrapped into a cell-array

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
