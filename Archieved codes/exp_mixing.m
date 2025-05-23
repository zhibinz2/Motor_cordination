 sca;
close all;
clc;
clear;
clear all; 
clearvars; 

addpath Bimanual
addpath DrawSquares
addpath DrawCurves

AssertOpenGL;
 
if ~IsLinux
  error('Sorry, this demo currently only works on a Linux.');    
end

%************************************************
% the set of conditions
conditionfuctions = {@A1_A2, @A1_A1, @A1_U2, @A1_U1, ...
@A2_A1, @A2_A2, @A2_U1, @A2_U2, ...
@U1_U2, @U1_U1, @U1_A2, @U1_A1, ...
@U2_U1, @U2_U2, @U2_A1, @U2_A2 };
    
%************************************** Randomization of the experiment
% set the random number seed as the date of today in formate such as 20210809
seed=input('enter the date in format YYYYMMDD:');
data.subjectnumber=seed;
rng(seed);

% number of trials
numTrials=50;
% number of blocks
numBlock=8;
% total trial number
numtotal=numTrials*numBlock; % should be 400
% num of conditions in the experiment
numconditions=16;
% how many semirandom permutation set in the experiment 
numPerm=numtotal/numconditions;
% create the whole set of random conditions for the experiment
allPerm=[];
for p=1:numPerm
    allPerm=[allPerm randperm(16)];
end
%*******************************************


try      
    
% Here we call some default settings for setting up Psychtoolbox 
PsychDefaultSetup(2);
% Start with black screen
Screen('Preference', 'VisualDebugLevel', 1); 
Screen('Preference', 'SkipSyncTests', 1);

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
% For help see`: help WhiteIndex and help BlackIndex
white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);

% Initialize some colors
red   = [1 0 0];
blue  = [0 0 1];
green = [0 1 0];
white = [white white white];
black = [black black black];
grey  = [0.5 0.5 0.5];

% Open an on screen window and color it black
% For help see: Screen Openwindow?
[windowPtr, windowRect] = PsychImaging('Openwindow', screenNumber, black);
% Get the size of the on screen windowPtr in pixels
% For help see: Screen windowSize?
[screenXpixels, screenYpixels] = Screen('windowSize', screenNumber);
% Get the centre coordinate of the window in pixels
% For help see: help RectCenter
[xCenter, yCenter] = RectCenter(windowRect); 
% Enable alpha blending for anti-aliasing
% For help see: Screen BlendFunction?
% Also see: Chapter 6 of the OpenGL programming guide
Screen('BlendFunction', windowPtr, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

% We can define a center for the dot coordinates to be relaitive to. Here
% we set the centre to be the centre of the screen
dotCenter = [xCenter yCenter];

% set the side length of square
steplength=screenXpixels/6;

% Numer of frames to wait when specifying good timing. Note: the use of
% wait frames is to show a generalisable coding. For example, by using
% waitframes = 2 one would flip on every other frame. See the PTB
% documentation for details. In what follows we flip every frame.  
waitframes = 1;

% Retreive the maximum priority number
topPriorityLevel = MaxPriority(windowPtr); 
% set Priority once at the start of a script after setting up onscreen window.
Priority(topPriorityLevel);

% Measure the vertical refresh rate of the monitor
ifi = Screen('GetFlipInterval', windowPtr);

% Set default connecting dot size
ConnectDotSize=60; 
% Set default connecting dot color to white
[ConnectDotColorL0, ConnectDotColorL1, ConnectDotColorL2, ConnectDotColorL3] = deal(white);
[ConnectDotColorR0, ConnectDotColorR1, ConnectDotColorR2, ConnectDotColorR3] = deal(white);

% Set size of the squares for photocell 
PhotosensorSize=30;
% Positions of the photocell
RightBottomSquare= [screenXpixels-PhotosensorSize*2 screenYpixels/2+230-PhotosensorSize screenXpixels screenYpixels/2+230+PhotosensorSize];
RightUpperSquare= [screenXpixels-PhotosensorSize*2 screenYpixels/2+110-PhotosensorSize screenXpixels screenYpixels/2+110+PhotosensorSize];
LeftBottomSquare= [0 screenYpixels/2+230-PhotosensorSize PhotosensorSize*2 screenYpixels/2+230+PhotosensorSize];
LeftUpperSquare= [0 screenYpixels/2+110-PhotosensorSize PhotosensorSize*2 screenYpixels/2+110+PhotosensorSize];


% DuoMice:
% Get handles for all virtual pointing devices, aka cursors:
typeOnly='masterPointer'; 
mice = GetMouseIndices(typeOnly);  

% DuoMice:
% Hide the system-generated cursors. We do this, because only the
% first mouse cursor is hardware-accelerated, ie., a GPU created
% hardware cursor. All other cursors are software-cursors, created
% by the Windowing system. These tend to flicker badly in our use
% case. Therefore we disable all system cursor images and draw our
% cursors ourselves for a more beautiful look:
HideCursor;

for block=1:numBlock
    
    % Initialize some values
    n = 1;
    % NumInside=[]; % To keep a record of the percentage of time inside the square
    %xLyL=[]; xRyR=[]; % to keep track of mouse trace
    data.dataBlock.dataTrialL.xLyL=[];
    data.dataBlock.dataTrialR.xRyR=[];
    %save block number info
    data.dataBlock(block).blockNumber=block;
    

    %************************************ show test and rest
    %Loads a window and waits for input to start recording - just for getting set up and ready
%     instructionStart = 'Press the space bar to begin when you are ready. It will start automatically in 5 minutes.';
%     Screen('DrawText', windowPtr, instructionStart, screenXpixels/3, screenYpixels/6, white); 
%     Screen(windowPtr, 'Flip');  
    
    instructionStart = ['Rest or hit a key to begin block ' num2str(numBlock)];
    DrawFormattedText2(instructionStart,'win',windowPtr,...
        'sx','center','sy','center','xalign','center','yalign','center','baseColor',white);
    Screen('Flip',windowPtr);
    % hit a key to continue
    KbStrokeWait;


%     KbName('UnifyKeyNames');
%     spaceKeyID = KbName('space'); 
%     %Waits for space bar       
%     [keyIsDown, secs, keyCode] = KbCheck;
%     while keyCode(spaceKeyID)~=1
%         [keyIsDown, secs, keyCode] = KbCheck;
%     end
    %*******************************************************

    % get a timestamp at the start of block
    vbl = Screen('Flip', windowPtr);
    
    

    %******************** Within block loop
    % how many trials 
    for t=1:numTrials 
        % Breakout by hitting esc
        [keyIsDown, keysecs, keyCode] = KbCheck;  
        if keyCode(KbName('escape'))
        Screen('CloseAll');
        break;
        end

        %************ Show trial number and rest
        restSecs =0.5; % rest 1 s to rest and look at trial number
        numFramesRest = round (restSecs/ifi);
        for Restframes=1:numFramesRest 
        %Screen('DrawText', windowPtr, ['trial ' num2str(t)], screenXpixels/3, screenYpixels/6, white); 
        DrawFormattedText2(['trial ' num2str(t) ' / ' num2str(numTrials)],'win',windowPtr,...
        'sx','center','sy','center','xalign','center','yalign','center','baseColor',white);
        % Flip to the screen   
        vbl  = Screen('Flip', windowPtr, vbl + (waitframes -0.5) * ifi);
        end      
        
        %*************************Randomized selection
        % pick a condition from randomized set allPerm
        conditionSelected = allPerm(numTrials*(block-1)+t);
        % produced the position parameters
            [xL,yL,PosL,PosL0,PosL1,PosL2,PosL3,ThicknessL,ColorL,...
                xR,yR,PosR,PosR0,PosR1,PosR2,PosR3,ThicknessR,ColorR] ...
                = conditionfuctions{conditionSelected}(steplength,yCenter,screenXpixels);
        %***********************************************


        % Setting default mouse Position for some time
        planSecs =1 ; % rest 1 s to look at trial number
        numFramesPlan = round (planSecs/ifi);

        % Length of time and number of frames we will use for each drawing trial
        moveSecs = 4; % 4 s to move
        numFramesMove = round(moveSecs / ifi);

        % total number of frames
        numFrames=numFramesPlan+numFramesMove; 

        % Run one trial 
        run run_trial.m 

        % reset n
        n=1;

        % DuoMice:
        % Show master cursors again:
        for mouse = mice
        ShowCursor('Arrow', [], mouse);
        end

        % Store xLyL
        %dataL(t).xLyL=xLyL; 
        %dataR(t).xRyR=xRyR;
        data.dataBlock(block).dataTrialL(t).xLyL=xLyL;
        data.dataBlock(block).dataTrialR(t).xRyR=xRyR;
        
        %save trial condition
        data.dataBlock(block).dataTrialL(t).condition=conditionSelected;
        data.dataBlock(block).dataTrialR(t).condition=conditionSelected;

    end
    


end



Priority(0);   
   
% Now we have drawn to the screen we wait for a keyboard button press (any
% key) to terminate the demo. For help see: help KbStrokeWait
% KbStrokeWait;                                 
 
% Clear the screen. "sca" is short hand for "Screen CloseAll". This clears
% all features related to PTB. Note: we leave the variables in the
% workspace so you can have a look at them if you want.
% For help see: help sca
sca;

catch
  sca
  psychrethrow(psychlasterror);
end  

