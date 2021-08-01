 sca;
close all;
clear all;
clearvars;
% Here we call some default settings for setting up Psychtoolbox 
PsychDefaultSetup(2);
Screen('Preference', 'SkipSyncTests', 1);
% Get the screen numbers. This gives us a number for each of the screens
% attached to our computer. For help see: Screen Screens?
screens = Screen('Screens');
% sca;
% Draw we select the maximum of these numbers. So in a situation where we
% have two screens attached to our monitor we will draw to the external
% screen. When only one screen is attached to the monitor we will draw to
% this. For help see: help max
screenNumber = max(screens);
% Define black and white (white will be 1 and black 0). This is because 
% luminace values are (in general) defined between 0 and 1.
% For help see`: help WhiteIndex and help BlackIndex
white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);
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
%  sca;
% We can define a center for the dot coordinates to be relaitive to. Here
% we set the centre to be the centre of the screen
dotCenter = [xCenter yCenter];
steplength=400;

[xL,yL]  = DrawSquare(steplength); % Create  the route
[xR,yR]  = DrawSquare(steplength); % Create  the route
[PosL,PosL0,PosL1,PosL2,PosL3,ThicknessL,ColorL] = DrawPosL(xL,yL,steplength,yCenter);
[PosR,PosR0,PosR1,PosR2,PosR3,ThicknessR,ColorR] = DrawPosR(xR,yR,steplength,yCenter,screenXpixels);

% Make a base Rect of 50 by 50 pixels
% baseRect = [0 0 50 50];

% Set the default color of the rect cue to red
% rectColor = [1 0 0];

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

% Initialize some values
n = 1;
red   = [1 0 0];
blue  = [0 0 1];
green = [0 1 0];
white = [white white white];
black = [black black black];
grey  = [0.5 0.5 0.5];

% Set default connecting dot size
ConnectDotSize=60; 
% Set default connecting dot color to blue
% ConnectDotColorR0 = blue;
% ConnectDotColorR1 = blue;
% ConnectDotColorR2 = blue;
% ConnectDotColorR3 = blue;
[ConnectDotColorL0, ConnectDotColorL1, ConnectDotColorL2, ConnectDotColorL3] = deal(white);
[ConnectDotColorR0, ConnectDotColorR1, ConnectDotColorR2, ConnectDotColorR3] = deal(white);


% Set size of the squares for photosenors 
PhotosensorSize=30;
% Positions of the four corners
RightBottomSquare= [screenXpixels-PhotosensorSize*2 screenYpixels-PhotosensorSize*2 screenXpixels screenYpixels];
RightUpperSquare= [screenXpixels-PhotosensorSize*2 0 screenXpixels PhotosensorSize*2];
LeftBottomSquare= [0 screenYpixels-PhotosensorSize*2 PhotosensorSize*2 screenYpixels];
LeftUpperSquare= [0 0 PhotosensorSize*2 PhotosensorSize*2];

% Loop the animation until a key is pressed
HideCursor

% DuoMice:
% Get handles for all virtual pointing devices, aka cursors:
typeOnly='masterPointer'; 
mice = GetMouseIndices(typeOnly); 
 
% NumInside=[]; % To keep a record of the percentage of time inside the square
xLyL=[]; xRyR=[]; % to keep track of mouse trace
dataR=[];dataL=[];

%*******************************************************************************************
%Loads a window and waits for input to start recording - just for getting set up and ready
instructions = 'Press the space bar to begin';
Screen('DrawText', windowPtr, instructions, screenXpixels/3, screenYpixels/2, white); 
Screen(windowPtr, 'Flip');


KbName('UnifyKeyNames');
spaceKeyID = KbName('space'); 
%Waits for space bar       
[keyIsDown, secs, keyCode] = KbCheck;
while keyCode(spaceKeyID)~=1
    [keyIsDown, secs, keyCode] = KbCheck;
end
%*******************************************************************************************
% get a timestamp at the start of trials
vbl = Screen('Flip', windowPtr);
% Flash once to mart the start of the trial planing
Screen('FillRect', windowPtr, white, RightBottomSquare);
Screen('FillRect', windowPtr, white, RightUpperSquare);
Screen('FillRect', windowPtr, white, LeftBottomSquare);
Screen('FillRect', windowPtr, white, LeftUpperSquare);
% Flash to mark the start of the trial planning   
vbl  = Screen('Flip', windowPtr, vbl + (waitframes -0.5) * ifi);   
%*******************************************************************************************
% Length of time and number of frames we will use for each drawing trial
numSecs = 3;
numFrames = round(numSecs / ifi);

% how many trials

for t=1:3
    
    % Breakout
    [keyIsDown, keysecs, keyCode] = KbCheck;
    if keyCode(KbName('escape'))
    Screen('CloseAll');
    break;
    end
    
%     % Show trial number
%     for frames=1:30 % planing phase
%     Screen('DrawText', windowPtr, ['trial ' num2str(t)], screenXpixels/3, screenYpixels/2, white); 
%     % Flip to the screen   
%     vbl  = Screen('Flip', windowPtr, vbl + (waitframes -0.5) * ifi);
%     end      
%     %Waits for space bar      
%     [keyIsDown, secs, keyCode] = KbCheck;
%     while keyCode(spaceKeyID)~=1 
%         [keyIsDown, secs, keyCode] = KbCheck;
%     end 

    % Setting default mouse Position
    SetMouse(PosL(1,1),PosL(2,1), window, mice(2));
    SetMouse(PosR(1,1),PosR(2,1), window, mice(1));

    % Run one trial 
    run BiOneTrial.m 
    % reset n
    n=1;
    
    % Store xLyL
    dataL(t).xLyL=xLyL; 
    dataR(t).xRyR=xRyR;
      
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

