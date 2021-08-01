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
% sca; 

% We can define a center for the dot coordinates to be relaitive to. Here
% we set the centre to be the centre of the screen
dotCenter = [xCenter yCenter];

steplength=400; 
[x,y]  = DrawSquare(steplength); % Create the route
Thickness   = 15 .* ones(1,length(x)); % thickness of the route
Color  = ones(3,length(x));
% Add an offset for center the arc
x      = flip(x + steplength); 
% Reverse the route
% y      = abs(y-max(y)); 
% Offset the middle of the arc to be center on the middle of the screen
y      = y + (yCenter - max(y) + (max(y) - min(y))/2); 
% the route
Pos    = [x; y];
% postion of connecting points
Pos0   = Pos(:,1);
Pos1   = Pos(:,steplength);
Pos2   = Pos(:,steplength*2-1);
Pos3   = Pos(:,steplength*3-2);

% Make a base Rect of 50 by 50 pixels
baseRect = [0 0 50 50];

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
ConnectDotColor0 = white;
ConnectDotColor1 = white;
ConnectDotColor2 = white;
ConnectDotColor3 = white; 

% Set size of the squares for photosenors 
PhotosensorSize=30;
% Positions of the four corners
RightBottomSquare= [screenXpixels-PhotosensorSize*2 screenYpixels-PhotosensorSize*2 screenXpixels screenYpixels];
RightUpperSquare= [screenXpixels-PhotosensorSize*2 0 screenXpixels PhotosensorSize*2];
LeftBottomSquare= [0 screenYpixels-PhotosensorSize*2 PhotosensorSize*2 screenYpixels];
LeftUpperSquare= [0 0 PhotosensorSize*2 PhotosensorSize*2];

% Loop the animation until a key is pressed
HideCursor

% NumInside=[]; % To keep a record of the percentage of time inside the square
xy=[]; % to keep track of mouse trace
data=[];

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
    
    
    % Show trial number
    for frames=1:30
    Screen('DrawText', windowPtr, ['trial ' num2str(t)], screenXpixels/3, screenYpixels/2, white); 
    % Flip to the screen   
    vbl  = Screen('Flip', windowPtr, vbl + (waitframes -0.5) * ifi);
    end      
%     %Waits for space bar      
%     [keyIsDown, secs, keyCode] = KbCheck;
%     while keyCode(spaceKeyID)~=1 
%         [keyIsDown, secs, keyCode] = KbCheck;
%     end 

    % Setting default mouse position
    SetMouse(Pos(1,1),Pos(2,1), screenNumber);
    
    % Run one trial 
    run OneTrial.m
      
    % Store xy
    dataUni(t).xy=xy; 
     
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

