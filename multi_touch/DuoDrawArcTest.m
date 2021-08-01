 
sca;
close all;
clearvars;


% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);
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
% luminace values are (in general) defined between 0 and 1.
% For help see`: help WhiteIndex and help BlackIndex
white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);

% Open an on screen window and color it black
% For help see: Screen OpenWindow?
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, black);

% Get the size of the on screen window in pixels
% For help see: Screen WindowSize?
[screenXpixels, screenYpixels] = Screen('WindowSize', window);

% Get the centre coordinate of the window in pixels
% For help see: help RectCenter
[xCenter, yCenter] = RectCenter(windowRect);

% Enable alpha blending for anti-aliasing
% For help see: Screen BlendFunction?
% Also see: Chapter 6 of the OpenGL programming guide
Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


% We can define a center for the dot coordinates to be relaitive to. Here
% we set the centre to be the centre of the screen
dotCenter = [xCenter yCenter];

N      = 400; % Define the number of dots, more will draw a clean arc
Radius = screenXpixels /4 ; % Define the radius of the arc, we use a radius that fit in the screen
[x1,y1]  = DrawArc(Radius,N); % Create the arc coordinates
Size   = 15 .* ones(1,N);
Color  = ones(3,N);
x1      = flip(x1 + Radius); % Add an offset for center the arc
y1      = abs(y1-max(y1)); % Reverse the arc
y1      = y1 + (yCenter - max(y1) + (max(y1) - min(y1))/2); % Offset the middle of the arc to be center on the middle of the screen
PosL    = [x1; y1];


[x2,y2]  = DrawArc(Radius,N); % Create second arc coordinates
Size   = 15 .* ones(1,N);
Color  = ones(3,N);
x2      = flip(x2 + Radius*3); % Add an offset for center the arc
y2      = abs(y2-max(y2)); % Reverse the arc
y2      = y2 + (yCenter - max(y2) + (max(y2) - min(y2))/2); % Offset the middle of the arc to be center on the middle of the screen
PosR    = [x2; y2];


% Make a base Rect of 50 by 50 pixels
baseRect = [0 0 50 50];

% Set the color of the rect to red
rectColor1 = [1 0 0];
rectColor2 = [1 0 0];

time = 0;

% Sync us and get a time stamp
vbl = Screen('Flip', window);
waitframes = 1;

% Maximum priority level
topPriorityLevel = MaxPriority(window);
Priority(topPriorityLevel);
ifi = Screen('GetFlipInterval', window);

SetMouse(PosL(1,1),PosL(2,1), window);

n = 1;
red = [1 0 0];
blue = [0 0 1];
% Loop the animation until a key is pressed
HideCursor

% To keep a record of the percentage of time inside the square
NumInside1=[] 
NumInside2=[]
   
while n<length(PosL)%~KbCheck
%     for mouse = mice
    [xM, yM, buttons] = GetMouse(window);
%     [xM1, yM1, buttons1] = GetMouse(window, mice(1));
%     [xM2, yM2, buttons2] = GetMouse(window, mice(2));
    
    centeredRect1 = CenterRectOnPointd(baseRect, PosL(1,n), PosL(2,n));
    centeredRect2 = CenterRectOnPointd(baseRect, PosR(1,n), PosR(2,n));
    
    % Draw the rect to the screen
    Screen('FillRect', window, rectColor1, centeredRect1);
    Screen('FillRect', window, rectColor2, centeredRect2);
    Screen('DrawDots', window, PosL,Size, Color, [0 0], 2);
    Screen('DrawDots', window, PosR,Size, Color, [0 0], 2);
    
    % See if the mouse cursor is inside the square
    inside1 = IsInRect(xM, yM, centeredRect1);
    inside2 = IsInRect(xM, yM, centeredRect2);
    
    % Set the color of the square to be red if the mouse is inside it or
    % blue if the mouse is outside it
    if inside1 == 1
        rectColor1 = red;
    elseif inside1 == 0
        rectColor1 = blue;
    end
    
    if inside2 == 1
        rectColor2 = red;
    elseif inside2 == 0
        rectColor2 = blue;
    end
    
   
    Screen('DrawDots', window, [xM yM], 16, [0.5 0.5 0.5], [], 2);
    
    
    % Flip to the screen
    vbl  = Screen('Flip', window, vbl + waitframes * ifi);
    
    NumInside1(n)=inside1;
    NumInside2(n)=inside2;
     
    xy(n,1) = xM;
    xy(n,2) = yM;
    n = n+1;
    
end

% the percentage of time in the square 
PercentIn1=sum(NumInside1)/length(PosL); 
PercentIn2=sum(NumInside2)/length(PosR); 

Screen('Flip', window);

% Now we have drawn to the screen we wait for a keyboard button press (any
% key) to terminate the demo. For help see: help KbStrokeWait
KbStrokeWait;

% Clear the screen. "sca" is short hand for "Screen CloseAll". This clears
% all features related to PTB. Note: we leave the variables in the
% workspace so you can have a look at them if you want.
% For help see: help sca
sca;