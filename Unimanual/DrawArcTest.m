    
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
Radius = screenYpixels /1.3 ; % Define the radius of the arc, we use a radius that fit in the screen
[x,y]  = DrawArc(Radius,N); % Create the arc coordinates
Size   = 15 .* ones(1,N);
Color  = ones(3,N);
x      = flip(x + Radius); % Add an offset for center the arc
y      = abs(y-max(y)); % Reverse the arc
y      = y + (yCenter - max(y) + (max(y) - min(y))/2); % Offset the middle of the arc to be center on the middle of the screen
Pos    = [x; y];


% Make a base Rect of 50 by 50 pixels
baseRect = [0 0 50 50];

% Set the color of the rect to red
rectColor = [1 0 0];


time = 0;

% Sync us and get a time stamp
vbl = Screen('Flip', window);
waitframes = 1;

% Maximum priority level
topPriorityLevel = MaxPriority(window);
Priority(topPriorityLevel);
ifi = Screen('GetFlipInterval', window);

SetMouse(Pos(1,1),Pos(2,1), window);

n = 1;
red = [1 0 0];
blue = [0 0 1];
% Loop the animation until a key is pressed
HideCursor

NumInside=[] % To keep a record of the percentage of time inside the square

while n<length(Pos)%~KbCheck
    
    [xM, yM, buttons] = GetMouse(window);
    
    centeredRect = CenterRectOnPointd(baseRect, Pos(1,n), Pos(2,n));
    
    if n == 1
        % Draw a white dot on the right-bottom corner of the screen for the photosensor
        RightBottom = [screenXpixels-15 screenYpixels-15]; % position of the right-bottom to draw the dot
        Screen('DrawDots', window, RightBottom, 15, ones(3,1), [0 0], 2);
    elseif n == length(Pos)-1
        % Draw a white dot on the right-bottom corner of the screen for the photosensor
        RightBottom = [screenXpixels-15 screenYpixels-15]; % position of the right-bottom to draw the dot
        Screen('DrawDots', window, RightBottom, 15, ones(3,1), [0 0], 2);
    end  
    
    
    % Draw the rect to the screen
    Screen('FillRect', window, rectColor, centeredRect);
    Screen('DrawDots', window, Pos,Size, Color, [0 0], 2);
    
    % See if the mouse cursor is inside the square
    inside = IsInRect(xM, yM, centeredRect);
    
    % Set the color of the square to be red if the mouse is inside it or
    % blue if the mouse is outside it
    if inside == 1
        rectColor = red;
    elseif inside == 0
        rectColor = blue;
    end
    
   
    Screen('DrawDots', window, [xM yM], 16, [0.5 0.5 0.5], [], 2);
    
    
    % Flip to the screen
    vbl  = Screen('Flip', window, vbl + waitframes * ifi);
    
    NumInside(n)=inside;
    
    xy(n,1) = xM;
    xy(n,2) = yM;
    n = n+1;
    
end

% the percentage of time in the square 
PercentIn=sum(NumInside)/length(Pos); 


Screen('Flip', window);

% Now we have drawn to the screen we wait for a keyboard button press (any
% key) to terminate the demo. For help see: help KbStrokeWait
KbStrokeWait;

% Clear the screen. "sca" is short hand for "Screen CloseAll". This clears
% all features related to PTB. Note: we leave the variables in the
% workspace so you can have a look at them if you want.
% For help see: help sca
sca;