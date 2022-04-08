%% Matlab function
I = checkerboard;
imshow(I)

n=20
I = checkerboard(n);
imshow(I)

n=50;p=2;q=3;
I = checkerboard(n,p,q);
imshow(I)

K = (checkerboard(90,5,5) > 0.5);
figure, imshow(K)

greenSquare=repmat(cat(3,0,1,0),10,10);
whiteSquare=repmat(cat(3,1,1,1),10,10);
QuadraSquare1=[greenSquare whiteSquare;whiteSquare greenSquare];
QuadraSquare2=[whiteSquare greenSquare;greenSquare whiteSquare];
K1=repmat(QuadraSquare1,5,5);
K2=repmat(QuadraSquare2,5,5);
figure, imshow(K1)
figure, imagesc(K1)

C = [0 2 4 6; 8 10 12 14; 16 18 20 22];
imagesc(C)
imshow([0 0.5])

K=repmat([green white;white green],5,5);
%% Psychtoolbox 
Screen MakeTexture?
% textureIndex=Screen('MakeTexture', WindowIndex, imageMatrix [, optimizeForDrawAngle=0] [, specialFlags=0] [, floatprecision] [, textureOrientation=0] [, textureShader=0]);
Screen DrawTexture?
% Screen('DrawTexture', windowPointer, texturePointer [,sourceRect] [,destinationRect] [,rotationAngle] [, filterMode] [, globalAlpha] [, modulateColor] [, textureShader] [, specialFlags] [, auxParameters]);

%% Tryout https://www.mathworks.com/matlabcentral/answers/349487-psychtoolbox-checkerboard-alternation-error
% Clear the workspace and the screen
sca;
close all;
clearvars;
%Screen('Preference', 'SkipSyncTests', 1);
% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);
% Get the screen numbers
screens = Screen('Screens');
% Draw to the external screen if avaliable
screenNumber = max(screens);
% Define black and white
white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);
grey = white / 2;
green = [0 1 0];
% Open an on screen window
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, grey, [0 0 400 400]);
ifi = Screen('GetFlipInterval', window);
% Query the maximum priority level
topPriorityLevel = MaxPriority(window);
% Get the size of the on screen window
[screenXpixels, screenYpixels] = Screen('WindowSize', window);
%%%%%%%%%Stripe information
% Get the centre coordinate of the window
[xCenter, yCenter] = RectCenter(windowRect);    
cx = (screenXpixels/2);
cy = (screenYpixels/2);
% Make a base Rect of 200 by 200 pixels
dimx = cx;
dimy = cy;
baseRect = [0 0 dimx dimy];
pos = [cx- dimx/2 ,cy - dimy/2,cx+ dimx/2 ,cy + dimy/2];
% MakeTexture: for help, see Screen DrawTexture?
Screen('ColorRange',window);
textureIndex=Screen('MakeTexture',window,K1);
Screen('DrawTexture', window, textureIndex, [],pos);
vbl = Screen('Flip', window);
vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
%% If screen crash, 
% To get out, press "Ctrl+0" to navigate to the invisible command window, 
% then "Ctrl + C" and "sca"
