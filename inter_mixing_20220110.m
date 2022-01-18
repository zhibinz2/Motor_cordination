% This edition change the 3/4 circle into 1/2 circle. And remove the connecting dots
sca;
clc;
close all;
clear all; 
clearvars; 

addpath Bimanual
addpath DrawSquares
addpath DrawCurves
addpath Conditions
addpath DrawReach
addpath PerformanceAnaylisis

AssertOpenGL;
 
if ~IsLinux
  error('Sorry, this demo currently only works on a Linux.');    
end

%**************************************************************************
% the set of conditions
% conditions = [1:1:11]*(pi/2/6);% 11 conditions from (0-180) degrees
% conditions = [0:1:12]*(pi/2/6);% 13 conditions from [0-180] degrees
% conditions = [0:1:10]*(pi/2/10);% block 1: 11 conditions from [0-90] degrees
% conditions = [1:1:10]*(pi/2/10)+pi/2;% block 2: 10 conditions from (90-180] degrees
% conditions = [0:1:10]*(pi/2/10)+pi/4;% 11 conditions from [45-135] degrees

% conditions in degree (7 condition with set movement ratio of 0:1, 1:4, 1:2, 1:1, 2:1, 4:1, 1:0)
% conditions = [0 2*atan(1/4)*180/pi 2*atan(1/2)*180/pi 2*atan(1)*180/pi 2*atan(2/1)*180/pi 2*atan(4/1)*180/pi 90];

% conditions in radian (7 condition with set movement ratio of 0:1, 1:4, 1:2, 1:1, 2:1, 4:1, 1:0)
conditions = [0 2*atan(1/4) 2*atan(1/2) 2*atan(1) 2*atan(2/1) 2*atan(4/1) pi];
conditionNames={'0:4' '1:4' '1:2' '1:1' 'log' 'exp' '4:0'};  % conditionNames{1}
% conditions = [2*atan(1/4) 2*atan(1/2) 2*atan(1) 2*atan(2/1) 2*atan(4/1)];

%************************************** Randomization of the experiment
% set the random number seed as the date of today in formate such as 20210809
seed=input('enter the date in format YYYYMMDD:');
behaviraldata.subjectnumber=seed;
ExptStartTime=char(datetime('now'));
behaviraldata.ExptStartTime=ExptStartTime;
rng(seed);

% *************************************************************************
% number of trials per block
numTrials=25;
% number of blocks
numBlock=length(conditions)*2;%length(conditions);
% total trial number
numtotal=numTrials*numBlock; 
% num of conditions in the experiment
numconditions=length(conditions);
% how many semirandom permutation set in the experiment 
numPerm=numtotal/numconditions;
% create the whole set of random conditions for the experiment
allPerm=[];
for p=1:numPerm
    allPerm=[allPerm randperm(numconditions)];
end



% *************************************************************************
% keep a record of the scores
TotalScore=0;% total score so far
ScoreLR=0;% initiate the score in each trial
GreenPixelsCovered=[]; %initiate the green pixels covered in each trial
TrialScores=[];% keep of a record of all trial scores

% set monetary reward
% fullBonusPerTrial=1;% $0.10 per trial if perfectly perfomed 
% fullBonus=fullBonusPerTrial*numtotal;
fullBonus=15;


% *************************************************************************



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

% DuoMice:****************************************************************
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
HideCursor(screenNumber,mice(2));
HideCursor(windowPtr,mice(1));
% *************************************************************************

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

% set the length of the reach ############################################
radius=screenYpixels/3;
% Set default connecting dot size
ConnectDotSize=120; 
% the thickness of the line
Thickness=ConnectDotSize/2;

% Create a fixation cross
FixCrX=[xCenter-Thickness/2:xCenter+Thickness/2 repmat(xCenter,1,Thickness+1)];
FixCrY=[repmat(yCenter,1,Thickness+1) yCenter-Thickness/2:yCenter+Thickness/2];

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

% Check if ifi=0.0167
if round(1/ifi)~=60
  error('Error: Screen flash frequency is not set at 60Hz.');    
end

% Set defaut joint dot color
Jcolor=green;

% Set size of the squares for photocell 
PhotosensorSize=30;
% Positions of the four corners
% RightBottomSquare= [screenXpixels-PhotosensorSize*2 screenYpixels-PhotosensorSize*2 screenXpixels screenYpixels];
% RightUpperSquare= [screenXpixels-PhotosensorSize*2 0 screenXpixels PhotosensorSize*2];
% LeftBottomSquare= [0 screenYpixels-PhotosensorSize*2 PhotosensorSize*2 screenYpixels];
% LeftUpperSquare= [0 0 PhotosensorSize*2 PhotosensorSize*2];
RightBottomSquare= [screenXpixels-PhotosensorSize*2 screenYpixels/2+230-PhotosensorSize screenXpixels screenYpixels/2+230+PhotosensorSize];
RightUpperSquare= [screenXpixels-PhotosensorSize*2 screenYpixels/2+110-PhotosensorSize screenXpixels screenYpixels/2+110+PhotosensorSize];
LeftBottomSquare= [0 screenYpixels/2+230-PhotosensorSize PhotosensorSize*2 screenYpixels/2+230+PhotosensorSize];
LeftUpperSquare= [0 screenYpixels/2+110-PhotosensorSize PhotosensorSize*2 screenYpixels/2+110+PhotosensorSize];

% *************************************************************************
% Time management here ####################################################
% #########################################################################
% #########################################################################
% #########################################################################
% #########################################################################
% Setting default mouse Position for some time
planSecs =0.5 ; % rest 1 s to look at trial number, visual evoked potential
numFramesPlan = round (planSecs/ifi);

% Length of time and number of frames we will use for each drawing trial
moveSecs = 1; %4; % 4 s to move
numFramesMove = round(moveSecs / ifi);

% total number of frames per trial
numFrames=numFramesPlan+numFramesMove; 

% *************************************************************************
% time duration to show bonus and trial number, and take rest before trial
Tintertrial=1;
numFramesRest = round (Tintertrial/ifi);

% time pause to smooth transition at the start of each trial
% Keep in mind: it takes about an additional 3 seconds to place the mice
Tpause=0.5;

% Estimate the time per trial
EstimateNavigateTime=3;
TimeTrial=planSecs+moveSecs+Tintertrial+Tpause+EstimateNavigateTime;
% Estimate total experiment time in minuts
TimeTotal=TimeTrial*numtotal/60;

% *************************************************************************
% time break between block
Tinterblock=3; 
% *************************************************************************

% Starting introduction
instructionStart=['You will be controlling two pens to do the task.'...
        '\n\n\n Each trial start when you have placed both pens at the center starting point.'...
        '\n\n\n Left pen will be shown in a red dot. Right shown in blue. Starting point in white.'...
        '\n\n\n Starting point will turn green, when you had placed both mice on it. Then the trial will start.'...
        '\n\n\n Movement of your left hand goes to the left, right hand goes to the right.'...
        '\n\n\n Move in horizontal direction to control the position of the dot.'...
        '\n\n\n You will not be able to see your pens in red and blue dots during the trials.'...
        '\n\n\n But the dot will turn green when comes close to a trajectory line shown on the screen.'...
        '\n\n\n Try to coordinate the moving speed of your hands in suggested ratio'...
        '\n\n\n Move the dot along the trajectory to the end of the line. Complete the movement within ' num2str(moveSecs) ' seconds.'...
        '\n\n\n You will achieve higer score if performed more quickly and closely to the trajectory.'...
        '\n\n\n Questions? If none, press any key to continue.']
DrawFormattedText2(instructionStart,'win',windowPtr,...
    'sx','center','sy','center','xalign','center','yalign','center','baseColor',white);
Screen('Flip',windowPtr);
% WaitSecs(Tinterblock);
% hit a key to continue
KbStrokeWait;

%*************************************************************************
% No resting EEG taking in session 2 and 3 ###############################
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................
% ........................................................................

%*******************************Loop through block
for block=1:numBlock
    
    % Initialize some values
    n = 1;
    
    % To keep a record of the percentage of time inside the square
    behaviraldata.dataBlock(block).dataTrialNumInside=[];
    %xLyL=[]; xRyR=[]; % to keep track of mouse trace
    behaviraldata.dataBlock(block).dataTrialL.xLyL=[];
    behaviraldata.dataBlock(block).dataTrialR.xRyR=[];
    behaviraldata.dataBlock(block).dataTrialJ.xJyJ=[];
    %save block number info
    behaviraldata.dataBlock(block).blockNumber=block;
    

    %************************************ show bonus before block and rest
    % Show performents
    if block ~= 1
        % ScoreLR=0.2;TotalScore=0.2 % just to test alignment
        Showbonus = ['Score: ' sprintf('%0.2f %%', ScoreLR*100) '  Average: ' sprintf('%0.2f %%', TotalScore*100)];
        DrawFormattedText2(Showbonus,'win',windowPtr,...
            'sx','center','sy', yCenter+screenYpixels/20,'xalign','center','yalign','top','baseColor',white);
    end
    
    % Show bonus criteria
%     ShowbonusCrit = ['Average  Bonus\n    >65%  $3   \n    >70%  $6   \n    >75%  $9   \n'...
%         '    >80%  $12   \n    >85%  $15   \n'];
    ShowbonusCrit = ['Average  Grade\n    <60%  F    \n    >60%  D    \n    >65%  C    \n    >70%  B-   \n    >75%  B    \n'...
        '    >80%  B+   \n    >85%  A    \n    >90%  A+   \n'];
    DrawFormattedText2(ShowbonusCrit,'win',windowPtr,...
            'sx','center','sy', yCenter+screenYpixels/10,'xalign','center','yalign','top','baseColor',white);
    
    % Show progress    
    Showblock = ['Hit a key to begin block ' num2str(block) ' / ' num2str(numBlock)]
    DrawFormattedText2(Showblock,'win',windowPtr,...
        'sx','center','sy', yCenter+screenYpixels/3,'xalign','center','yalign','top','baseColor',white);
    Screen('Flip',windowPtr);

%     WaitSecs(Tinterblock); % to separate blocks in the photocell signal
    % hit a key to continue
    KbStrokeWait;

    %*******************************************************

    % get a timestamp at the start of block
      vbl = Screen('Flip', windowPtr);
    
    %#################### Loop throought trials Within block 
    %******************** Within block loop
    % how many trials 
    for t=1:numTrials 
        % Breakout by hitting esc
        [keyIsDown, keysecs, keyCode] = KbCheck;  
        if keyCode(KbName('escape'))
        Screen('CloseAll');
        break;
        end

         %************ Show bonus of previous trial
        if t ~= 1
            for Restframes=1:numFramesRest
                Showbonus=['Score: ' sprintf('%0.2f %%', ScoreLR*100) ';      Average: ' sprintf('%0.2f %%', TotalScore*100)];
                DrawFormattedText2(Showbonus,'win',windowPtr,...
                'sx','center','sy', yCenter+screenYpixels/20,'xalign','center','yalign','top','baseColor',white); 
            
%                 Showtrial=['Beginning trial ' num2str(t) ' / ' num2str(numTrials) ', in block ' num2str(block) ' / ' num2str(numBlock)];
%                 DrawFormattedText2(Showtrial,'win',windowPtr,...
%                 'sx','center','sy', yCenter+screenYpixels/5,'xalign','center','yalign','top','baseColor',white);
                %Flip to the screen   
                vbl  = Screen('Flip', windowPtr, vbl + (waitframes -0.5) * ifi);
            end
        end            
        %*************************Randomized selection
        % pick a condition from randomized set allPerm
        conditionSelected = allPerm(numTrials*(block-1)+t);
        % produced the position parameters
            rad_ang=conditions(conditionSelected);
            [x,y] = drawReach(radius,rad_ang, xCenter, yCenter);
        % produced the polyshape area for score calculation
            % polyTrajatory = polyshape([x(1)-3*Thickness/2 x(1)+3*Thickness/2 x(end)+3*Thickness/2 x(end)-3*Thickness/2],...
            %   [y(1)+3*Thickness/2 y(1)-3*Thickness/2 y(end)-3*Thickness/2 y(end)+3*Thickness/2]);
            
            % the circle version
            polyTrajatory=polycircle(x(1),y(1),Thickness/2);
            for cir=2:length(x)
                polyTrajatory=union(polyTrajatory,polycircle(x(cir),y(cir),Thickness/2));
            end
        
        % initialize a small green pixel polyshap area for score calculation
            %PolyshapGreenPixelInitial=polyshape([x(1)-Thickness/2 x(1)+Thickness/2 x(1)+Thickness/2 x(1)-Thickness/2],...
            %   [y(1)-Thickness/2 y(1)-Thickness/2 y(1)+Thickness/2 y(1)+Thickness/2]);
            
             % the circle version
            PolyshapGreenPixelInitial=polycircle(x(1),y(1),Thickness/2); 
        
            %figure;plot(polyTrajatory);hold on; plot(PolyshapGreenPixelInitial); hold off;
            %ylim([0 screenYpixels]);xlim([0 screenXpixels]);set(gca, 'YDir', 'reverse');
            PolyshapeUnion=intersect(polyTrajatory,PolyshapGreenPixelInitial);
            %figure;plot(PolyshapeUnion);ylim([0 screenYpixels]);xlim([0 screenXpixels]);set(gca, 'YDir', 'reverse');area(PolyshapeUnion)
            
        %***********************************************
        
        % Trial start only when both mice placed at the starting point
        Insidestart = 0;
        while Insidestart == 0
            % If esc is press, break out of the while loop and close the screen
            [keyIsDown, keysecs, keyCode] = KbCheck;
            if keyCode(KbName('escape'))
                Screen('CloseAll');
                break;
            end
            
            % Draw the route
            % Screen('DrawDots', windowPtr, [x;y], Thickness, white, [0 0], 2);

            % Display the starting point 
            Screen('DrawDots', windowPtr, [x(1) y(1)], Thickness, white, [], 2);
            
            % Display instruction
            textPlan=['Place both pens at the starting points'];
            DrawFormattedText2(textPlan,'win',windowPtr,...
                'sx','center','sy', yCenter+screenYpixels/20,'xalign','center','yalign','top','baseColor',white);
            
            % Show trial and block number
            Showtrial=['Beginning trial ' num2str(t) ' / ' num2str(numTrials) ', in block ' num2str(block) ' / ' num2str(numBlock)];
            DrawFormattedText2(Showtrial,'win',windowPtr,...
            'sx','center','sy', yCenter+screenYpixels/3,'xalign','center','yalign','top','baseColor',white);
            
            % Get mouse location
            [xML0, yML0] = GetMouse(windowPtr,mice(2));
            [xMR0, yMR0] = GetMouse(windowPtr,mice(1));

            % Shift the mouse location to map the tablets
%             xML=xML0/2; % the upper left quatrand of the scrren
%             yML=yML0/2; % the upper left quatrand of the scrren
%             xMR=xCenter+xMR0/2; % the upper right quatrand of the screen
%             yMR=yMR0/2;% the upper right quatrand of the screen
            % xML=xML0/2; % the left side of the screen
            xML=xML0*(11/20); % a bit more than half of the left side of the scrren
            yML=yML0; % the left side of the scrren
            % xMR=xCenter+xMR0*(11/20); % the right side of the screen
            xMR=(xCenter-screenXpixels/20)+xMR0*(11/20); % a bit more than half of the right side of the scrren
            yMR=yMR0;% the right side of the screen
     
            % Display the cursor as a dot
            Screen('DrawDots', windowPtr, [xML yML], Thickness, red, [], 2);
            Screen('DrawDots', windowPtr, [xMR yMR], Thickness, blue, [], 2);
            
            % initial the before values
            xMLbefore=xML;
            yMLbefore=yML;
            xMRbefore=xMR;
            yMRbefore=yMR;
            
            xJbefore=[]; yJbefore=[];% initialize the joint before value
            
            % flip to screen
            vbl  = Screen('Flip', windowPtr, vbl + (waitframes -0.5) * ifi);
            
            % update the while loop
            Insidestart = (sqrt((xML-xCenter).^2+(yML-yCenter).^2) < Thickness/2/5) & ...
                (sqrt((xMR-xCenter).^2+(yMR-yCenter).^2)<Thickness/2/5);
            
            if Insidestart == 1
                Screen('DrawDots', windowPtr, [x(1) y(1)], Thickness, green, [], 2);
                % flip to screen
                vbl  = Screen('Flip', windowPtr, vbl + (waitframes -0.5) * ifi);
                % pause to smooth transition
                WaitSecs(Tpause);
            end
            
        end
        
        % reset the trial score to zero
        ScoreLR=0;
        GreenPixelsCovered=[];

        % get a timestamp at the start of the trial
        vbl = Screen('Flip', windowPtr);
        % Run one trial 
        run inter_run_trial_20220110.m 


        % reset n
        n=1;

        % DuoMice:
        % Show master cursors again:
%         for mouse = mice
%         ShowCursor('Arrow', [], mouse);
%         end


        % To keep a record of the percentage of time inside the square
        behaviraldata.dataBlock(block).dataTrialNumInside(t).NumInside=NumInside;
        
        % Store behaviral data
        behaviraldata.dataBlock(block).dataTrialL(t).xLyL=xLyL;
        behaviraldata.dataBlock(block).dataTrialR(t).xRyR=xRyR;
        behaviraldata.dataBlock(block).dataTrialJ(t).xJyJ=xJyJ;
        
        %save trial condition
        behaviraldata.dataBlock(block).dataTrialNumInside(t).condition=conditionSelected;
        behaviraldata.dataBlock(block).dataTrialL(t).condition=conditionSelected;
        behaviraldata.dataBlock(block).dataTrialR(t).condition=conditionSelected;
        behaviraldata.dataBlock(block).dataTrialJ(t).condition=conditionSelected;
        
        % calculate the score in this trial
        % Count interset of polyshape areas
%         GreenPixelsCovered(1,:)=[]; % remove the initial []
%         GreenPixelsCoveredxJ=GreenPixelsCovered(:,1);GreenPixelsCoveredxJ=GreenPixelsCoveredxJ';
%         GreenPixelsCoveredyJ=GreenPixelsCovered(:,1);GreenPixelsCoveredyJ=GreenPixelsCoveredyJ';
%         polyGreenPixelsCovered=polyshape(GreenPixelsCoveredxJ,GreenPixelsCoveredyJ);
        
    
   
%         figure;
%         plot(polyGreenPixelsCovered);
%         hold on;
%         plot(polyTrajatory);
%         hold off;
            
        % update the scores
        % ScoreLR=fullBonusPerTrial*area(PolyshapeUnion)/area(polyTrajatory);
        ScoreLR=1;%area(PolyshapeUnion)/area(polyTrajatory); % use percentile reward
        % Save trial scores
        TrialScores=[TrialScores ScoreLR];
        
        % TotalScore=TotalScore+ScoreLR;
        if TotalScore==0;
            TotalScore=ScoreLR;
        else
            TotalScore=mean(TrialScores); % show average score instead
        end
        
        % reset trial score values
        clear PolyshapGreenPixelInitial; 
        clear polyTrajatory; 
        clear PolyshapeUnion; 

    end 

end


% Show The End ************************************************************

% TotalReward=fullBonus*TotalScore;

% Record the end of experiment
ExptEndTime=char(datetime('now'));
behaviraldata.ExptEndTime=ExptEndTime;

% Show performance
Showbonus = ['Score: ' sprintf('%0.2f %%', ScoreLR*100) '  Average: ' sprintf('%0.2f %%', TotalScore*100)];
DrawFormattedText2(Showbonus,'win',windowPtr,...
            'sx','center','sy', yCenter+screenYpixels/20,'xalign','center','yalign','top','baseColor',white);
        
% Show bonus criteria
% ShowbonusCrit = ['Average  Bonus\n    >65%  $3   \n    >70%  $6   \n    >75%  $9   \n'...
%         '    >80%  $12   \n    >85%  $15   \n'];
ShowbonusCrit = ['Average  Grade\n    <60%  F    \n    >60%  D    \n    >65%  C    \n    >70%  B-   \n    >75%  B    \n'...
        '    >80%  B+   \n    >85%  A    \n    >90%  A+   \n'];
DrawFormattedText2(ShowbonusCrit,'win',windowPtr,...
            'sx','center','sy', yCenter+screenYpixels/10,'xalign','center','yalign','top','baseColor',white);
    
% Show the end
TheEnd = ['The End'];
DrawFormattedText2(TheEnd,'win',windowPtr,...
    'sx','center','sy', yCenter+screenYpixels/3,'xalign','center','yalign','top','baseColor',white);
Screen('Flip',windowPtr);
WaitSecs(3)
% hit a key to continue
KbStrokeWait;


%*************************************


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

% Save the behavioral data
run SaveBehaviraldata.m
cd ..
run git_control.m