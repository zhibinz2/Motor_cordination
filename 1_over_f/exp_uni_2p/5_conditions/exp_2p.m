%% This version is 1 over f tapping for 2 person over an extended 3 monitors-window
% There are 5 conditions
% This is the full length version
cd /home/hnl/Documents/GitHub/Motor_cordination/1_over_f/exp_uni_2p/5_conditions
sca; clc; close all; clear all; clearvars; 

%% set keyboards
mydir  = pwd; % or cd
idcs   = strfind(mydir,'/');
newdir = mydir(1:idcs(end-2)-1);
cd (newdir);
addpath Cedrus_Keyboard/
run setCedrusRB.m % plug in the left RB pad first then the right RB pad
run SuppressWarning.m 
cd (mydir);

%% Check OpenGL
% Break and issue an error message if the installed Psychtoolbox is not
% based on OpenGL or Screen() is not working properly.
AssertOpenGL;

% if ~IsLinux
%   error('Sorry, this demo currently only works on a Linux.');    
% end

%% Set trial conditions ****************************************************
conditions = [1 2 3 4 5 6];
conditionNames={'uncoupled' 'L-lead' 'R-lead' 'mutual-2Hz' 'mutual-3Hz' 'mutual-faster'};  % conditionNames{1}
% conditionNames={'uncoupled-2Hz' 'L-lead-2Hz' 'R-lead-2Hz' 'mutual-2Hz' 'mutual-3Hz' 'mutual-3Hz-faster'};  % conditionNames{1}
% Instruction for the experimentor
ConditionInstructions={'uncoupled-2Hz' 'L-lead-2Hz' 'R-lead-2Hz' 'mutual-2Hz' 'mutual-3Hz' 'mutual-3Hz-faster'}; 
% Instruction for the left player
ConditionInstructionsL={'Tap on your own.\nPacer start at 2Hz.'... %  \nPacer 2Hz.
    'Tap on your own.\nPacer start at 2Hz.'... %  \nPacer 2Hz
    'Syn with player R in blue.'... % \nPacer 2Hz
    'Pacer start at 2Hz.\nThen syn with each other.'... % 
    'Pacer start at 3Hz.\nThen syn with each other.'... % \nPacer 2Hz
    'Pacer start at 3Hz.\nThen syn with each other as fast as you can.'}; 
% Instruction for the right player
ConditionInstructionsR={'Tap on your own.\nPacer start at 2Hz.'... %  \nPacer 2Hz.
    'Syn with player L in red.'... %  \nPacer 2Hz
    'Tap on your own.\nPacer start at 2Hz.'... % \nPacer 2Hz
    'Pacer start at 2Hz.\nThen syn with each other.'... % 
    'Pacer start at 3Hz.\nThen syn with each other.'... % \nPacer 2Hz
    'Pacer start at 3Hz.\nThen syn with each other as fast as you can.'}; 

% Block & Trial number of the experiment **********************************
% number of taps per trial/condition
numTaps=600; % 600 @@@@@@@@ 
% number of trials per block
numTrials=6;
% number of blocks
numBlock=1;
% total trial number
numtotal=numTrials*numBlock; 
% num of conditions in the experiment
numconditions=length(conditions);
% **********************************Block & Trial number of the experiment

% Randomization of the experiment #########################################
% set the random number seed as the date of today in formate such as 20210809
seed=input('enter the date in format YYYYMMDD:');
rng(seed);

% % how many semirandom permutation set in the experiment 
% numPerm=numtotal/numconditions;
% % create the whole set of random conditions for the experiment
% allPerm=[];
% for p=1:numPerm
%     allPerm=[allPerm randperm(numconditions)];
% end

% if not randomized
allPerm=[1 2 3 4 5 6];
% ######################################### Randomization of the experiment 


%% ########################################################################
try      
 
    %% Here we call some default settings for setting up Psychtoolbox 
    PsychDefaultSetup(2);

    % Start with black screen
    % Removes the blue screen flash and minimize extraneous warnings.
    Screen('Preference', 'VisualDebugLevel', 0); % disable all visual alerts
    Screen('Preference', 'SkipSyncTests', 1); %  shorten the maximum duration of the sync tests to 3 seconds worst case
    Screen('Preference', 'SuppressAllWarnings', 1); %disable all output to the command window 
    % The use of “mirror mode”, or “clone mode”, where multiple displays show the
    % same content, will almost always cause  timing and performance problems.

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
    megenta = [1 0 1];% fill([0 1 1 0],[0 0 1 1],megenta)
    cyan = [0 1 1]; % fill([0 1 1 0],[0 0 1 1],cc)
    purple = [0.6 0.1 0.9];

    % Open an on screen window and color it black
    % For help see: Screen Openwindow?
    % This will draw on a black backgroud with a size of [0 0 500 1000] and
    % return a window pointer windowPtr
%     [windowPtr, windowRect] = PsychImaging('Openwindow', screenNumber, black, [0 0 590*3*60/59 330*40/33]); 
    [windowPtr, windowRect] = PsychImaging('Openwindow', screenNumber, black); 

    % Get the size of the on screen windowPtr in pixels
    % For help see: Screen windowSize?
    [screenXpixels, screenYpixels] = Screen('windowSize', windowPtr);

    % Get the centre coordinate of the window in pixels
    % For help see: help RectCenter
    [xCenter, yCenter] = RectCenter(windowRect); 

    % Get centre coordinate of the left and right monitors in the display  window in pixels
    xCenterL=screenXpixels/6; xCenterR=5*screenXpixels/6;
    
    % Stimulus design ***********************************************
    % Just a cross
    
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
    
    % set 60-144 Hz Screen rate????????????????????????????????????????????
    % Measure the vertical refresh rate of the monitor ????????????????????
    ifi = Screen('GetFlipInterval', windowPtr);

    % Check if ifi=0.0167
%     if round(1/ifi)~=60
%       error('Error: Screen flash frequency is not set at 60 Hz.');    
%     end

    % Check if ifi=0.0167
    if round(1/ifi)~=100
      error('Error: Screen flash frequency is not set at 100 Hz.');    
    end

    % Check if ifi=0.0083
%     if round(1/ifi)~=120
%       error('Error: Screen flash frequency is not set at 144 Hz.');    
%     end

    % Check if ifi=0.0069
%     if round(1/ifi)~=144
%       error('Error: Screen flash frequency is not set at 144Hz.');    
%     end
    % ?????????????????????????????????????????????????????????????????????

    %% Set size of the squares for photocell ###############################
    PhotosensorSize=30;
    % Positions of the four corners
    RightBottomSquare= [screenXpixels-PhotosensorSize*2 screenYpixels-PhotosensorSize*2 screenXpixels screenYpixels];
    RightUpperSquare= [screenXpixels-PhotosensorSize*2 0 screenXpixels PhotosensorSize*2];
    LeftBottomSquare= [0 screenYpixels-PhotosensorSize*2 PhotosensorSize*2 screenYpixels];
    LeftUpperSquare= [0 0 PhotosensorSize*2 PhotosensorSize*2];

    %% Numer of frames to wait when specifying good timing. 
    % Note: the use of wait frames is to show a generalisable coding. 
    % For example, by using waitframes = 2 one would flip on every other frame. See the PTB
    % documentation for details. In what follows we flip every frame. 
    % In order to flip with feedback at 50 Hz
    % We have to filp every (1/50)/ifi frames or 1/50 ms/frame
    waitframes=2;

    %% Randomization of the pace stimuli in each condition ****************************************
    % Mean stimulus interval for 2Hz pacing
    MeanTapInterval2Hz=0.5; % second
    NumFramesInterval2Hz=round(MeanTapInterval2Hz/(ifi*waitframes));  
    
    % condition 1-4 (paced the frist 30 taps, set the rest of the frames with value zeros)
    Showframes1=[1:NumFramesInterval2Hz:NumFramesInterval2Hz*30 zeros(1,570)]; 
    Showframes2=Showframes1;
    Showframes3=Showframes1;
    Showframes4=Showframes1;
    
    % Mean stimulus interval for 3Hz pacing
    MeanTapInterval3Hz=1/3; % second
    NumFramesInterval3Hz=round(MeanTapInterval3Hz/(ifi*waitframes));  % on average 72 frames per stimulus 
    % condition 5
    Showframes5=[1:NumFramesInterval3Hz:NumFramesInterval3Hz*30 zeros(1,570)]; 
    Showframes6=Showframes5;
    
    % combine all 6 conditions
    Showframes=[Showframes1;Showframes2;Showframes3;Showframes4;Showframes5;Showframes6];

    %% Setting time variables**********************************************
    % Length of one minute baseline
    BaselineLength = 2; % in seconds
    numFramesBaseline = round(BaselineLength / ifi / waitframes);
    
    %% Hide Mice:****************************************************************
    % Get handles for all virtual pointing devices, aka cursors:
    typeOnly='masterPointer'; 
    mice = GetMouseIndices(typeOnly);  
    HideCursor(windowPtr,mice);

    %% Create a fixation cross
    FixCrX=[xCenter-round(screenYpixels/200):xCenter+round(screenYpixels/200)-1 repmat(xCenter,1,round(screenYpixels/100)+1)];
    FixCrY=[repmat(yCenter,1,round(screenYpixels/100)+1) yCenter-round(screenYpixels/200)+1:yCenter+round(screenYpixels/200)];
   
    %% Take resting EEG
    run restingEEG.m

    %% ############################### Loop through block
    % *******************************Loop through block
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

            % pick a condition from randomized set allPerm
            conditionSelected = allPerm(numTrials*(block-1)+t);
            Showframeselected=Showframes(conditionSelected,:);
            
            % select color
            if conditionSelected == 1
                color = green;
            elseif conditionSelected == 2
                color = red;
            elseif conditionSelected == 3
                color = blue;
            elseif conditionSelected == 4
                color = megenta;
            else (conditionSelected == 5) | (conditionSelected == 6) 
                color = purple;
            end

            % show instruction for each trial / condition
            ShowCondition=[ConditionInstructions{conditionSelected}]; % Center monitor
            DrawFormattedText2(ShowCondition,'win',windowPtr,...
            'sx','center','sy', 'center','xalign','center','yalign','top','baseColor',color);
            ShowCondition=[ConditionInstructionsL{conditionSelected}]; % Left monitor
            DrawFormattedText2(ShowCondition,'win',windowPtr,...
            'sx',xCenterL,'sy', 'center','xalign','center','yalign','top','baseColor',color);
            ShowCondition=[ConditionInstructionsR{conditionSelected}]; % Right monitor
            DrawFormattedText2(ShowCondition,'win',windowPtr,...
            'sx',xCenterR,'sy', 'center','xalign','center','yalign','top','baseColor',color);
            
            Screen('TextSize',windowPtr, 10); % smaller size font
            Showtrial=['Read to tap?'...
                 '\nExperimenter hit a key to continue.']; 
            DrawFormattedText2(Showtrial,'win',windowPtr,...
            'sx','center','sy', screenYpixels*0.8,'xalign','center','yalign','top','baseColor',white);
            DrawFormattedText2(Showtrial,'win',windowPtr,...
            'sx',xCenterL,'sy', screenYpixels*0.8,'xalign','center','yalign','top','baseColor',white);
            DrawFormattedText2(Showtrial,'win',windowPtr,...
            'sx',xCenterR,'sy', screenYpixels*0.8,'xalign','center','yalign','top','baseColor',white);
            Screen('TextSize',windowPtr, 18); % reset default font size
            Screen('Flip', windowPtr);
            pause(1); 
            % hit a key to continue
            KbStrokeWait;
            Screen('Flip',windowPtr);% Reponse with a black screen
            pause(1); % given time to return the position of arm after key press (no need for 2p, experimentor control)


            % flip to screen and pause for 1 sec
%             vbl=Screen('Flip', windowPtr);
%             pause(2); % or WaitSecs(1);

            % LSL markers to the local network
            % send data into the outlet, sample by sample
        %     outlet1.push_sample(1); % send data number 1
            % send markers into the outlet
%             mrk = markers{3};
%             outlet2.push_sample({mrk});   % note that the string is wrapped into a cell-array

            
            % get a timestamp and begin the trial
            vbl = Screen('Flip', windowPtr);
            % pause 1 second so that the first condition don't come up
            % too suddenly
            % pause(1)

            % Run the trial #################
            run trial_2p.m

            % send markers into the outlet
%             mrk = markers{4};
%             outlet2.push_sample({mrk});   % note that the string is wrapped into a cell-array

            if block ~= numBlock | t ~= numTrials % only bypass the last trial
                run IntervalResting.m
%                 % Show Resting
%                 between_block_rest= 3; % in seconds
%                 Resting = ['Take a rest for at least ' num2str(between_block_rest) ' s. \n Experimentor hit a key to continue.'];
%                 DrawFormattedText2(Resting,'win',windowPtr,...
%                     'sx', 'center','sy', 'center','xalign','center','yalign','top','baseColor',white);
%                 DrawFormattedText2(Resting,'win',windowPtr,...
%                     'sx', xCenterL,'sy', 'center','xalign','center','yalign','top','baseColor',white);
%                 DrawFormattedText2(Resting,'win',windowPtr,...
%                     'sx', xCenterR,'sy', 'center','xalign','center','yalign','top','baseColor',white);
%                 vbl=Screen('Flip',windowPtr); 
%                 % Rest 1 sec
%                 pause(between_block_rest);
%                 % hit a key to continue
%                 KbStrokeWait;   
            end

        end
        
    end
% initials trial frame number
% Show The End
TheEnd = ['The End. /nThank you!'];
DrawFormattedText2(TheEnd,'win',windowPtr,...
    'sx','center','sy', 'center','xalign','center','yalign','top','baseColor',white);
DrawFormattedText2(TheEnd,'win',windowPtr,...
    'sx',xCenterL,'sy', 'center','xalign','center','yalign','top','baseColor',white);
DrawFormattedText2(TheEnd,'win',windowPtr,...
    'sx',xCenterR,'sy', 'center','xalign','center','yalign','top','baseColor',white);
vbl=Screen('Flip',windowPtr);
WaitSecs(3)
% hit a key to continue
KbStrokeWait;     

%*************************************
Priority(0);   
sca;

catch
    sca;
    psychrethrow(psychlasterror);
end  


%% save data
% cd /home/hnl/Documents/GitHub/1overf/stimulus_data_storage
cd /home/hnl/Acquisition/zhibin
filename=[num2str(seed) '.mat'];
save(filename);
% cd /home/hnl/Documents/GitHub/Motor_cordination/1_over_f/exp_uni_2p/5_conditions