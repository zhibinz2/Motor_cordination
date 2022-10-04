cd /home/hnl/Documents/GitHub/Motor_cordination/1_over_f/exp_uni_2p/Syncopation_2P_1.3Hz_230taps

%% This version is syncopated tapping for 2 person over an extended 3 monitors-window
% There are 4 conditions
% This is the full length version
sca; clc; close all; clear all; clearvars; 

%% load previous synchronization experiment results
last_seed=20221003;
% cd /home/hnl/Acquisition/zhibin/20220609_2P/Segmented_data/1_50Hz_ICAautomized/
% [BPCondi_dataName]=uigetfile('BP*.mat');% select the Left player EEG
% Path_filename=[pwd '/' BPCondi_dataName];
% load(Path_filename);
% cd /home/hnl/Acquisition/zhibin/20220713_2P/Segmented_data/
% load('data_variables20220713.mat'); save('BP20220713.mat','BP');
% load('/home/hnl/Acquisition/zhibin/20220713_2P/Segmented_data/BP20220713.mat');
% cd /home/hnl/Acquisition/zhibin/20220804_2P/Segmented_data/
% load('data_variables20220804.mat'); save('BP20220804.mat','BP');
% load('/home/hnl/Acquisition/zhibin/20220804_2P/Segmented_data/BP20220804.mat');
% save('BP20220810.mat','BP');
% load('/home/hnl/Acquisition/zhibin/20220810_2P/Segmented_data/BP20220810.mat');
% save('BP20220815.mat','BP');
% load('/home/hnl/Acquisition/zhibin/20220815_2P/Segmented_data/BP20220815.mat');
% load(['/home/hnl/Acquisition/zhibin/20221003_2P/Segmented_data/BP20221003.mat']);
load(['/home/hnl/Acquisition/zhibin/' num2str(last_seed) '_2P/Segmented_data/BP' num2str(last_seed) '.mat']);
% load data_structure20220713.mat; condiSeq=data(1).condiSeq; save('condiSeq20220713.mat','condiSeq');
% load 'condiSeq20220713.mat';
% load data_structure20220804.mat; save('condiSeq20220804.mat','condiSeq');
% load '/home/hnl/Acquisition/zhibin/20220804_2P/Segmented_data/condiSeq20220804.mat';
% save('condiSeq20220810.mat','condiSeq');
% load '/home/hnl/Acquisition/zhibin/20220810_2P/Segmented_data/condiSeq20220810.mat';
% save('condiSeq20220815.mat','condiSeq');
% load '/home/hnl/Acquisition/zhibin/20220815_2P/Segmented_data/condiSeq20220815.mat';
load(['/home/hnl/Acquisition/zhibin/20220815_2P/Segmented_data/condiSeq' num2str(last_seed) '.mat']);
%% Go back to experimental directory
cd /home/hnl/Documents/GitHub/Motor_cordination/1_over_f/exp_uni_2p/Syncopation_2P_1.3Hz_230taps

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
conditions = [1 2 3 4];
conditionNames={'uncoupled' 'L-lead' 'R-lead' 'mutual-1.3Hz'};  % conditionNames{1}
% conditionNames={'uncoupled-2Hz' 'L-lead-2Hz' 'R-lead-2Hz' 'mutual-2Hz' 'mutual-3Hz' 'mutual-3Hz-faster'};  % conditionNames{1}
% Instruction for the experimentor
ConditionInstructions={'uncoupled-1.3Hz' 'L-lead-1.3Hz' 'R-lead-1.3Hz' 'mutual-1.3Hz'}; 
% Instruction for the left player
ConditionInstructionsL={'Tap between red.'... %  \nPacer 2Hz.
    'Tap between red.'... %  \nPacer 2Hz
    'Tap between blue.'... % \nPacer 2Hz
    'Tap between blue.'}; 
% Instruction for the right player
ConditionInstructionsR={'Tap between blue.'... %  \nPacer 2Hz.
    'Tap between red.'... %  \nPacer 2Hz
    'Tap between blue.'... % \nPacer 2Hz
    'Tap between red.'}; 


% Block & Trial number of the experiment **********************************
% % number of taps per trial/condition
% numTaps=600; % 600 @@@@@@@@ 
% number of trials per block
numTrials=12;
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

% if not randomized
allPerm=[1:12];
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
    % Mean stimulus interval for 1.3Hz pacing
    MeanTapInterval13Hz=1/1.3; % about 0.75 second
    NumFramesInterval13Hz=round(MeanTapInterval13Hz/(ifi*waitframes));  
    % Resample previous tapping recording
    sr=2000; nTrials=size(BP(1).BP,2);
    for i=1:nTrials
        Showframes(1).Showframes{i} = resamBP(BP(1).BP{i},ifi,sr,waitframes);
        Showframes(2).Showframes{i} = resamBP(BP(2).BP{i},ifi,sr,waitframes);
    end
    % save(['Showframes' num2str(seed) '.mat'],'Showframes');

    % combine all 4 conditions
%     Showframes=zeros(4,MaxTaps);
%     Showframes(1,1:length(ShowframesL))= ShowframesL; 
%     Showframes(2,1:length(ShowframesR))= ShowframesR;
%     Showframes(3,1:length(ShowframesR))= ShowframesR;
%     Showframes(4,1:length(ShowframesR))= ShowframesR;
    % figure;plot(Showframes(1,:),'r.');hold on;plot(Showframes(2,:),'b.');

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
    run restingEEG.m %&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


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
            % Showframeselected=Showframes(conditionSelected,:);
            
            % select color
            if condiSeq(conditionSelected) == 1;
                color = green;
            elseif condiSeq(conditionSelected) == 2;
                color = red;
            elseif condiSeq(conditionSelected) == 3;
                color = blue;
            else condiSeq(conditionSelected) == 4;
                color = megenta;
            end

            % show instruction for each trial / condition
            % for experimenter
            ShowCondition=[ConditionInstructions{condiSeq(conditionSelected)}]; % Center monitor
            DrawFormattedText2([ShowCondition ' -- trial ' num2str(t)],'win',windowPtr,...
            'sx','center','sy', 'center','xalign','center','yalign','top','baseColor',color);
            % for left player
            if condiSeq(conditionSelected) == 1;
                ShowCondition=[ConditionInstructionsL{condiSeq(conditionSelected)}]; % Left monitor
                DrawFormattedText2(ShowCondition,'win',windowPtr,...
                'sx',xCenterL,'sy', 'center','xalign','center','yalign','top','baseColor',red);
            elseif condiSeq(conditionSelected) == 4;
                ShowCondition=[ConditionInstructionsL{condiSeq(conditionSelected)}]; % Left monitor
                DrawFormattedText2(ShowCondition,'win',windowPtr,...
                'sx',xCenterL,'sy', 'center','xalign','center','yalign','top','baseColor',blue);
            else
                ShowCondition=[ConditionInstructionsL{condiSeq(conditionSelected)}]; % Left monitor
                DrawFormattedText2(ShowCondition,'win',windowPtr,...
                'sx',xCenterL,'sy', 'center','xalign','center','yalign','top','baseColor',color);
            end
            % for right player
            if condiSeq(conditionSelected) == 1;
                ShowCondition=[ConditionInstructionsR{condiSeq(conditionSelected)}]; % Right monitor
                DrawFormattedText2(ShowCondition,'win',windowPtr,...
                'sx',xCenterR,'sy', 'center','xalign','center','yalign','top','baseColor',blue);
            elseif condiSeq(conditionSelected) == 4;
                ShowCondition=[ConditionInstructionsR{condiSeq(conditionSelected)}]; % Right monitor
                DrawFormattedText2(ShowCondition,'win',windowPtr,...
                'sx',xCenterR,'sy', 'center','xalign','center','yalign','top','baseColor',red);
            else
                ShowCondition=[ConditionInstructionsR{condiSeq(conditionSelected)}]; % Right monitor
                DrawFormattedText2(ShowCondition,'win',windowPtr,...
                'sx',xCenterR,'sy', 'center','xalign','center','yalign','top','baseColor',color);
            end
            
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
            
            % select the showframes for L and R players
            ShowframesL=Showframes(1).Showframes{t};
            ShowframesR=Showframes(2).Showframes{t};
            
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
                run IntervalResting.m %&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
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
TheEnd = ['The End. \nThank you!'];
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