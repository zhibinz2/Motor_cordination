%% go to saving directory
cd /ssd/zhibin/1overf/20220515_2P/Segmented_data/1_50Hz_ICAautomized
cd /ssd/zhibin/1overf/20220517_2P/Segmented_data/1_50Hz_ICAautomized

cd /ssd/zhibin/1overf/20220713_2P/Segmented_data
cd /ssd/zhibin/1overf/20220721_2P/Segmented_data
cd /ssd/zhibin/1overf/20220804_2P/Segmented_data/
cd /ssd/zhibin/1overf/20220808_2P/Segmented_data/
cd /ssd/zhibin/1overf/20220810_2P/Segmented_data/
cd /ssd/zhibin/1overf/20220811_2P/Segmented_data/
cd /ssd/zhibin/1overf/20220815_2P/Segmented_data/
cd /ssd/zhibin/1overf/20220816_2P/Segmented_data/
cd /ssd/zhibin/1overf/20221003_2P/Segmented_data/
%% Collect all the variables
% experimental data extracted =====================================
allPerm;
condiSeq=allPerm; % for synchronization
conditionNames;seed;ifi;

% select the session type
sessionTypes={'synch','synco'};
session=1; % synch
session=2; % synco

% add to mat file ***
filename_variables=['data_variables' num2str(seed) '.mat'];
% save(filename_variables,'allPerm'); % for synchronization
save(filename_variables,'condiSeq'); % for syncopation
save(filename_variables,'conditionNames','ifi','seed','-append');
save(filename_variables,'sessionTypes','session','-append');

% organize into data structure *** (skip)
data=[]; % first dimension has to be 1 in order to explore field names in GUI
% data.condiSeq=allPerm;
% data.conditionNames=conditionNames; 
% data.ifi=ifi; 
% data(1,1).condiSeq=allPerm; data(2,1).condiSeq=allPerm;
% data(1,2).conditionNames=conditionNames; data(2,2).conditionNames=conditionNames;
for i=1:2 % first dimension has to be 1 in order to show field names when explore in GUI
    % data(1,i).condiSeq=allPerm; % same as data(i).condiSeq=allPerm;
    data(1,i).condiSeq=condiSeq; % for syncopation
    data(1,i).seed=seed;
    data(1,i).session=session;
    data(1,i).condiNames=conditionNames;
    data(1,i).refreshRate=1/ifi;
end

% Recording variables extracted =====================================
timeL;samplesL;TRIGGERindL;srL;channels_infoL;labelsL;
timeR;samplesR;TRIGGERindR;srR;channels_infoR;labelsR;
if srL==srR; sr=srL; end;

% add to mat file ***
save(filename_variables,'sr','timeL','TRIGGERindL','srL','channels_infoL','labelsL',...
    'timeR','TRIGGERindR','srR','channels_infoR','labelsR','-append');

% organize into data structure ***(skip)
% data=[];
data(1,1).sr=srL; data(1,2).sr=srR;
% data(2,1).samples=samplesL; data(2,2).samples=samplesR;
% data=[];
% data(1,1).time=timeL; data(2,1).time=timeR;
data(1,1).time=timeL; data(1,2).time=timeR;
data=rmfield(data,'time');% try remove the field
% data(1,2).samples=samplesL; data(2,2).samples=samplesR;

% EEG variables extracted ==========================================
labels;% for plotting
chaninfo; % channel coordinates
mixedsigL;mixedsigR;

% add to mat file ***
save(filename_variables,'labels','chaninfo','-append');

% organize into data structure ***(skip)
for i=1:2
    data(1,i).labels=labels; 
    data(1,i).chaninfo=chaninfo; 
end


% Behaviral variables extracted =====================================
% BP from trigger chan (Response Pad)
BottonPresTimeIndL;BottonPresTimeL01;numPresL;
BottonPresTimeIndR;BottonPresTimeR01;numPresR;
% FB from trigger chan (light sensor 2)
FeedbTimeIndL;FeedbTimeL01;numFbL;
FeedbTimeIndR;FeedbTimeR01;numFbR;
% Pacers and Makers (light sensor 1)
PacerTimeIndL;PacerTimeL01;numPacerL;
PacerTimeIndR;PacerTimeR01;numPacerR;
% FB from light detector (ISO aux chan)
locsL; % Index in datatime; should be same as FeedbTimeIndL
locsR; % Index in datatime; should be same as FeedbTimeIndR
% add to mat file ***
% save(filename_variables,'','-append');
% organize into data structure ***

% EMG variables needed =============================================
filtered_EMGL;filtered_EMGR;
loL;loR; % Hilbert envelop
% add to mat file ***
% save(filename_variables,'','-append');
% organize into data structure ***

%% Segment EEG - synchronization
% mixedsigL=mixedsigL';
% mixedsigR=mixedsigR';
SegIndL;
SegIndR;

% organize into data structure *** 
% EEGtemp1=[];EEGtemp2=[];
for i=1:12
    data(1).EEG{i}=mixedsigL(SegIndL(2*i-1):SegIndL(2*i),:); % in cells
    data(2).EEG{i}=mixedsigR(SegIndR(2*i-1):SegIndR(2*i),:);
%     EEGtemp1(i).EEG=mixedsigL(SegIndL(2*i-1):SegIndL(2*i),:); % in sub-structure
%     EEGtemp2(i).EEG=mixedsigR(SegIndR(2*i-1):SegIndR(2*i),:);
end
% data(1,1).EEG=EEGtemp1;
% data(1,2).EEG=EEGtemp2;

% add to mat file ***
EEG=[];
EEG(1,1).EEG=data(1).EEG;
EEG(2,1).EEG=data(2).EEG;
save(filename_variables,'EEG','-append');

%% Segment EEG of the pacing period
% mixedsigL=mixedsigL';
% mixedsigR=mixedsigR';
Seg30IndL;
Seg30IndR;

% organize into data structure ***
for i=1:12
    data(1).EEG_pacing{i}=mixedsigL(Seg30IndL(2*i-1):Seg30IndL(2*i),:); % in cells
    data(2).EEG_pacing{i}=mixedsigR(Seg30IndR(2*i-1):Seg30IndR(2*i),:);
end

% add to mat file ***
EEG_pacing=[];
EEG_pacing(1,1).EEG_pacing=data(1).EEG_pacing;
EEG_pacing(2,1).EEG_pacing=data(2).EEG_pacing;
save(filename_variables,'EEG_pacing','-append');

%% segment button presses - synchronization
BottonPresTimeL01; BottonPresTimeR01;
SegIndL; SegIndR;

% data=rmfield(data,'BP');% try remove the field
% organize into data structure ***
for i=1:12
    data(1).BP{i}=BottonPresTimeL01(SegIndL(2*i-1):SegIndL(2*i),:); % in cells
    data(2).BP{i}=BottonPresTimeR01(SegIndR(2*i-1):SegIndR(2*i),:);
end

% add to mat file ***
BP=[];
BP(1,1).BP=data(1).BP;
BP(2,1).BP=data(2).BP;
save(filename_variables,'BP','-append');

%% segment feedbacks from the other  - synchronization
FeedbTimeL01; FeedbTimeR01;
SegIndL; SegIndR;

% organize into data structure ***
for i=1:12
    data(1).FB{i}=FeedbTimeL01(SegIndL(2*i-1):SegIndL(2*i),:); % in cells
    data(2).FB{i}=FeedbTimeR01(SegIndR(2*i-1):SegIndR(2*i),:);
end
 
% add to mat file ***
FB=[];
FB(1,1).FB=data(1).FB;
FB(2,1).FB=data(2).FB;
save(filename_variables,'FB','-append');

%% %% segment EMG - synchronization
filtered_EMGL; filtered_EMGR; % before hilbert transform
upL;upR; % after hilbert transform
SegIndL; SegIndR;

% organize into data structure ***
for i=1:12
    data(1).EMG_env{i}=upL(SegIndL(2*i-1):SegIndL(2*i),:); % in cells
    data(2).EMG_env{i}=upR(SegIndR(2*i-1):SegIndR(2*i),:);
    data(1).EMG{i}=filtered_EMGL(SegIndL(2*i-1):SegIndL(2*i),:); % in cells
    data(2).EMG{i}=filtered_EMGR(SegIndR(2*i-1):SegIndR(2*i),:);
end

% add to mat file ***
EMG=[];
EMG(1,1).EMG_env=data(1).EMG_env;
EMG(2,1).EMG_env=data(2).EMG_env;
EMG(1,1).EMG=data(1).EMG;
EMG(2,1).EMG=data(2).EMG;
save(filename_variables,'EMG','-append');

%%  segment Photocell - synchronization
locsL; locsR;
locsL01=zeros(1,length(timeL)); 
locsR01=zeros(1,length(timeR));
locsL01(locsL)=1; % plot(locsL01);
locsR01(locsR)=1; % plot(locsR01);
SegIndL; SegIndR;

% organize into data structure ***
for i=1:12
    data(1).FBlocs{i}=locsL01(SegIndL(2*i-1):SegIndL(2*i)); % in cells
    data(2).FBlocs{i}=locsR01(SegIndR(2*i-1):SegIndR(2*i));
end

% add to mat file ***
FBlocs=[];
FBlocs(1,1).FBlocs=data(1).FBlocs;
FBlocs(2,1).FBlocs=data(2).FBlocs;
save(filename_variables,'FBlocs','-append');

%% Save data structure (skip);
tic
filename_structure=['data_structure' num2str(seed) '.mat'];
save(filename_structure,'data');
toc
% takes about 45s

%% syncopation segmentation below (for 20220610_2P)

%% EEG - syncopation (for 20220610_2P)
% mixedsigL=mixedsigL';
% mixedsigR=mixedsigR';
% SegtimeIndL; SegtimeIndR;

EEGOpenEyeRestingL=mixedsigL(SegtimeIndL(1):SegtimeIndL(2),:);
EEGCloseEyeRestingL=mixedsigL(SegtimeIndL(3):SegtimeIndL(4),:);
EEGCondi1L=mixedsigL(SegtimeIndL(5):SegtimeIndL(6),:);
EEGResting2L=mixedsigL(SegtimeIndL(6):SegtimeIndL(7),:);
EEGCondi2L=mixedsigL(SegtimeIndL(8):SegtimeIndL(9),:);
EEGResting3L=mixedsigL(SegtimeIndL(9):SegtimeIndL(10),:);
EEGCondi3L=mixedsigL(SegtimeIndL(11):SegtimeIndL(12),:);
EEGResting4L=mixedsigL(SegtimeIndL(12):SegtimeIndL(13),:);
EEGCondi4L=mixedsigL(SegtimeIndL(14):end,:);

EEGOpenEyeRestingR=mixedsigR(SegtimeIndR(1):SegtimeIndR(2),:);
EEGCloseEyeRestingR=mixedsigR(SegtimeIndR(3):SegtimeIndR(4),:);
EEGCondi1R=mixedsigR(SegtimeIndR(5):SegtimeIndR(6),:);
EEGResting2R=mixedsigR(SegtimeIndR(6):SegtimeIndR(7),:);
EEGCondi2R=mixedsigR(SegtimeIndR(8):SegtimeIndR(9),:);
EEGResting3R=mixedsigR(SegtimeIndR(9):SegtimeIndR(10),:);
EEGCondi3R=mixedsigR(SegtimeIndR(11):SegtimeIndR(12),:);
EEGResting4R=mixedsigR(SegtimeIndR(12):SegtimeIndR(13),:);
EEGCondi4R=mixedsigR(SegtimeIndR(14):end,:);

save(['EEG' num2str(seed) '.mat'],'EEGCondi1L','EEGCondi2L','EEGCondi3L',...
    'EEGCondi4L',...
    'EEGCondi1R','EEGCondi2R','EEGCondi3R',...
    'EEGCondi4R');

%% BP - syncopation (for 20220610_2P)
BottonPresTimeL01; BottonPresTimeR01;

BPOpenEyeRestingL=BottonPresTimeL01(SegtimeIndL(1):SegtimeIndL(2),:);
BPCloseEyeRestingL=BottonPresTimeL01(SegtimeIndL(3):SegtimeIndL(4),:);
BPCondi1L=BottonPresTimeL01(SegtimeIndL(5):SegtimeIndL(6),:);
BPResting2L=BottonPresTimeL01(SegtimeIndL(6):SegtimeIndL(7),:);
BPCondi2L=BottonPresTimeL01(SegtimeIndL(8):SegtimeIndL(9),:);
BPResting3L=BottonPresTimeL01(SegtimeIndL(9):SegtimeIndL(10),:);
BPCondi3L=BottonPresTimeL01(SegtimeIndL(11):SegtimeIndL(12),:);
BPResting4L=BottonPresTimeL01(SegtimeIndL(12):SegtimeIndL(13),:);
BPCondi4L=BottonPresTimeL01(SegtimeIndL(14):end,:);

BPOpenEyeRestingR=BottonPresTimeR01(SegtimeIndR(1):SegtimeIndR(2),:);
BPCloseEyeRestingR=BottonPresTimeR01(SegtimeIndR(3):SegtimeIndR(4),:);
BPCondi1R=BottonPresTimeR01(SegtimeIndR(5):SegtimeIndR(6),:);
BPResting2R=BottonPresTimeR01(SegtimeIndR(6):SegtimeIndR(7),:);
BPCondi2R=BottonPresTimeR01(SegtimeIndR(8):SegtimeIndR(9),:);
BPResting3R=BottonPresTimeR01(SegtimeIndR(9):SegtimeIndR(10),:);
BPCondi3R=BottonPresTimeR01(SegtimeIndR(11):SegtimeIndR(12),:);
BPResting4R=BottonPresTimeR01(SegtimeIndR(12):SegtimeIndR(13),:);
BPCondi4R=BottonPresTimeR01(SegtimeIndR(14):end,:);

save(['BP' num2str(seed) '.mat'],'BPCondi1L','BPCondi2L','BPCondi3L',...
    'BPCondi4L',...
    'BPCondi1R','BPCondi2R','BPCondi3R',...
    'BPCondi4R');

%% FB - syncopation (for 20220610_2P)
FeedbTimeL01; FeedbTimeR01;

FBOpenEyeRestingL=FeedbTimeL01(SegtimeIndL(1):SegtimeIndL(2),:);
FBCloseEyeRestingL=FeedbTimeL01(SegtimeIndL(3):SegtimeIndL(4),:);
FBCondi1L=FeedbTimeL01(SegtimeIndL(5):SegtimeIndL(6),:);
FBResting2L=FeedbTimeL01(SegtimeIndL(6):SegtimeIndL(7),:);
FBCondi2L=FeedbTimeL01(SegtimeIndL(8):SegtimeIndL(9),:);
FBResting3L=FeedbTimeL01(SegtimeIndL(9):SegtimeIndL(10),:);
FBCondi3L=FeedbTimeL01(SegtimeIndL(11):SegtimeIndL(12),:);
FBResting4L=FeedbTimeL01(SegtimeIndL(12):SegtimeIndL(13),:);
FBCondi4L=FeedbTimeL01(SegtimeIndL(14):end,:);

FBOpenEyeRestingR=FeedbTimeR01(SegtimeIndR(1):SegtimeIndR(2),:);
FBCloseEyeRestingR=FeedbTimeR01(SegtimeIndR(3):SegtimeIndR(4),:);
FBCondi1R=FeedbTimeR01(SegtimeIndR(5):SegtimeIndR(6),:);
FBResting2R=FeedbTimeR01(SegtimeIndR(6):SegtimeIndR(7),:);
FBCondi2R=FeedbTimeR01(SegtimeIndR(8):SegtimeIndR(9),:);
FBResting3R=FeedbTimeR01(SegtimeIndR(9):SegtimeIndR(10),:);
FBCondi3R=FeedbTimeR01(SegtimeIndR(11):SegtimeIndR(12),:);
FBResting4R=FeedbTimeR01(SegtimeIndR(12):SegtimeIndR(13),:);
FBCondi4R=FeedbTimeR01(SegtimeIndR(14):end,:);

save(['FB' num2str(seed) '.mat'],'FBCondi1L','FBCondi2L','FBCondi3L',...
    'FBCondi4L',...
    'FBCondi1R','FBCondi2R','FBCondi3R',...
    'FBCondi4R');

%% EMG - syncopation (for 20220610_2P)
filtered_EMGL; filtered_EMGR; % before hilbert transform
upL;upR; % after hilbert transform

EMGOpenEyeRestingL=upL(SegtimeIndL(1):SegtimeIndL(2),:);
EMGCloseEyeRestingL=upL(SegtimeIndL(3):SegtimeIndL(4),:);
EMGCondi1L=upL(SegtimeIndL(5):SegtimeIndL(6),:);
EMGResting2L=upL(SegtimeIndL(6):SegtimeIndL(7),:);
EMGCondi2L=upL(SegtimeIndL(8):SegtimeIndL(9),:);
EMGResting3L=upL(SegtimeIndL(9):SegtimeIndL(10),:);
EMGCondi3L=upL(SegtimeIndL(11):SegtimeIndL(12),:);
EMGResting4L=upL(SegtimeIndL(12):SegtimeIndL(13),:);
EMGCondi4L=upL(SegtimeIndL(14):end,:);

EMGOpenEyeRestingR=upR(SegtimeIndR(1):SegtimeIndR(2),:);
EMGCloseEyeRestingR=upR(SegtimeIndR(3):SegtimeIndR(4),:);
EMGCondi1R=upR(SegtimeIndR(5):SegtimeIndR(6),:);
EMGResting2R=upR(SegtimeIndR(6):SegtimeIndR(7),:);
EMGCondi2R=upR(SegtimeIndR(8):SegtimeIndR(9),:);
EMGResting3R=upR(SegtimeIndR(9):SegtimeIndR(10),:);
EMGCondi3R=upR(SegtimeIndR(11):SegtimeIndR(12),:);
EMGResting4R=upR(SegtimeIndR(12):SegtimeIndR(13),:);
EMGCondi4R=upR(SegtimeIndR(14):end,:);

save(['EMG' num2str(seed) '.mat'],'EMGCondi1L','EMGCondi2L','EMGCondi3L',...
    'EMGCondi4L',...
    'EMGCondi1R','EMGCondi2R','EMGCondi3R',...
    'EMGCondi4R');

%% segment Photocell -  syncopation (for 20220610_2P)
locsL; locsR;
locsL01=zeros(1,length(timeL)); 
locsR01=zeros(1,length(timeR));
locsL01(locsL)=1; % plot(locsL01);
locsR01(locsR)=1; % plot(locsR01);

FBlocsOpenEyeRestingL=locsL01(SegtimeIndL(1):SegtimeIndL(2));
FBlocsCloseEyeRestingL=locsL01(SegtimeIndL(3):SegtimeIndL(4));
FBlocsCondi1L=locsL01(SegtimeIndL(5):SegtimeIndL(6));
FBlocsResting2L=locsL01(SegtimeIndL(6):SegtimeIndL(7));
FBlocsCondi2L=locsL01(SegtimeIndL(8):SegtimeIndL(9));
FBlocsResting3L=locsL01(SegtimeIndL(9):SegtimeIndL(10));
FBlocsCondi3L=locsL01(SegtimeIndL(11):SegtimeIndL(12));
FBlocsResting4L=locsL01(SegtimeIndL(12):SegtimeIndL(13));
FBlocsCondi4L=locsL01(SegtimeIndL(14):end);

FBlocsOpenEyeRestingR=locsR01(SegtimeIndR(1):SegtimeIndR(2));
FBlocsCloseEyeRestingR=locsR01(SegtimeIndR(3):SegtimeIndR(4));
FBlocsCondi1R=locsR01(SegtimeIndR(5):SegtimeIndR(6));
FBlocsResting2R=locsR01(SegtimeIndR(6):SegtimeIndR(7));
FBlocsCondi2R=locsR01(SegtimeIndR(8):SegtimeIndR(9));
FBlocsResting3R=locsR01(SegtimeIndR(9):SegtimeIndR(10));
FBlocsCondi3R=locsR01(SegtimeIndR(11):SegtimeIndR(12));
FBlocsResting4R=locsR01(SegtimeIndR(12):SegtimeIndR(13));
FBlocsCondi4R=locsR01(SegtimeIndR(14):end);

save(['FBlocs' num2str(seed) '.mat'],'FBlocsCondi1L','FBlocsCondi2L','FBlocsCondi3L',...
    'FBlocsCondi4L',...
    'FBlocsCondi1R','FBlocsCondi2R','FBlocsCondi3R',...
    'FBlocsCondi4R');

%% segment Pacers (for calculating errors) -  syncopation (for 20220610_2P)
PacerTimeIndL; PacerTimeIndR;
PacerTimeL01; PacerTimeR01;

PacerOpenEyeRestingL=PacerTimeL01(SegtimeIndL(1):SegtimeIndL(2));
PacerCloseEyeRestingL=PacerTimeL01(SegtimeIndL(3):SegtimeIndL(4));
PacerCondi1L=PacerTimeL01(SegtimeIndL(5):SegtimeIndL(6));
PacerResting2L=PacerTimeL01(SegtimeIndL(6):SegtimeIndL(7));
PacerCondi2L=PacerTimeL01(SegtimeIndL(8):SegtimeIndL(9));
PacerResting3L=PacerTimeL01(SegtimeIndL(9):SegtimeIndL(10));
PacerCondi3L=PacerTimeL01(SegtimeIndL(11):SegtimeIndL(12));
PacerResting4L=PacerTimeL01(SegtimeIndL(12):SegtimeIndL(13));
PacerCondi4L=PacerTimeL01(SegtimeIndL(14):end);

PacerOpenEyeRestingR=PacerTimeR01(SegtimeIndR(1):SegtimeIndR(2));
PacerCloseEyeRestingR=PacerTimeR01(SegtimeIndR(3):SegtimeIndR(4));
PacerCondi1R=PacerTimeR01(SegtimeIndR(5):SegtimeIndR(6));
PacerResting2R=PacerTimeR01(SegtimeIndR(6):SegtimeIndR(7));
PacerCondi2R=PacerTimeR01(SegtimeIndR(8):SegtimeIndR(9));
PacerResting3R=PacerTimeR01(SegtimeIndR(9):SegtimeIndR(10));
PacerCondi3R=PacerTimeR01(SegtimeIndR(11):SegtimeIndR(12));
PacerResting4R=PacerTimeR01(SegtimeIndR(12):SegtimeIndR(13));
PacerCondi4R=PacerTimeR01(SegtimeIndR(14):end);

save(['Pacer' num2str(seed) '.mat'],'PacerCondi1L','PacerCondi2L','PacerCondi3L',...
    'PacerCondi4L',...
    'PacerCondi1R','PacerCondi2R','PacerCondi3R',...
    'PacerCondi4R');

%% save BP and condiSeq for Syncopation experiment
seed;
% cd(['/home/hnl/Acquisition/zhibin/' num2str(last_seed) '_2P/Segmented_data/']); % on hnlstim2
cd(['/ssd/zhibin/1overf/' num2str(seed) '_2P/Segmented_data/']);
save(['BP' num2str(last_seed) '.mat'],'BP');
save(['condiSeq' num2str(last_seed) '.mat'],'condiSeq');

%% saves all variables from the current workspace
tic
save([num2str(seed) 'workspace.mat']);
toc
% takes about 11 min

%% Get variable names from mat file (faster)
load("20220721.mat")

tic
load([num2str(seed) 'workspace.mat']);
toc
% takes about 3 min
varlist1 = who();

% varlist1 = who('-file','20220811.mat');

tic
% matObj1 = matfile('20220721workspace.mat');
% varlist1 = who(matObj1);
toc

varlist2 = who();

varlist_diff = setdiff(varlist2,varlist1);

% clear varlist1

% append new variables to the old mat files
% tic
% varlist_diff_string=string(varlist_diff);
% for i=1:length(varlist_diff_string)
%     save('20220811workspace.mat',varlist_diff_string(i),'-append');
% end
% toc

tic 
varlist_diff_string=string(varlist_diff);
save([num2str(seed) 'workspace.mat'],varlist_diff_string{:},'-append');
toc
% takes about 1.5 min

% clear all new added variables
clear(varlist_diff{:})

%% moving files 

seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816];
cd /ssd/zhibin/1overf/
for s=1:8
    mkdir(['/ssd/zhibin/1overf/' num2str(seeds(s,:)) '_2P/Cleaned_data/Plots'])
    cd(['/ssd/zhibin/1overf/' num2str(seeds(s,:)) '_2P/Cleaned_data'])
    movefile cleanup*.png Plots
    cd /ssd/zhibin/1overf/
end

%%