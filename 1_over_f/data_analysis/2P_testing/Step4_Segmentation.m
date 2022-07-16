%% go to saving directory
cd /ssd/zhibin/1overf/20220515_2P/Segmented_data/1_50Hz_ICAautomized
cd /ssd/zhibin/1overf/20220517_2P/Segmented_data/1_50Hz_ICAautomized

cd /ssd/zhibin/1overf/20220713_2P/Segmented_data/1.5_50Hz_ICAautomized

%% Collect all the variables
% experimental data extracted =====================================
allPerm;conditionNames;seed;
ifi;
% add to mat file ***
filename_variables=['data_variables' num2str(seed) '.mat'];
save(filename_variables,'allPerm');
save(filename_variables,'conditionNames','ifi','-append');
% organize into data structure ***
data=[]; % first dimension has to be 1 in order to explore field names in GUI
% data.condiSeq=allPerm;
% data.conditionNames=conditionNames; 
% data.ifi=ifi; 
% data(1,1).condiSeq=allPerm; data(2,1).condiSeq=allPerm;
% data(1,2).conditionNames=conditionNames; data(2,2).conditionNames=conditionNames;
for i=1:2 % first dimension has to be 1 in order to show field names when explore in GUI
    data(1,i).condiSeq=allPerm;
    data(1,i).condiNames=conditionNames;
    data(1,i).refreshRate=ifi;
end

% Recording variables extracted =====================================
timeL;samplesL;TRIGGERindL;srL;channels_infoL;labelsL;
timeR;samplesR;TRIGGERindR;srR;channels_infoR;labelsR;
if srL==srR; sr=srL; end;
% add to mat file ***
save(filename_variables,'sr','timeL','TRIGGERindL','srL','channels_infoL','labelsL',...
    'timeR','TRIGGERindR','srR','channels_infoR','labelsR','-append');
% organize into data structure ***
% data=[];
data(1,1).sr=srL; data(1,2).sr=srR;
% data(2,1).samples=samplesL; data(2,2).samples=samplesR;
% data=[];
% data(1,1).time=timeL; data(2,1).time=timeR;
% data(1,2).samples=samplesL; data(2,2).samples=samplesR;

% EEG variables extracted ==========================================
labels;% for plotting
chaninfo; % channel coordinates
mixedsigL;mixedsigR;
% add to mat file ***
save(filename_variables,'labels','chaninfo','-append');
% organize into data structure ***
data(1,1).labels=labels; data(1,2).labels=labels;
data(1,1).chaninfo=chaninfo; data(1,2).chaninfo=chaninfo;


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
save(filename_variables,'','-append');
% organize into data structure ***

% EMG variables needed =============================================
filtered_EMGL;filtered_EMGR;
loL;loR; % Hilbert envelop
% add to mat file ***
save(filename_variables,'','-append');
% organize into data structure ***

%% Select segmentation indicies to use
% Select segmentation points based on pacers and phtotocells on light senor 1
% 50 Pacer markers selected (2 resting + 12 blocks + 11 resting)
% SegIndL=PacerTimeIndL([1:2 repelem([2+[32*([1:12]-1)]],4)+repmat([1 2 3 32],1,12)]);
% 12 segmentation points of the last pacer (30th tap)
% SegIndL1=PacerTimeIndL([4+[32*([1:12]-1)]+30]);
% SegIndR1=PacerTimeIndR([4+[32*([1:12]-1)]+30]);
SegIndL1=PacersL([6+[4*([1:12]-1)]]);
SegIndR1=PacersR([6+[4*([1:12]-1)]]);

% Select segmentation points of last BP in each trial
% determine a threshold of numbers of samples in between trials
if srL==srR; sr=srL; end; threshold = sr*30; 
% 12 segmentation points of the last tap in each trial
SegIndL2=[BottonPresTimeIndL(find([diff(BottonPresTimeIndL')]>threshold)); BottonPresTimeIndL(end)]; 
SegIndR2=[BottonPresTimeIndR(find([diff(BottonPresTimeIndR')]>threshold)); BottonPresTimeIndR(end)]; 
figure('units','normalized','outerposition',[0 0 1 0.6]);
subplot(2,1,1);
plot(BottonPresTimeIndL,ones(1,length(BottonPresTimeIndL)),'ro');hold on;
plot(SegIndL2,ones(length(SegIndL2),11,1).*1.1,'ro'); ylim([0 2]);
subplot(2,1,2);
plot(BottonPresTimeIndR,ones(1,length(BottonPresTimeIndR)),'bo');hold on;
plot(SegIndR2,ones(length(SegIndR2),11,1).*1.1,'bo'); ylim([0 2]);

% Combine last pacer and last tap indicies
% SegIndL=unique([SegIndL1;SegIndL2]);
% SegIndR=unique([SegIndR1;SegIndR2]);
SegIndL=sort([SegIndL1;SegIndL2]);
SegIndR=sort([SegIndR1;SegIndR2]);

%% Segment EEG - synchronization
% mixedsigL=mixedsigL';
% mixedsigR=mixedsigR';
% SegIndL;
% SegIndR;
for i=1:12
    data(1,5,1,i)=mixedsigL(SegIndL(2*i-1):SegIndL(2*i),:);
    data(1,5,2,i)=mixedsigR(SegIndR(i));
data

% PacersL=PacerTimeIndL([1 2 3 4 5 34 35 36 37 66 67 68 69 98 99 100 101 130 131 132 133 162 163 164 165 194]);
% PacersR=PacerTimeIndR([1 2 3 4 5 34 35 36 37 66 67 68 69 98 99 100 101 130 131 132 133 162 163 164 165 194]);
EEGOpenEyeRestingL=mixedsigL(PacersL(1):PacersL(2),:);
EEGCloseEyeRestingL=mixedsigL(PacersL(3):PacersL(4),:);
EEGCondi1L=mixedsigL(PacersL(6):PacersL(7),:);
EEGResting2L=mixedsigL(PacersL(7):PacersL(8),:);
EEGCondi2L=mixedsigL(PacersL(10):PacersL(11),:);
EEGResting3L=mixedsigL(PacersL(11):PacersL(12),:);
EEGCondi3L=mixedsigL(PacersL(14):PacersL(15),:);
EEGResting4L=mixedsigL(PacersL(15):PacersL(16),:);
EEGCondi4L=mixedsigL(PacersL(18):PacersL(19),:);
EEGResting5L=mixedsigL(PacersL(19):PacersL(20),:);
EEGCondi5L=mixedsigL(PacersL(22):PacersL(23),:);% EEGCondi5L=mixedsigL(PacersL(22):end,:); % if 5 conditions
EEGResting6L=mixedsigL(PacersL(23):PacersL(24),:);
EEGCondi6L=mixedsigL(PacersL(26):end,:);

EEGOpenEyeRestingR=mixedsigR(PacersR(1):PacersR(2),:);
EEGCloseEyeRestingR=mixedsigR(PacersR(3):PacersR(4),:);
EEGCondi1R=mixedsigR(PacersR(6):PacersR(7),:);
EEGResting2R=mixedsigR(PacersR(7):PacersR(8),:);
EEGCondi2R=mixedsigR(PacersR(10):PacersR(11),:);
EEGResting3R=mixedsigR(PacersR(11):PacersR(12),:);
EEGCondi3R=mixedsigR(PacersR(14):PacersR(15),:);
EEGResting4R=mixedsigR(PacersR(15):PacersR(16),:);
EEGCondi4R=mixedsigR(PacersR(18):PacersR(19),:);
EEGResting5R=mixedsigR(PacersR(19):PacersR(20),:);
EEGCondi5R=mixedsigR(PacersR(22):PacersR(23),:); % EEGCondi5R=mixedsigR(PacersR(22):end,:);
EEGResting6R=mixedsigR(PacersR(23):PacersR(24),:);
EEGCondi6R=mixedsigR(PacersR(26):end,:);

save(['EEG' num2str(seed) '.mat'],'EEGCondi1L','EEGCondi2L','EEGCondi3L',...
    'EEGCondi4L','EEGCondi5L','EEGCondi6L',...
    'EEGCondi1R','EEGCondi2R','EEGCondi3R',...
    'EEGCondi4R','EEGCondi5R','EEGCondi6R');

% if 5 conditions
save(['EEG' num2str(seed) '.mat'],'EEGCondi1L','EEGCondi2L','EEGCondi3L',...
    'EEGCondi4L','EEGCondi5L',...
    'EEGCondi1R','EEGCondi2R','EEGCondi3R',...
    'EEGCondi4R','EEGCondi5R');

%% segment button presses - synchronization
BottonPresTimeL01; BottonPresTimeR01;

BPOpenEyeRestingL=BottonPresTimeL01(PacersL(1):PacersL(2),:);
BPCloseEyeRestingL=BottonPresTimeL01(PacersL(3):PacersL(4),:);
BPCondi1L=BottonPresTimeL01(PacersL(6):PacersL(7),:);
BPResting2L=BottonPresTimeL01(PacersL(7):PacersL(8),:);
BPCondi2L=BottonPresTimeL01(PacersL(10):PacersL(11),:);
BPResting3L=BottonPresTimeL01(PacersL(11):PacersL(12),:);
BPCondi3L=BottonPresTimeL01(PacersL(14):PacersL(15),:);
BPResting4L=BottonPresTimeL01(PacersL(15):PacersL(16),:);
BPCondi4L=BottonPresTimeL01(PacersL(18):PacersL(19),:);
BPResting5L=BottonPresTimeL01(PacersL(19):PacersL(20),:);
BPCondi5L=BottonPresTimeL01(PacersL(22):PacersL(23),:); % BPCondi5L=BottonPresTimeL01(PacersL(22):end,:); % if 5 conditions
BPResting6L=BottonPresTimeL01(PacersL(23):PacersL(24),:);
BPCondi6L=BottonPresTimeL01(PacersL(26):end,:);

BPOpenEyeRestingR=BottonPresTimeR01(PacersR(1):PacersR(2),:);
BPCloseEyeRestingR=BottonPresTimeR01(PacersR(3):PacersR(4),:);
BPCondi1R=BottonPresTimeR01(PacersR(6):PacersR(7),:);
BPResting2R=BottonPresTimeR01(PacersR(7):PacersR(8),:);
BPCondi2R=BottonPresTimeR01(PacersR(10):PacersR(11),:);
BPResting3R=BottonPresTimeR01(PacersR(11):PacersR(12),:);
BPCondi3R=BottonPresTimeR01(PacersR(14):PacersR(15),:);
BPResting4R=BottonPresTimeR01(PacersR(15):PacersR(16),:);
BPCondi4R=BottonPresTimeR01(PacersR(18):PacersR(19),:);
BPResting5R=BottonPresTimeR01(PacersR(19):PacersR(20),:);
BPCondi5R=BottonPresTimeR01(PacersR(22):PacersR(23),:); % BPCondi5R=BottonPresTimeR01(PacersR(22):end,:); % if 5 conditions
BPResting6R=BottonPresTimeR01(PacersR(23):PacersR(24),:);
BPCondi6R=BottonPresTimeR01(PacersR(26):end,:);

save(['BP' num2str(seed) '.mat'],'BPCondi1L','BPCondi2L','BPCondi3L',...
    'BPCondi4L','BPCondi5L','BPCondi6L',...
    'BPCondi1R','BPCondi2R','BPCondi3R',...
    'BPCondi4R','BPCondi5R','BPCondi6R');
% if 5 conditions
save(['BP' num2str(seed) '.mat'],'BPCondi1L','BPCondi2L','BPCondi3L',...
    'BPCondi4L','BPCondi5L',...
    'BPCondi1R','BPCondi2R','BPCondi3R',...
    'BPCondi4R','BPCondi5R');

%% segment feedbacks from the other  - synchronization
FeedbTimeL01; FeedbTimeR01;

FBOpenEyeRestingL=FeedbTimeL01(PacersL(1):PacersL(2),:);
FBCloseEyeRestingL=FeedbTimeL01(PacersL(3):PacersL(4),:);
FBCondi1L=FeedbTimeL01(PacersL(6):PacersL(7),:);
FBResting2L=FeedbTimeL01(PacersL(7):PacersL(8),:);
FBCondi2L=FeedbTimeL01(PacersL(10):PacersL(11),:);
FBResting3L=FeedbTimeL01(PacersL(11):PacersL(12),:);
FBCondi3L=FeedbTimeL01(PacersL(14):PacersL(15),:);
FBResting4L=FeedbTimeL01(PacersL(15):PacersL(16),:);
FBCondi4L=FeedbTimeL01(PacersL(18):PacersL(19),:);
FBResting5L=FeedbTimeL01(PacersL(19):PacersL(20),:);
FBCondi5L=FeedbTimeL01(PacersL(22):PacersL(23),:); % FBCondi5L=FeedbTimeL01(PacersL(22):end,:); % if 5 conditions
FBResting6L=FeedbTimeL01(PacersL(23):PacersL(24),:);
FBCondi6L=FeedbTimeL01(PacersL(26):end,:);

FBOpenEyeRestingR=FeedbTimeR01(PacersR(1):PacersR(2),:);
FBCloseEyeRestingR=FeedbTimeR01(PacersR(3):PacersR(4),:);
FBCondi1R=FeedbTimeR01(PacersR(6):PacersR(7),:);
FBResting2R=FeedbTimeR01(PacersR(7):PacersR(8),:);
FBCondi2R=FeedbTimeR01(PacersR(10):PacersR(11),:);
FBResting3R=FeedbTimeR01(PacersR(11):PacersR(12),:);
FBCondi3R=FeedbTimeR01(PacersR(14):PacersR(15),:);
FBResting4R=FeedbTimeR01(PacersR(15):PacersR(16),:);
FBCondi4R=FeedbTimeR01(PacersR(18):PacersR(19),:);
FBResting5R=FeedbTimeR01(PacersR(19):PacersR(20),:);
FBCondi5R=FeedbTimeR01(PacersR(22):PacersR(23),:); % FBCondi5R=FeedbTimeR01(PacersR(22):end,:); % if 5 conditions
FBResting6R=FeedbTimeR01(PacersR(23):PacersR(24),:);
FBCondi6R=FeedbTimeR01(PacersR(26):end,:);

save(['FB' num2str(seed) '.mat'],'FBCondi1L','FBCondi2L','FBCondi3L',...
    'FBCondi4L','FBCondi5L','FBCondi6L',...
    'FBCondi1R','FBCondi2R','FBCondi3R',...
    'FBCondi4R','FBCondi5R','FBCondi6R');
% if 5 conditions
save(['FB' num2str(seed) '.mat'],'FBCondi1L','FBCondi2L','FBCondi3L',...
    'FBCondi4L','FBCondi5L',...
    'FBCondi1R','FBCondi2R','FBCondi3R',...
    'FBCondi4R','FBCondi5R');

%% %% segment EMG - synchronization
filtered_EMGL; filtered_EMGR; % before hilbert transform
upL;upR; % after hilbert transform

EMGOpenEyeRestingL=upL(PacersL(1):PacersL(2),:);
EMGCloseEyeRestingL=upL(PacersL(3):PacersL(4),:);
EMGCondi1L=upL(PacersL(6):PacersL(7),:);
EMGResting2L=upL(PacersL(7):PacersL(8),:);
EMGCondi2L=upL(PacersL(10):PacersL(11),:);
EMGResting3L=upL(PacersL(11):PacersL(12),:);
EMGCondi3L=upL(PacersL(14):PacersL(15),:);
EMGResting4L=upL(PacersL(15):PacersL(16),:);
EMGCondi4L=upL(PacersL(18):PacersL(19),:);
EMGResting5L=upL(PacersL(19):PacersL(20),:);
EMGCondi5L=upL(PacersL(22):PacersL(23),:); % EMGCondi5L=upL(PacersL(22):end,:);  %  if 5 conditions
EMGResting6L=upL(PacersL(23):PacersL(24),:);
EMGCondi6L=upL(PacersL(26):end,:);

EMGOpenEyeRestingR=upR(PacersR(1):PacersR(2),:);
EMGCloseEyeRestingR=upR(PacersR(3):PacersR(4),:);
EMGCondi1R=upR(PacersR(6):PacersR(7),:);
EMGResting2R=upR(PacersR(7):PacersR(8),:);
EMGCondi2R=upR(PacersR(10):PacersR(11),:);
EMGResting3R=upR(PacersR(11):PacersR(12),:);
EMGCondi3R=upR(PacersR(14):PacersR(15),:);
EMGResting4R=upR(PacersR(15):PacersR(16),:);
EMGCondi4R=upR(PacersR(18):PacersR(19),:);
EMGResting5R=upR(PacersR(19):PacersR(20),:);
EMGCondi5R=upR(PacersR(22):PacersR(23),:); % EMGCondi5R=upR(PacersR(22):end,:);  % if 5 conditions
EMGResting6R=upR(PacersR(23):PacersR(24),:);
EMGCondi6R=upR(PacersR(26):end,:);

save(['EMG' num2str(seed) '.mat'],'EMGCondi1L','EMGCondi2L','EMGCondi3L',...
    'EMGCondi4L','EMGCondi5L','EMGCondi6L',...
    'EMGCondi1R','EMGCondi2R','EMGCondi3R',...
    'EMGCondi4R','EMGCondi5R','EMGCondi6R');
% if 5 conditions
save(['EMG' num2str(seed) '.mat'],'EMGCondi1L','EMGCondi2L','EMGCondi3L',...
    'EMGCondi4L','EMGCondi5L',...
    'EMGCondi1R','EMGCondi2R','EMGCondi3R',...
    'EMGCondi4R','EMGCondi5R');

%%  segment Photocell - synchronization
locsL; locsR;
locsL01=zeros(1,length(timeL)); 
locsR01=zeros(1,length(timeL));
locsL01(locsL)=1; % plot(locsL01);
locsR01(locsR)=1; % plot(locsR01);

FBlocsOpenEyeRestingL=locsL01(PacersL(1):PacersL(2));
FBlocsCloseEyeRestingL=locsL01(PacersL(3):PacersL(4));
FBlocsCondi1L=locsL01(PacersL(6):PacersL(7));
FBlocsResting2L=locsL01(PacersL(7):PacersL(8));
FBlocsCondi2L=locsL01(PacersL(10):PacersL(11));
FBlocsResting3L=locsL01(PacersL(11):PacersL(12));
FBlocsCondi3L=locsL01(PacersL(14):PacersL(15));
FBlocsResting4L=locsL01(PacersL(15):PacersL(16));
FBlocsCondi4L=locsL01(PacersL(18):PacersL(19));
FBlocsResting5L=locsL01(PacersL(19):PacersL(20));
FBlocsCondi5L=locsL01(PacersL(22):PacersL(23)); % FBlocsCondi5L=locsL01(PacersL(22):end); % if 5 conditions
FBlocsResting6L=locsL01(PacersL(23):PacersL(24));
FBlocsCondi6L=locsL01(PacersL(26):end);

FBlocsOpenEyeRestingR=locsR01(PacersR(1):PacersR(2));
FBlocsCloseEyeRestingR=locsR01(PacersR(3):PacersR(4));
FBlocsCondi1R=locsR01(PacersR(6):PacersR(7));
FBlocsResting2R=locsR01(PacersR(7):PacersR(8));
FBlocsCondi2R=locsR01(PacersR(10):PacersR(11));
FBlocsResting3R=locsR01(PacersR(11):PacersR(12));
FBlocsCondi3R=locsR01(PacersR(14):PacersR(15));
FBlocsResting4R=locsR01(PacersR(15):PacersR(16));
FBlocsCondi4R=locsR01(PacersR(18):PacersR(19));
FBlocsResting5R=locsR01(PacersR(19):PacersR(20));
FBlocsCondi5R=locsR01(PacersR(22):PacersR(23)); % FBlocsCondi5R=locsR01(PacersR(22):end); % if 5 conditions
FBlocsResting6R=locsR01(PacersR(23):PacersR(24));
FBlocsCondi6R=locsR01(PacersR(26):end);

save(['FBlocs' num2str(seed) '.mat'],'FBlocsCondi1L','FBlocsCondi2L','FBlocsCondi3L',...
    'FBlocsCondi4L','FBlocsCondi5L','FBlocsCondi6L',...
    'FBlocsCondi1R','FBlocsCondi2R','FBlocsCondi3R',...
    'FBlocsCondi4R','FBlocsCondi5R','FBlocsCondi6R');
% if 5 conditions
save(['FBlocs' num2str(seed) '.mat'],'FBlocsCondi1L','FBlocsCondi2L','FBlocsCondi3L',...
    'FBlocsCondi4L','FBlocsCondi5L',...
    'FBlocsCondi1R','FBlocsCondi2R','FBlocsCondi3R',...
    'FBlocsCondi4R','FBlocsCondi5R');

%% EEG - syncopation
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

%% BP - syncopation
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

%% FB - syncopation
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

%% EMG - syncopation
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

%% segment Photocell -  syncopation
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

%% segment Pacers (for calculating errors) -  syncopation
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