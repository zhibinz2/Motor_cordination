%% tryout with calculation
A=[0 0 0 0 1 1 0 0 1 1 1 1 0 0 0 0 1 1 1 1 1 0 0 0 ];
pressInd=find([0 diff(A)]==1)
pressIntervals=[pressInd(1) diff(pressInd)]

%% extract time points for key presses (maybe create a function to do the same for pacers,feedbacks)
% look for values 
PresInd=find(TriggersL==239); % extract Index of real key presses in the values
plot(PresInd,ones(1,length(PresInd)),'ro'); % look at the above Index (one press produced several indices)
% ifi=1/100;waitframes=2;sr=2000;
% determine a threshold of numbers of frames in the button press interval
threshold = 1*1.5; % interval larger than 1 consecutive samples
% button presses recorded via DP port into stimtracker is sampled at 2000 Hz
% feedback via photocell to stimtracker has minial interval every other waitframes (minmal of 1000 Hz)
BottonPresTimeInd=PresInd(find([0 diff(PresInd')]>threshold)); % exact index of key press onset in datatimes (reduce several indices into one)

% create a time series that assign botton presses as 1, all other as 0
BottonPresTimeL01=zeros(size(timeL));
BottonPresTimeL01(BottonPresTimeInd)=1;
plot(BottonPresTimeL01)

%% extract time points for feedbacks
% look for values 
FeedbInd=find(TriggersR==251); % extract Index of real key presses in the values
plot(FeedbInd,ones(1,length(FeedbInd)),'ro'); % look at the above Index (one press produced several indices)
% ifi=1/100;waitframes=2;sr=2000;
% determine a threshold of numbers of frames in the button press interval
threshold = ifi*waitframes*srR*1.5; % number of samples within a interval larger than the waitframes
% button presses recorded via DP port into stimtracker is sampled at 2000 Hz
% feedback via photocell to stimtracker has minial interval every other waitframes (minmal of 1000 Hz)
FeedbTimeInd=FeedbInd(find([0 diff(FeedbInd')]>threshold)); % exact index of key press onset in datatimes (reduce several indices into one)

% create a time series that assign botton presses as 1, all other as 0
FeedbTimeL01=zeros(size(timeL));
FeedbTimeL01(FeedbTimeInd)=1;
plot(FeedbTimeL01)


%% extract time points for pacers
% look for values 
PacerInd=find(TriggersL==127); % extract Index of real key presses in the values
plot(PacerInd,ones(1,length(PacerInd)),'ro'); % look at the above Index (one press produced several indices)
% Pacer at 2Hz has intervals of 500ms (1000 samples)
% Pacer at 3Hz has interval of 333ms (666.6 samples)
% determine a threshold of numbers of samples in pacer interval
threshold = srL/3*0.75; % number of samples within a interval larger than the waitframes
% button presses recorded via DP port into stimtracker is sampled at 2000 Hz
% feedback via photocell to stimtracker has minial interval every other waitframes (minmal of 1000 Hz)
PacerTimeInd=PacerInd(find([0 diff(PacerInd')]>threshold)); % exact index of key press onset in datatimes (reduce several indices into one)
% 2+2+30*6+2*5 = 194

% create a time series that assign botton presses as 1, all other as 0
PacerTimeL01=zeros(size(timeL));
PacerTimeL01(PacerTimeInd)=1;
plot(PacerTimeL01);
% xlim([-1 length(PacerTimeL01)]);

%% Reset time in play L and play R in allignment
