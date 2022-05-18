A=[0 0 0 0 1 1 0 0 1 1 1 1 0 0 0 0 1 1 1 1 1 0 0 0 ];
pressInd=find([0 diff(A)]==1)
pressIntervals=[pressInd(1) diff(pressInd)]

%% extrat time points for key presses, pacers,feedbacks


% look for values 
PresInd=find(A == 1); % extract Index of real key presses in the values
plot(PresInd,ones(1,length(PresInd)),'ro'); % look at the above Index (one press produced several indices)
ifi=1/100;waitframes=2;sr=2000;
threshold = ifi*waitframes*sr; % determine a threshold of numbers of frames in the button press interval
BottonPresTimeInd=PresInd(find([1 diff(PresInd')>threshold])); % exact index of key press onset in datatimes (reduce several indices into one)

% create a time series that assign botton presses as 1, all other as 0
BottonPresTime01=zeros(size(time));
BottonPresTime01(BottonPresTimeInd)=1;
plot(BottonPresTime01)

%% reset time in play R 
