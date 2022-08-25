%% Select segmentation indicies to use 
% Work for 20220713_2P; 20220721_2P; 20220810_2P; 20220815_2P

%  For synchronization
% Work for 20220713_2P; 20220810_2P:
% Select segmentation points based on pacers and phtotocells on light senor 1
% 50 Pacer markers selected (2 resting + 12 blocks + 11 resting)
% SegIndL=PacerTimeIndL([1:2 repelem([2+[32*([1:12]-1)]],4)+repmat([1 2 3 32],1,12)]);
% SegIndL1=PacerTimeIndL([4+[32*([1:12]-1)]+30]);% 12 segmentation points of the last pacer (30th pacer)
% SegIndR1=PacerTimeIndR([4+[32*([1:12]-1)]+30]);
SegIndL1=PacersL([6+[4*([1:12]-1)]]);% 12 segmentation points of the last pacer (30th pacer)
SegIndR1=PacersR([6+[4*([1:12]-1)]]);

% For syncopation 
% Work for 20220721_2P; 20220811_2P; 20220816_2P;
Seg30IndL1=PacerTimeIndL(SegPacerIndL([[1:12]*3+2]));% the 1st pacer
Seg30IndR1=PacerTimeIndR(SegPacerIndR([[1:12]*3+2]));% the 1st pacer
SegIndL1=PacerTimeIndL(SegPacerIndL([[1:12]*3+2])+29);% the 30th pacer
SegIndR1=PacerTimeIndR(SegPacerIndR([[1:12]*3+2])+29);% the 30th pacer
% Work for 20220808_2P;
Seg30IndL1=SegtimeIndL([[1:12]*3-2]);% the 1st pacer
Seg30IndR1=SegtimeIndR([[1:12]*3-2]);% the 1st pacer
SegIndL1=PacerTimeIndL(SegPacerIndL([[1:12]*3]-2)+29);% the 30th pacer
SegIndR1=PacerTimeIndR(SegPacerIndR([[1:12]*3]-2)+29);% the 30th pacer
% Examine 
% figure('units','normalized','outerposition',[0 0 1 0.6]);
canvas(0.5,0.25);
plot(PacerTimeIndL-PacerTimeIndL(1),2*ones(1,length(PacerTimeIndL)),'ro');hold on;
plot(SegIndL1-PacerTimeIndL(1),2.1*ones(length(SegIndL1),11,1),'ro'); 
plot(PacerTimeIndR-PacerTimeIndR(1),ones(1,length(PacerTimeIndR)),'bo');hold on;
plot(SegIndR1-PacerTimeIndR(1),1.1*ones(length(SegIndR1),11,1),'bo'); 
ylim([-3 6]);

% Select segmentation points of last BP in each trial
% determine a threshold of numbers of samples in between trials
if srL==srR; sr=srL; end; threshold = sr*30; 
% threshold = sr*15; % for 20220811_2P
% 12 segmentation points of the last tap in each trial
SegIndL2=[BottonPresTimeIndL(find([diff(BottonPresTimeIndL')]>threshold)); BottonPresTimeIndL(end)]; 
SegIndR2=[BottonPresTimeIndR(find([diff(BottonPresTimeIndR')]>threshold)); BottonPresTimeIndR(end)]; 
% Exam
% figure('units','normalized','outerposition',[0 0 1 0.6]);
canvas(0.5,0.25);
subplot(2,1,1);
plot(BottonPresTimeIndL,ones(1,length(BottonPresTimeIndL)),'ro');hold on;
plot(SegIndL2,ones(length(SegIndL2),11,1).*1.1,'ro'); ylim([0 2]);
subplot(2,1,2);
plot(BottonPresTimeIndR,ones(1,length(BottonPresTimeIndR)),'bo');hold on;
plot(SegIndR2,ones(length(SegIndR2),11,1).*1.1,'bo'); ylim([0 2]);
% For 20220810_2P, select the last 12 elements from L side
SegIndL2=SegIndL2(6:end);
% For 20220811_2P, remove the 3rd elements from R side
SegIndR2(3)=[];

% Combine last pacer and last tap indicies
% SegIndL=unique([SegIndL1;SegIndL2]);
% SegIndR=unique([SegIndR1;SegIndR2]);
SegIndL=sort([SegIndL1;SegIndL2]);
SegIndR=sort([SegIndR1;SegIndR2]);
% Examine 
% figure('units','normalized','outerposition',[0 0 1 0.6]);
canvas(0.5,0.25);
plot(PacerTimeIndL-PacerTimeIndL(1),2*ones(1,length(PacerTimeIndL)),'ro');hold on;
plot(SegIndL-PacerTimeIndL(1),2.1*ones(1,length(SegIndL)),'ro'); 
plot(PacerTimeIndR-PacerTimeIndR(1),ones(1,length(PacerTimeIndR)),'bo');hold on;
plot(SegIndR-PacerTimeIndR(1),1.1*ones(1,length(SegIndR)),'bo'); 
ylim([-2 5]);

% figure('units','normalized','outerposition',[0 0 1 0.6]);
canvas(0.5,0.25);
plot(PacerTimeIndL-PacerTimeIndL(1),2*ones(1,length(PacerTimeIndL)),'ro');hold on;
plot(SegIndL1-PacerTimeIndL(1),2.1*ones(1,length(SegIndL1)),'r.'); 
plot(SegIndL2-PacerTimeIndL(1),2.1*ones(1,length(SegIndL2)),'rx'); 
plot(PacerTimeIndR-PacerTimeIndR(1),ones(1,length(PacerTimeIndR)),'bo');hold on;
plot(SegIndR1-PacerTimeIndR(1),1.1*ones(1,length(SegIndR1)),'b.'); 
plot(SegIndR2-PacerTimeIndR(1),1.1*ones(1,length(SegIndR2)),'bx'); 
ylim([-1 4]);


%% Select segmentation indicies to use (202200804_2P)
%  For synchronization (20220804_2P):
% Select segmentation points based on pacers and phtotocells on light senor 1
% 50 Pacer markers selected (2 resting + 12 blocks + 11 resting)
SegIndL1=PacersL([[4+[4*([1:10]-1)]] 46 50]);% 12 segmentation points of the last pacer (30th pacer)
SegIndR1=PacersR([[4+[4*([1:10]-1)]] 46 50]);
% Examine 
figure('units','normalized','outerposition',[0 0 1 0.6]);
plot(finalPacerTimeIndL-finalPacerTimeIndL(1),2*ones(1,length(finalPacerTimeIndL)),'ro');hold on;
plot(SegIndL1-finalPacerTimeIndL(1),2.1*ones(length(SegIndL1),11,1),'ro'); 
plot(finalPacerTimeIndR-finalPacerTimeIndR(1),ones(1,length(finalPacerTimeIndR)),'bo');hold on;
plot(SegIndR1-finalPacerTimeIndR(1),1.1*ones(length(SegIndR1),11,1),'bo'); 
ylim([-3 6]);

% Select segmentation points of last BP in each trial
% determine a threshold of numbers of samples in between trials
if srL==srR; sr=srL; end; threshold = sr*30; 
% 12 segmentation points of the last tap in each trial
SegIndL2=[BottonPresTimeIndL(find([diff(BottonPresTimeIndL')]>threshold)); BottonPresTimeIndL(end)]; 
SegIndR2=[BottonPresTimeIndR(find([diff(BottonPresTimeIndR')]>threshold)); BottonPresTimeIndR(end)]; 
% Exam
figure('units','normalized','outerposition',[0 0 1 0.6]);
subplot(2,1,1);
plot(BottonPresTimeIndL,ones(1,length(BottonPresTimeIndL)),'ro');hold on;
plot(SegIndL2,ones(length(SegIndL2),11,1).*1.1,'ro'); ylim([0 2]);
subplot(2,1,2);
plot(BottonPresTimeIndR,ones(1,length(BottonPresTimeIndR)),'bo');hold on;
plot(SegIndR2,ones(length(SegIndR2),11,1).*1.1,'bo'); ylim([0 2]);
% remove the first ending BP of L and first 2 ending BP of R
SegIndL2=SegIndL2(2:end);
SegIndR2=SegIndR2(3:end);

% Combine last pacer and last tap indicies
% SegIndL=unique([SegIndL1;SegIndL2]);
% SegIndR=unique([SegIndR1;SegIndR2]);
SegIndL=sort([SegIndL1;SegIndL2]);
SegIndR=sort([SegIndR1;SegIndR2]);
% Move the last 2 trials to the front (for 20220804_2P)
SegIndL=SegIndL([end-3:end 1:end-4]);
SegIndR=SegIndR([end-3:end 1:end-4]);

% Examine 
figure('units','normalized','outerposition',[0 0 1 0.6]);
plot(finalPacerTimeIndL-finalPacerTimeIndL(1),2*ones(1,length(finalPacerTimeIndL)),'ro');hold on;
plot(SegIndL-finalPacerTimeIndL(1),2.1*ones(1,length(SegIndL)),'ro'); 
plot(finalPacerTimeIndR-finalPacerTimeIndR(1),ones(1,length(finalPacerTimeIndR)),'bo');hold on;
plot(SegIndR-finalPacerTimeIndR(1),1.1*ones(1,length(SegIndR)),'bo'); 
ylim([-2 5]);

% figure('units','normalized','outerposition',[0 0 1 0.6]);
canvas(0.5,0.25);
plot(finalPacerTimeIndL-finalPacerTimeIndL(1),2*ones(1,length(finalPacerTimeIndL)),'ro');hold on;
plot(SegIndL1-finalPacerTimeIndL(1),2.1*ones(1,length(SegIndL1))','r.'); 
plot(SegIndL2-finalPacerTimeIndL(1),2.1*ones(1,length(SegIndL2)),'rx'); 
plot(finalPacerTimeIndR-finalPacerTimeIndR(1),ones(1,length(finalPacerTimeIndR)),'bo');hold on;
plot(SegIndR1-finalPacerTimeIndR(1),1.1*ones(1,length(SegIndR1)),'b.'); 
plot(SegIndR2-finalPacerTimeIndR(1),1.1*ones(1,length(SegIndR2)),'bx'); 
legend('L pacers','L seg start','L seg end','R pacers','R seg start','R seg end')
ylim([-1 4]);

%% select the first 30 pacers intervals (synchronization)
% work for 202200713_2P;

PacersL;PacersR; % each has 50 elements (4+2*11+2*12=50)
Seg30IndL1=PacersL([5+[4*([1:12]-1)]]);% 12 segmentation points of the first pacer (30th pacer)
Seg30IndR1=PacersR([5+[4*([1:12]-1)]]);
Seg30IndL2=PacersL([6+[4*([1:12]-1)]]);% 12 segmentation points of the last pacer (30th pacer)
Seg30IndR2=PacersR([6+[4*([1:12]-1)]]);
% Combine first pacer and last pacer indicies
Seg30IndL=sort([Seg30IndL1;Seg30IndL2]);
Seg30IndR=sort([Seg30IndR1;Seg30IndR2]);
% Examine 
% figure('units','normalized','outerposition',[0 0 1 0.6]);
canvas(0.5,0.25);
plot(PacerTimeIndL-PacerTimeIndL(1),2*ones(1,length(PacerTimeIndL)),'ro');hold on;
plot(Seg30IndL1-PacerTimeIndL(1),2.1*ones(length(Seg30IndL1),11,1),'r.'); 
plot(Seg30IndL2-PacerTimeIndL(1),2.1*ones(1,length(SegIndL2)),'rx'); 
plot(PacerTimeIndR-PacerTimeIndR(1),ones(1,length(PacerTimeIndR)),'bo');hold on;
plot(Seg30IndR1-PacerTimeIndR(1),1.1*ones(length(Seg30IndR1),11,1),'b.'); 
plot(Seg30IndR2-PacerTimeIndR(1),1.1*ones(1,length(SegIndR2)),'bx'); 
ylim([-3 6]);

%% select the first 30 pacers intervals (syncopation)
% work for 20220808_2P; 20220811_2P; 202200816_2P;
Seg30IndL1;% the 1st pacer
Seg30IndR1;% the 1st pacer
SegIndL1;% the 30th pacer
SegIndR1;% the 30th pacer

% Combine first pacer and last pacer indicies
Seg30IndL=sort([Seg30IndL1;SegIndL1]);
Seg30IndR=sort([Seg30IndR1;SegIndR1]);
% Examine 
% figure('units','normalized','outerposition',[0 0 1 0.6]);
canvas(0.5,0.25);
plot(PacerTimeIndL-PacerTimeIndL(1),2*ones(1,length(PacerTimeIndL)),'ro');hold on;
plot(Seg30IndL1-PacerTimeIndL(1),2.1*ones(length(Seg30IndL1),11,1),'r.'); 
plot(SegIndL1-PacerTimeIndL(1),2.1*ones(1,length(SegIndL1)),'rx'); 
plot(PacerTimeIndR-PacerTimeIndR(1),ones(1,length(PacerTimeIndR)),'bo');hold on;
plot(Seg30IndR1-PacerTimeIndR(1),1.1*ones(length(Seg30IndR1),11,1),'b.'); 
plot(SegIndR1-PacerTimeIndR(1),1.1*ones(1,length(SegIndR1)),'bx'); 
ylim([-3 6]);

%% select the first 30 pacers intervals (202200804_2P)
PacersL;PacersR; % each has 50 elements (2*10+2*10+4+2+2+2=50)
Seg30IndL1=PacersL([[3+[4*([1:10]-1)]] 45 49]);% 12 segmentation points of the first pacer (30th pacer)
Seg30IndR1=PacersR([[3+[4*([1:10]-1)]] 45 49]);
Seg30IndL2=PacersL([[4+[4*([1:10]-1)]] 46 50]);% 12 segmentation points of the last pacer (30th pacer)
Seg30IndR2=PacersR([[4+[4*([1:10]-1)]] 46 50]);

% Combine first pacer and last pacer indicies
Seg30IndL=sort([Seg30IndL1;Seg30IndL2]);
Seg30IndR=sort([Seg30IndR1;Seg30IndR2]);
% Move the last 2 trials to the front (for 20220804_2P)
Seg30IndL=Seg30IndL([end-3:end 1:end-4]);
Seg30IndR=Seg30IndR([end-3:end 1:end-4]);
% Examine 
% figure('units','normalized','outerposition',[0 0 1 0.6]);
canvas(0.5,0.25);
plot(finalPacerTimeIndL-finalPacerTimeIndL(1),2*ones(1,length(finalPacerTimeIndL)),'ro');hold on;
plot(Seg30IndL1-finalPacerTimeIndL(1),2.1*ones(length(Seg30IndL1),11,1),'r.'); 
plot(Seg30IndL2-finalPacerTimeIndL(1),2.1*ones(1,length(SegIndL2)),'rx'); 
plot(finalPacerTimeIndR-finalPacerTimeIndR(1),ones(1,length(finalPacerTimeIndR)),'bo');hold on;
plot(Seg30IndR1-finalPacerTimeIndR(1),1.1*ones(length(Seg30IndR1),11,1),'b.'); 
plot(Seg30IndR2-finalPacerTimeIndR(1),1.1*ones(1,length(SegIndR2)),'bx'); 
ylim([-3 6]);


%% select the resting interval (maybe not, they are too noisy)
%% select the eye open/close resting eeg (later if needed)