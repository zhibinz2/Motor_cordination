cd C:\Users\NIRS\Desktop\USB documents\Matlab\matlab file converter
addpath Matlab

oxysoft2matlab

oxysoft2matlab ({'20211214.oxy5', 'Measurement Examples.oxyproj'})

%%
DataFNIR=ans;
DataFNIR=nirs_data;
RighthandData=ans;
LefthandData=ans;

oxyvals=DataFNIR.oxyvals;
dxyvals=DataFNIR.dxyvals;
time=DataFNIR.time;
SampleNo=DataFNIR.sampleNo;
wavelengths=DataFNIR.wavelengths;
DPF=DataFNIR.DPF;
diatance=DataFNIR.distance;
Fs=DataFNIR.Fs;
label=DataFNIR.label;
RxLabel=DataFNIR.RxLabel;
TxLabel=DataFNIR.TxLabel;
Rx_TxId=DataFNIR.Rx_TxId;
chanPos=DataFNIR.chanPos;
transPos=DataFNIR.transPos;
receiPos=DataFNIR.receiPos;
ADlabel=DataFNIR.ADlabel;
ADvalues=DataFNIR.ADvalues;
events=DataFNIR.events;
%% DataFNIR
oxyvalsR=DataFNIR.oxyvals;
dxyvalsR=DataFNIR.dxyvals;
timeR=DataFNIR.time;

oxyvalsL=DataFNIR.oxyvals;
dxyvalsL=DataFNIR.dxyvals;
timeL=DataFNIR.time;
%% nirs_data
oxyvalsR=nirs_data.oxyvals;
dxyvalsR=nirs_data.dxyvals;
timeR=nirs_data.time;

oxyvalsL=nirs_data.oxyvals;
dxyvalsL=nirs_data.dxyvals;
timeL=nirs_data.time;

%% Just plot as they are
figure;
subplot(2,1,1);
plotx(timeR,oxyvalsR,'r'); hold on;
plotx(timeR,dxyvalsR,'b'); hold off;
xlabel('time (second)');title('Right arm movements: oxyval (red); dxyvals (blue)');xlim([1 50]);
subplot(2,1,2);
plotx(timeL,oxyvalsL,'r'); hold on;
plotx(timeL,dxyvalsL,'b'); hold off;
xlabel('time (second)');title('Left arm movements: oxyval (red); dxyvals (blue)');xlim([1 50]);

Badchannels=[5 6];
Goodchannels=[1:1:22];
Goodchannels(Badchannels)=[];


%% plot baseline corrected
figure;
subplot(2,1,1);
oxyvalsRb=oxyvalsR-repmat(mean(oxyvalsR(50:100,:),1),2766,1);
dxyvalsRb=dxyvalsR-repmat(mean(dxyvalsR(50:100,:),1),2766,1);
plotx(timeR,oxyvalsRb(:,Goodchannels),'r'); hold on;
plotx(timeR,dxyvalsRb(:,Goodchannels),'b'); hold off;
xlabel('time (second)');title('Right arm movements: oxyval (red); dxyvals (blue)');xlim([1 50]);ylim([-2 2]);
subplot(2,1,2);
oxyvalsLb=oxyvalsL-repmat(mean(oxyvalsL(50:100,:),1),3216,1);
dxyvalsLb=dxyvalsL-repmat(mean(dxyvalsL(50:100,:),1),3216,1);
plotx(timeL,oxyvalsLb(:,Goodchannels),'r'); hold on;
plotx(timeL,dxyvalsLb(:,Goodchannels),'b'); hold off;
xlabel('time (second)');title('Left arm movements: oxyval (red); dxyvals (blue)');xlim([1 50]);ylim([-2 2]);