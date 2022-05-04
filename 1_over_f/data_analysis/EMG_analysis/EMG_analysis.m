% View raw EMG
close;figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(time,EMG2,'k');
%% EMG preprocessing 
% detrend
detrend_data=detrend(EMG2,2); % even better
hold on;
plot(time,detrend_data,'r'); % move to zero level
legend('raw-EMG','detrend-EMG');
% high pass filter % a little bit edge effect at the end
Hd = makefilter(sr,0.2,0.15,6,20,0); 
filtered_data=filtfilthd(Hd,detrend_data);
plot(time,filtered_data,'b');
legend('raw-EMG','detrend-EMG','high-pass-EMG');

%% Hilbert transform
% method 1 ****** 
EnveEMG=abs(hilbert(filtered_data));
close;figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(time,filtered_data,'b');
hold on; plot(time',EnveEMG,'r');

% method 2 ****** adjust Hilbert filter length
fl = 50*sr;
[up1,lo1] = envelope(filtered_data,fl,'analytic');
close;figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(time,filtered_data,'b');
hold on; plot(time,up1,'r');

% method 3 ******** adjust length of window
wl = 100;
[up1,lo1] = envelope(filtered_data,wl,'rms');
close;figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(time,filtered_data,'b');
hold on; plot(time,up1,'r');

% method 4 *******  spline interpolation over local maxima (best)
np = 1000; % N-tap Hilbert filter
[up1,lo1] = envelope(filtered_data,np,'peak');
close;figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(time,filtered_data,'b');
hold on; plot(time,up1,'r');
ylim([-1000 1000]);title(['np = ' num2str(np)]);
xlim([1000 1100]);
hold off;
%% plot all together and zoom in to examine
figure('units','normalized','outerposition',[0 0 1 0.6]);
subplot(3,1,1);
plot(Photocell);
ylabel('photocell signal');title('photocell signal');
hold on;
plot(BottonPresTimeInd,ones(1,length(BottonPresTimeInd)),'r.'); title('botton press');
legend({'photocell signal','botton press'});
subplot(3,1,2);
plot(up1,'k'); title('EMG'); % EMG channel
ylim([-1500 1500]);
subplot(3,1,3);
plot(filtered_data4);
xlabel('time (sample)');title('EEG');
ylim([-500 500]);

[x, ~] = ginput(2); % read two mouse clicks on the plot % x were index, y were real values
% get the proximate index
string(x)
ind1=round(x(1))
ind2=round(x(2))

subplot(3,1,1);
plot(Photocell');
ylabel('photocell signal');title('photocell signal');
hold on;
plot(BottonPresTimeInd,ones(1,length(BottonPresTimeInd)),'r.'); title('botton press');
legend({'photocell signal','botton press'});
xlim([ind1 ind2]);

subplot(3,1,2);
plot(up1,'r'); title('EMG'); % EMG channel
hold on; plot(filtered_data,'b');
xlim([ind1 ind2]);ylim([-1500 1500]);hold off;

subplot(3,1,3);
plotx(filtered_data4);
xlabel('time (sample)');title('EEG');
ylim([-500 500]);
xlim([ind1 ind2]);

suptitle(num2str(seed));

%% extract peaks from EMG evelops
% Examine peaks detection in analog1
[~, y] = ginput(1); % select a value for 'MinPeakProminence'
Halfhigh1=y;
% Check if need to adjust the Halfhigh cutoff
close;figure('units','normalized','outerposition',[0 0 1 0.3]);
findpeaks(up1,1:length(up1),'MinPeakProminence',Halfhigh1,'Annotate','extents');
yline(Halfhigh1,'m','MinPeakProminence');
ylim([min(up1),max(up1)+Halfhigh1]);

% locate the trial sessions % pks=value of the peak % locs=time of the peak
[pks,locs] = findpeaks(up1,time,'MinPeakProminence',Halfhigh1,'Annotate','extents');

%% plot all together and zoom in to examine
figure('units','normalized','outerposition',[0 0 1 0.6]);
subplot(3,1,1);
plot(Photocell);
ylabel('photocell signal');title('photocell signal');
hold on;
plot(BottonPresTimeInd,ones(1,length(BottonPresTimeInd)),'r.'); title('botton press');
legend({'photocell signal','botton press'});
xlim([0 4e6]);
subplot(3,1,2);
findpeaks(up1,1:length(up1),'MinPeakProminence',Halfhigh1,'Annotate','extents');
yline(Halfhigh1,'m','MinPeakProminence');
ylim([-500 1500]);
xlim([0 4e6]);
subplot(3,1,3);
plot(filtered_data4);
xlabel('time (sample)');title('EEG');
ylim([-500 500]);
xlim([0 4e6]);

[x, ~] = ginput(2); % read two mouse clicks on the plot % x were index, y were real values
% get the proximate index
string(x)
ind1=round(x(1))
ind2=round(x(2))

subplot(3,1,1);
plot(Photocell');
ylabel('photocell signal');title('photocell signal');
hold on;
plot(BottonPresTimeInd,ones(1,length(BottonPresTimeInd)),'r.'); title('botton press');
legend({'photocell signal','botton press'});
ylim([-500 1500]);
xlim([ind1 ind2]);

subplot(3,1,2);
findpeaks(up1,1:length(up1),'MinPeakProminence',Halfhigh1,'Annotate','extents');
title('EMG');
yline(Halfhigh1,'m','MinPeakProminence');
ylim([-500 1500]);
xlim([ind1 ind2]);

subplot(3,1,3);
plotx(filtered_data4);
xlabel('time (sample)');title('EEG');
ylim([-500 500]);
xlim([ind1 ind2]);

% zoom in on the plot in the horizontal direction.
% h = zoom;
% h.Motion = 'vertical';
% h.Enable = 'on';

suptitle([num2str(seed) ' RT condition']);

%% Compute error based on EMG
BottonPressTime;
BottonPresTimeInd;
EMGtime=locs;
[pks,locs] = findpeaks(up1,1:length(up1),'MinPeakProminence',Halfhigh1,'Annotate','extents');
EMGInd=locs;

% % Select a section to include only the Buttom press EMG
% [x, y] = ginput(2); 
% Halfhigh1=y;

% compute the error for condition 1
clear Error1 minValue closetIndex
Error1=[]; % for synchronization (condition 1)
StartStim1=4+240*(find(allPerm==1)-1)+1; % the starting photocell of condition 1
% find the time difference with the closest botton press 
for i=StartStim1:StartStim1+239 % i=StartStim1
% for i=StartStim1:length(locs) % i=StartStim1
    [minValue,closetIndex]=min(abs(EMGtime-PhotocellTime(i))); % closetIndex in BottonPressTime
    Error1(i-StartStim1+1)=EMGtime(closetIndex)-PhotocellTime(i);
end
plot(Error1,'r.');title('synchronization error');xlabel('taps');ylabel('timing error (s)');
% remove error > 1 s
indices=find(Error1>1);
Error1(indices)=NaN;
Error1(indices)=[];
% remove error < -1 s
indices=find(Error1<-1);
Error1(indices)=NaN;
Error1(indices)=[];
% Error1(indices)=2;


% compute the error for condition 2
clear Error2
Error2=[]; % for syncopation (condition 2)
StartStim2=4+240*(find(allPerm==2)-1)+1; % the starting photocell of condition 2
% find the time difference with the closest botton press 
for i=StartStim2:StartStim2+238 % i=StartStim2 % syncopation has one press less than condition 1 and 3
    Midpoints=(PhotocellTime(i)+PhotocellTime(i+1))/2;
    [minValue,closetIndex]=min(abs(EMGtime-Midpoints)); % closetIndex in BottonPressTime
    Error2(i-StartStim2+1)=EMGtime(closetIndex)-Midpoints;
end
plot(Error2,'r.');title('syncopation error');xlabel('taps');ylabel('timing error (s)');
% remove error > 1 s and < -1 s
indices=find(Error2>1);
Error2(indices)=NaN;
Error2(indices)=[];
% Error2(indices)=2;
indices=find(Error2<-1);
Error2(indices)=NaN;
Error2(indices)=[];
% Error2(indices)=2;

% compute the error for condition 3
clear Error3
Error3=[]; % for syncopation (condition 2)
StartStim3=4+240*(find(allPerm==3)-1)+1; % the starting photocell of condition 3
% find the time difference with the next botton press
for i=StartStim3:StartStim3+239 % i=StartStim3 % syncopation has one press less than condition 1 and 3
    diffInd=(find((EMGtime-PhotocellTime(i))>0));
    Error3(i-StartStim3+1)=EMGtime(diffInd(1))-PhotocellTime(i);
end
plot(Error3,'r.');title('randomization error');xlabel('taps');ylabel('timing error (s)');
% large errors are misses
% remove error > 1 s (In a simple RT task, anything longer than 1 second should be thrown out.)
indices=find(Error3>1);
Error3(indices)=NaN;
Error3(indices)=[];
% Error3(indices)=3;

