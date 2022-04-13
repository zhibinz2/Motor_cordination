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
np = 100;
[up1,lo1] = envelope(filtered_data,np,'peak');
close;figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(time,filtered_data,'b');
hold on; plot(time,up1,'r');
ylim([-1000 1000]);title(['np = ' num2str(np)]);
xlim([1000 1100]);

%% plot all together
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
