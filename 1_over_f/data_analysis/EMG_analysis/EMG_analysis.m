% View raw EMG
close;figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(time,EMG,'k');
%% EMG preprocessing 
% detrend
detrend_data=detrend(EMG,2); % even better
hold on;
plot(time,detrend_data,'r'); % move to zero level
legend('raw-EMG','detrend-EMG');
% high pass filter
Hd = makefilter(sr,0.2,0.15,6,20,0); 
filtered_data=filtfilthd(Hd,detrend_data);
plot(time,filtered_data,'b');
legend('raw-EMG','detrend-EMG','high-pass-EMG');

%% Hilbert transform
EnveEMG=abs(hilbert(filtered_data));

close;figure('units','normalized','outerposition',[0 0 1 0.3]);
plot(time,filtered_data,'b');
hold on; plot(time',EnveEMG,'r');



