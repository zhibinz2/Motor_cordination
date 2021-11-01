
C3EEG=reRef_data(:,28,1);
C4EEG=reRef_data(:,32,1);



r=xcorr(squeeze(C3EEG), squeeze(C3EEG));
plot(r);

r=xcorr(C3EEG, C3EEG);
plot(r);

r=xcorr(C3EEG);
plot(r);

subplot(1,2,1);
[r,lags]=xcorr(C3EEG, C3EEG);
plot(lags./2,r);xlabel('time [ms]');ylabel('cross correlation');xlim([-250 250]);
title('auto-correlation of C3 itself')

subplot(1,2,2);
[r,lags]=xcorr(C3EEG, C4EEG);
plot(lags./2,r);xlabel('time [ms]');ylabel('cross correlation');xlim([-250 250]);
title('cross-correlation between C3 and C4')
sgtitle('Subplot Grid Title') 
suptitle('I am a super title')


PlotXcorr(C3EEG,C4EEG,Fs);
testdataC3=C3EEG;
testdataC4=C3EEG;
%% try out xcorr exmaples
n = 0:15;
x = 0.84.^n;
y = circshift(x,5);
[c,lags] = xcorr(x,y);
stem(lags,c);


plot(x,'r');
hold on;
plot(y,'b');
plot(lags,'g')

%% MOTOR cortex
FC5EEG=reRef_data(:,18,1);
FC3EEG=reRef_data(:,19,1);
FC1EEG=reRef_data(:,20,1);
FC2EEG=reRef_data(:,22,1);
FC4EEG=reRef_data(:,23,1);
FC6EEG=reRef_data(:,24,1);
FCC5hEEG=reRef_data(:,81,1);
FCC3hEEG=reRef_data(:,82,1);
FCC1hEEG=reRef_data(:,83,1);
FCC2hEEG=reRef_data(:,84,1);
FCC4hEEG=reRef_data(:,85,1);
FCC6hEEG=reRef_data(:,86,1);
C5EEG=reRef_data(:,27,1);
C3EEG=reRef_data(:,28,1); % center of left motor cortex
C1EEG=reRef_data(:,29,1);
C2EEG=reRef_data(:,31,1);
C4EEG=reRef_data(:,32,1); % center of right motor cortex
C6EEG=reRef_data(:,33,1);
CCP5hEEG=reRef_data(:,101,1);
CCP3hEEG=reRef_data(:,102,1);
CCP1hEEG=reRef_data(:,103,1);
CCP2hEEG=reRef_data(:,104,1);
CCP4hEEG=reRef_data(:,105,1);
CCP6hEEG=reRef_data(:,106,1);
P5EEG=reRef_data(:,47,1);
P3EEG=reRef_data(:,48,1);
P1EEG=reRef_data(:,49,1);
P2EEG=reRef_data(:,51,1);
P4EEG=reRef_data(:,52,1);
P6EEG=reRef_data(:,53,1);
%% 

