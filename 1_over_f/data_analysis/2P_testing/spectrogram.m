%% 1overf of EEG - wavelet beta/alpha
wfreq = [2 4 6 8 10 14 18 24 30 40];
wfc = 1.5;
sr=2000;

UniCondi=unique(CondiData);

AllchanNames=labels;

%% Padding 500 ms of zeros before and after baselinecorrected_trial
EEGseries=EEGCondi1L;
padding=zeros(1000,size(EEGseries,2));
EEGseries_paddings=cat(1,padding,EEGseries,padding);

%%
figure;
subplot(1,2,1)
plot(Calinterval(BP_L')./2,'r');
hold on;
plot(Calinterval(BP_R')./2,'b');