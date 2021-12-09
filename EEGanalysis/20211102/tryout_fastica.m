cd /home/zhibin/Documents/Acquisition/bimanual_Reach_Jack_2021111802/
EEGfileName='bimanual_Reach_Jack_20211118'; EEG=loadcurry([pwd '/' EEGfileName '.cdt']);
cd /home/zhibin/Documents/GitHub/Motor_cordination/EEGanalysis/20211102
tic
run organize_EEG_filter_step2.m

[icasig, A, W] = fastica(EEGdata);
% Each row of matrix mixedsig is one observed signal.
% the rows of icasig contain the estimated independent components.


mixedsig=EEGdata;


[icasig] = fastica(mixedsig); % the rows of icasig contain the
% estimated independent components.
%
[icasig, A, W] = fastica(mixedsig); % outputs the estimated separating
% matrix W and the corresponding mixing matrix A. 
% Each row of A and W corresponding to one channel
% Each colum of A and W correponding to the weight distribution of one componenet

[A, W] = fastica(mixedsig); % gives only the estimated mixing matrix
% A and the separating matrix W.

EEGdata=A*icasig;
EEGdataRestore=A*icasig;

icasig=W*EEGdata;
icasigRestore=W*EEGdata;

%%
fasticag
%%

figure;
plotx(icasig(:,10000:20000));
figure;
plotx(icasigRestore(:,10000:20000));

figure;
plotx(EEGdata(:,10000:20000));
figure;
plotx(EEGdataRestore(:,10000:20000));


imagesc(A);
topoplot(A(1,:),EEG.chanlocs(1:128),'conv','on');
topoplot(A(:,1),EEG.chanlocs(1:128),'conv','on');

imagesc(W);
topoplot(W(1,:),EEG.chanlocs(1:128),'conv','on');
topoplot(W(:,1),EEG.chanlocs(1:128),'conv','on');

imagesc(A);
topoplot(A(1,:),EEG.chanlocs(1:128),'conv','on');
topoplot(A(:,1),EEG.chanlocs(1:128),'conv','on');

Chanlocs_topolot=EEG.chanlocs(1:128);

topoplot(A(1,:),Neuroscan_topoplot_chans);

topoplot(A(1,:),chanlocs(1:128),'conv','on');
%%
% 139 channel labels
chanlocs = struct('labels', {'FP1' 'FPZ' 'FP2' 'AF3' 'AF4' 'F11' 'F7' 'F5' 'F3' 'F1' 'FZ' 'F2' 'F4' 'F6' 'F8' 'F12' 'FT11' 'FC5' 'FC3' 'FC1' 'FCZ' 'FC2' 'FC4' 'FC6' 'FT12' 'T7' 'C5' 'C3' 'C1' 'CZ' 'C2' 'C4' 'C6' 'T8' 'TP7' 'CP5' 'CP3' 'CP1' 'CPZ' 'CP2' 'CP4' 'CP6' 'TP8' 'M1' 'M2' 'P7' 'P5' 'P3' 'P1' 'PZ' 'P2' 'P4' 'P6' 'P8' 'PO7' 'PO3' 'POZ' 'PO4' 'PO8' 'O1' 'OZ' 'O2' 'CB1' 'CB2' 'AFP1' 'AFP2' 'AF7' 'AF5' 'AFZ' 'AF6' 'AF8' 'AFF5H' 'AFF3H' 'AFF1H' 'AFF2H' 'AFF4H' 'AFF6H' 'F9' 'F10' 'FFT7H' 'FFC5H' 'FFC3H' 'FFC1H' 'FFC2H' 'FFC4H' 'FFC6H' 'FFT8H' 'FT9' 'FT7' 'FT8' 'FT10' 'FTT7H' 'FCC5H' 'FCC3H' 'FCC1H' 'FCC2H' 'FCC4H' 'FCC6H' 'FTT8H' 'TTP7H' 'CCP5H' 'CCP3H' 'CCP1H' 'CCP2H' 'CCP4H' 'CCP6H' 'TTP8H' 'TPP7H' 'CPP5H' 'CPP3H' 'CPP1H' 'CPP2H' 'CPP4H' 'CPP6H' 'TPP8H' 'P9' 'P10' 'PPO3H' 'PPO1H' 'PPO2H' 'PPO4H' 'PO9' 'PO5' 'PO1' 'PO2' 'PO6' 'PO10' 'CBZ' 'VEOG' 'HEOG' 'EMG1' 'EMG2' 'HL 1' 'HL 2' 'EMG3' 'EMG4' 'EMG5' 'EMG6' 'TRIGGER'});
% 128 channel labels
chanlocs = struct('labels', {'FP1' 'FPZ' 'FP2' 'AF3' 'AF4' 'F11' 'F7' 'F5' 'F3' 'F1' 'FZ' 'F2' 'F4' 'F6' 'F8' 'F12' 'FT11' 'FC5' 'FC3' 'FC1' 'FCZ' 'FC2' 'FC4' 'FC6' 'FT12' 'T7' 'C5' 'C3' 'C1' 'CZ' 'C2' 'C4' 'C6' 'T8' 'TP7' 'CP5' 'CP3' 'CP1' 'CPZ' 'CP2' 'CP4' 'CP6' 'TP8' 'M1' 'M2' 'P7' 'P5' 'P3' 'P1' 'PZ' 'P2' 'P4' 'P6' 'P8' 'PO7' 'PO3' 'POZ' 'PO4' 'PO8' 'O1' 'OZ' 'O2' 'CB1' 'CB2' 'AFP1' 'AFP2' 'AF7' 'AF5' 'AFZ' 'AF6' 'AF8' 'AFF5H' 'AFF3H' 'AFF1H' 'AFF2H' 'AFF4H' 'AFF6H' 'F9' 'F10' 'FFT7H' 'FFC5H' 'FFC3H' 'FFC1H' 'FFC2H' 'FFC4H' 'FFC6H' 'FFT8H' 'FT9' 'FT7' 'FT8' 'FT10' 'FTT7H' 'FCC5H' 'FCC3H' 'FCC1H' 'FCC2H' 'FCC4H' 'FCC6H' 'FTT8H' 'TTP7H' 'CCP5H' 'CCP3H' 'CCP1H' 'CCP2H' 'CCP4H' 'CCP6H' 'TTP8H' 'TPP7H' 'CPP5H' 'CPP3H' 'CPP1H' 'CPP2H' 'CPP4H' 'CPP6H' 'TPP8H' 'P9' 'P10' 'PPO3H' 'PPO1H' 'PPO2H' 'PPO4H' 'PO9' 'PO5' 'PO1' 'PO2' 'PO6' 'PO10' 'CBZ'});
pop_chanedit( chanlocs );

clear
[azimuth,elevation,r] = cart2sph(x,y,z)

%% https://stackoverflow.com/questions/32212968/how-to-project-data-onto-independent-components-using-fastica-in-matlab

% Generate data (N=1000) from distribution
X = gendata(1000);

% Estimate ICs and projections of X  
[Y, A, W] = fastica(X, 'approach', 'defl');    

% New points from the same distribution  
Xn = gendata(50);

% Project new point on ICA estimated independent components  
Yn = W*Xn;

%%
eeglab

%% test out

open fastica

x=linspace(0.1,10,2000);
y1=2*sin(x);y2=sin(2*x); % two component
figure;plot(x,y1,'r');hold on; plot(x,y2,'b');legend('y1','y2'); title('2 source');

ymix1=y1+0.1*y2;ymix2=0.5*y1+0.5*y2;ymix3=0.1*y1+y2;
figure;plot(x,ymix1,'r');hold on;plot(x,ymix2,'g');plot(x,ymix3,'b');legend('ymix1','ymix2','ymix3');title('mix signals');ylim([-2 2]);

mixedsig=[ymix1;ymix2;ymix3];

[icasig, A, W] = fastica(mixedsig); % mixing matrix A and the separating matrix W

figure; plot(x,mixedsig);legend('ymix1','ymix2','ymix3');ylim([-2 2]);

figure;plot(x,icasig(1,:));hold on; plot(x,icasig(2,:));legend('icasig1','icasig2');title('plot 2 ica signal');ylim([-2 2]);

% separating matrix W; restore ica signal
figure;plot(x,W(1,:)*mixedsig);hold on; plot(x,W(2,:)*mixedsig);
figure;plot(x,W*mixedsig);
legend('icasig1restore','icasig2restore');title('restore 2 ica signal');

% mixing matrix A; restore mix signal
figure;plot(x,A(1,:)*icasig);hold on; plot(x,A(2,:)*icasig);plot(x,A(3,:)*icasig);
figure;plot(x,A*icasig);
legend('ymix1restore','ymix2restore','ymix3restore');title('restore mix signal');

% Reconstruct the mixedsig
figure;plot(x,A(:,1)*icasig(1,:));legend('reconstruct_ymix1','reconstruct_ymix2','reconstruct_ymix3');title('keep component1'); % Keep component 1 remove 2
figure;plot(x,A(:,2)*icasig(2,:));legend('reconstruct_ymix1','reconstruct_ymix2','reconstruct_ymix3');title('keep component2'); % Keep component 2 remove 1

% % Doesn't make sense
% figure;plot(x,A(:,2)*icasig(1,:));legend('reconstruct_ymix1','reconstruct_ymix2','reconstruct_ymix3');title('keep component1'); % Keep component 1 remove 2
% figure;plot(x,A(:,1)*icasig(2,:));legend('reconstruct_ymix1','reconstruct_ymix2','reconstruct_ymix3');title('keep component2'); % Keep component 2 remove 1