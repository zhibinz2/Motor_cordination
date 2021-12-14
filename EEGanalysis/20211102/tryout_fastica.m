cd /home/zhibin/Documents/Acquisition/bimanual_Reach_Jack_2021111802/
EEGfileName='bimanual_Reach_Jack_20211118'; EEG=loadcurry([pwd '/' EEGfileName '.cdt']);
cd /home/zhibin/Documents/GitHub/Motor_cordination/EEGanalysis/20211102
tic
run organize_EEG_filter_step2.m


% shorten EEG to try out
PlotEnd=350; % how many trial to plot
EEGdataShort=filtered_data(1:EventInd(PlotEnd),:);
TimesdataShort=datatimes(1:EventInd(PlotEnd));
% 
plot(TimesdataShort, EEGdataShort);hold on;
for i=1:PlotEnd
    xline(EventInd(i));
end
hold off;

tic
[icasig, A, W] = fastica(EEGdataShort');
toc
% Each row of matrix mixedsig is one observed signal.
% the rows of icasig contain the estimated independent components.

% other EEGLAB method (tryout)
tic
[weights,sphere,compvars,bias,signs,lrates,activations] = runica(EEGdataShort'); % train using defaults 
toc


%% tryout
imagesc(sphere)
plot(activations');
plot(icasig);

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

%% tryout
fasticag
%% tryout

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
topoplot(A(:,1),test,'conv','on','nosedir','+Y');
%% tryout
for i=1:12
    subplot(1,12,i);
    topoplot(A(:,i),test,'nosedir','+Y'); % plot each colum
end

for i=1:12
    subplot(1,12,i);
    topoplot(A(i,:),test,'nosedir','+Y');
end
%% tryout
imagesc(W);
topoplot(W(1,:),EEG.chanlocs(1:128),'conv','on');
topoplot(W(:,1),EEG.chanlocs(1:128),'conv','on');

imagesc(A);
topoplot(A(1,:),EEG.chanlocs(1:128),'conv','on');
topoplot(A(:,1),EEG.chanlocs(1:128),'conv','on');

Chanlocs_topolot=EEG.chanlocs(1:128);

test=EEG.chanlocs;

topoplot(A(1,:),Neuroscan_topoplot_chans);

topoplot(A(1,:),chanlocs(1:128),'conv','on');

%% Plot ICA component

for i=1:128
    SqA(i)=sumsqr(A(:,i));
end
figure;
plot(1:128,SqA,'ro');ylabel('sum of square of column in A');xlabel('ICs');
[B,I]=sort(SqA,'descend');

% Component 56 55 3 have the highest weights in the mixing matrics A
% topoplot to examine them
cd /home/zhibin/Documents/GitHub/Motor_cordination/EEGanalysis
load('Neuroscan_spherical__topoplot_chans.mat')
figure;
subplot(3,1,1);topoplot(A(:,56),test,'nosedir','+Y');title('component 56');colorbar;
subplot(3,1,2);topoplot(A(:,52),test,'nosedir','+Y');title('component 52');colorbar;
subplot(3,1,3);topoplot(A(:,3),test,'nosedir','+Y');title('component 3'); colorbar;

figure;
subplot(3,1,1);topoplot(A(:,102),test,'nosedir','+Y');title('component 102');colorbar;
subplot(3,1,2);topoplot(A(:,105),test,'nosedir','+Y');title('component 105');colorbar;
subplot(3,1,3);topoplot(A(:,116),test,'nosedir','+Y');title('component 116'); colorbar;

%% Display signal before ICA
PlotStart=95;PlotEnd=100;DisplayChan=[3 55 56];
figure;
subplot(2,1,2);
plotx(TimesdataShort(EventInd(PlotStart):EventInd(PlotEnd)),EEGdataShort((EventInd(PlotStart):EventInd(PlotEnd)),DisplayChan));
hold on;
for i=PlotStart:PlotEnd
    xline(EventInd(i),'r',{'End of trial'});
    xline(EventInd(i)-1000,'g',{'Go signal'});
end
hold off;
xlim([EventInd(PlotStart+1)-2000 EventInd(PlotEnd)-2000]);
legend({'3','55','56'});
title('Channel Singal Before ICA');xlabel('time samples');ylabel('uV');
subplot(2,1,1);
plotx(TimesdataShort(EventInd(PlotStart):EventInd(PlotEnd)),EEGdataShort((EventInd(PlotStart):EventInd(PlotEnd)),:));
hold on;
for i=PlotStart:PlotEnd
    xline(EventInd(i),'r',{'End of trial'});
    xline(EventInd(i)-1000,'g',{'Go signal'});
end
hold off;
xlim([EventInd(PlotStart+1)-2000 EventInd(PlotEnd)-2000]);
title('Channel Singal Before ICA');

%% Plot ICA signal
PlotStart=95;PlotEnd=100;DisplayChan=[3 55 56];
figure;
subplot(2,1,2);
plotx(TimesdataShort(EventInd(PlotStart):EventInd(PlotEnd)),icasig(DisplayChan,(EventInd(PlotStart):EventInd(PlotEnd))));
hold on;
for i=PlotStart:PlotEnd
    xline(EventInd(i),'r',{'End of trial'});
    xline(EventInd(i)-1000,'g',{'Go signal'});
end
hold off;
xlim([EventInd(PlotStart+1)-2000 EventInd(PlotEnd)-2000]);
legend({'3','55','56'});
title('IC time series');xlabel('time samples');ylabel('uV');
subplot(2,1,1);
plotx(TimesdataShort(EventInd(PlotStart):EventInd(PlotEnd)),icasig(:,(EventInd(PlotStart):EventInd(PlotEnd))));
hold on;
for i=PlotStart:PlotEnd
    xline(EventInd(i),'r',{'End of trial'});
    xline(EventInd(i)-1000,'g',{'Go signal'});
end
hold off;
xlim([EventInd(PlotStart+1)-2000 EventInd(PlotEnd)-2000]);
title('IC time series');


%% Calculate Correlation
% test
% [RHO3,PVAL3] = corr(EEGdataShort(:,3),EEGdataShort(:,3));

% FP1 and FP2 are channel 1 and 3;
% compute correlation between FP1 and FP2;

[RHO1,PVAL1] = corr(EEGdataShort(:,1),icasig');
figure
subplot(1,2,1);
% yyaxis left
plot(1:128,RHO1,'bo');hold on;ylabel('correlation coefficient');
% yyaxis right
plot(1:128,PVAL1,'ro');ylabel('p value');xlabel('ICs');
legend({'correlation coefficient' 'p-values '});title('correlation with FP1');hold off;
subplot(1,2,2);
[RHO3,PVAL3] = corr(EEGdataShort(:,3),icasig');
% yyaxis left
plot(1:128,RHO3,'bo');hold on;ylabel('correlation coefficient');
% yyaxis right
plot(1:128,PVAL3,'ro');ylabel('p value');xlabel('ICs');
legend({'correlation coefficient' 'p-values '});title('correlation with FP2');hold off;

% component 55,56,3 seem to be highly correlated with FP1 and FP2
% Examine them
load('Neuroscan_spherical__topoplot_chans.mat')
topoplot(A(:,55),test,'nosedir','+Y');
topoplot(A(:,56),test,'nosedir','+Y');
topoplot(A(:,3),test,'nosedir','+Y');



%% Display signal (tryout)

PlotEnd=5;
figure
plotx(TimesdataShort(1:EventInd(PlotEnd)),EEGdataShort((1:EventInd(PlotEnd)),:));

hold on;
for i=1:PlotEnd
    xline(EventInd(i));
end
hold off;

figure
plotx(TimesdataShort(1:EventInd(PlotEnd)),icasig(:,(1:EventInd(PlotEnd))));

%% Remove components and display
% A(:,55)=0;A(:,56)=0;A(:,3)=0;
% icasig(65,:)=0;icasig(52,:)=0;icasig(3,:)=0;

% test
% OO=ones(128,128);OO(DisplayChan,:)=0;imagesc(OO);

A(:,DisplayChan)=0; icasig(DisplayChan,:)=0;

mixedsig=A*icasig;

figure
% Plot ICA signal
PlotStart=95;PlotEnd=100;DisplayChan=[3 55 56];
figure;
subplot(2,1,2);
plotx(TimesdataShort(EventInd(PlotStart):EventInd(PlotEnd)),mixedsig(DisplayChan,(EventInd(PlotStart):EventInd(PlotEnd))));
hold on;
for i=PlotStart:PlotEnd
    xline(EventInd(i),'r',{'End of trial'});
    xline(EventInd(i)-1000,'g',{'Go signal'});
end
hold off;
xlim([EventInd(PlotStart+1)-2000 EventInd(PlotEnd)-2000]);
legend({'3','55','56'});
title('Channel 3 55 56 with ICs removed');xlabel('time samples');ylabel('uV');
subplot(2,1,1);
plotx(TimesdataShort(EventInd(PlotStart):EventInd(PlotEnd)),mixedsig(:,(EventInd(PlotStart):EventInd(PlotEnd))));
hold on;
for i=PlotStart:PlotEnd
    xline(EventInd(i),'r',{'End of trial'});
    xline(EventInd(i)-1000,'g',{'Go signal'});
end
hold off;
xlim([EventInd(PlotStart+1)-2000 EventInd(PlotEnd)-2000]);
title('Mixed Signal with ICs removed');



















%% tryout


for i=1:128
    SqA(i)=sumsqr(W(:,i));
end
plot(1:128,SqA,'ro');ylabel('sum of square of column in A');

for i=1:128
    SqA(i)=sumsqr(A(i,:));
end
plot(1:128,SqA,'ro');ylabel('sum of square of rows in A');

RemoveICA1sig=A(:,2:128)*icasig(2:128,:);%remove #1 component and get back 128 channel signal

[icasig2, A2, W2] = fastica(RemoveICA1sig);

for i=1:12
    subplot(1,12,i);
    topoplot(A2(:,i),test,'nosedir','+Y'); % plot each colum
end

plotx(icasig2(1:127,1000:5000)');



%% plot channels in 3d
% 139 channel labels
chanlocs = struct('labels', {'FP1' 'FPZ' 'FP2' 'AF3' 'AF4' 'F11' 'F7' 'F5' 'F3' 'F1' 'FZ' 'F2' 'F4' 'F6' 'F8' 'F12' 'FT11' 'FC5' 'FC3' 'FC1' 'FCZ' 'FC2' 'FC4' 'FC6' 'FT12' 'T7' 'C5' 'C3' 'C1' 'CZ' 'C2' 'C4' 'C6' 'T8' 'TP7' 'CP5' 'CP3' 'CP1' 'CPZ' 'CP2' 'CP4' 'CP6' 'TP8' 'M1' 'M2' 'P7' 'P5' 'P3' 'P1' 'PZ' 'P2' 'P4' 'P6' 'P8' 'PO7' 'PO3' 'POZ' 'PO4' 'PO8' 'O1' 'OZ' 'O2' 'CB1' 'CB2' 'AFP1' 'AFP2' 'AF7' 'AF5' 'AFZ' 'AF6' 'AF8' 'AFF5H' 'AFF3H' 'AFF1H' 'AFF2H' 'AFF4H' 'AFF6H' 'F9' 'F10' 'FFT7H' 'FFC5H' 'FFC3H' 'FFC1H' 'FFC2H' 'FFC4H' 'FFC6H' 'FFT8H' 'FT9' 'FT7' 'FT8' 'FT10' 'FTT7H' 'FCC5H' 'FCC3H' 'FCC1H' 'FCC2H' 'FCC4H' 'FCC6H' 'FTT8H' 'TTP7H' 'CCP5H' 'CCP3H' 'CCP1H' 'CCP2H' 'CCP4H' 'CCP6H' 'TTP8H' 'TPP7H' 'CPP5H' 'CPP3H' 'CPP1H' 'CPP2H' 'CPP4H' 'CPP6H' 'TPP8H' 'P9' 'P10' 'PPO3H' 'PPO1H' 'PPO2H' 'PPO4H' 'PO9' 'PO5' 'PO1' 'PO2' 'PO6' 'PO10' 'CBZ' 'VEOG' 'HEOG' 'EMG1' 'EMG2' 'HL 1' 'HL 2' 'EMG3' 'EMG4' 'EMG5' 'EMG6' 'TRIGGER'});
% 128 channel labels
chanlocs = struct('labels', {'FP1' 'FPZ' 'FP2' 'AF3' 'AF4' 'F11' 'F7' 'F5' 'F3' 'F1' 'FZ' 'F2' 'F4' 'F6' 'F8' 'F12' 'FT11' 'FC5' 'FC3' 'FC1' 'FCZ' 'FC2' 'FC4' 'FC6' 'FT12' 'T7' 'C5' 'C3' 'C1' 'CZ' 'C2' 'C4' 'C6' 'T8' 'TP7' 'CP5' 'CP3' 'CP1' 'CPZ' 'CP2' 'CP4' 'CP6' 'TP8' 'M1' 'M2' 'P7' 'P5' 'P3' 'P1' 'PZ' 'P2' 'P4' 'P6' 'P8' 'PO7' 'PO3' 'POZ' 'PO4' 'PO8' 'O1' 'OZ' 'O2' 'CB1' 'CB2' 'AFP1' 'AFP2' 'AF7' 'AF5' 'AFZ' 'AF6' 'AF8' 'AFF5H' 'AFF3H' 'AFF1H' 'AFF2H' 'AFF4H' 'AFF6H' 'F9' 'F10' 'FFT7H' 'FFC5H' 'FFC3H' 'FFC1H' 'FFC2H' 'FFC4H' 'FFC6H' 'FFT8H' 'FT9' 'FT7' 'FT8' 'FT10' 'FTT7H' 'FCC5H' 'FCC3H' 'FCC1H' 'FCC2H' 'FCC4H' 'FCC6H' 'FTT8H' 'TTP7H' 'CCP5H' 'CCP3H' 'CCP1H' 'CCP2H' 'CCP4H' 'CCP6H' 'TTP8H' 'TPP7H' 'CPP5H' 'CPP3H' 'CPP1H' 'CPP2H' 'CPP4H' 'CPP6H' 'TPP8H' 'P9' 'P10' 'PPO3H' 'PPO1H' 'PPO2H' 'PPO4H' 'PO9' 'PO5' 'PO1' 'PO2' 'PO6' 'PO10' 'CBZ'});
pop_chanedit( chanlocs );

clear
[azimuth,elevation,r] = cart2sph(x,y,z)

%% tryout
% https://stackoverflow.com/questions/32212968/how-to-project-data-onto-independent-components-using-fastica-in-matlab

% Generate data (N=1000) from distribution
X = gendata(1000);

% Estimate ICs and projections of X  
[Y, A, W] = fastica(X, 'approach', 'defl');    

% New points from the same distribution  
Xn = gendata(50);

% Project new point on ICA estimated independent components  
Yn = W*Xn;

%% tryout
eeglab

%% test out with my own signal combination

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
