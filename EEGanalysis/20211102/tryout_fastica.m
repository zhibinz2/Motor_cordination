[icasig, A, W] = fastica(EEGdata);
% Each row of matrix mixedsig is one observed signal.
% the rows of icasig contain the estimated independent components.


mixedsig=EEGdata;


[icasig] = FASTICA (mixedsig); % the rows of icasig contain the
% estimated independent components.
%
[icasig, A, W] = FASTICA (mixedsig); % outputs the estimated separating
% matrix W and the corresponding mixing matrix A. 
% Each row of A and W corresponding to one channel
% Each colum of A and W correponding to the weight distribution of one componenet

[A, W] = FASTICA (mixedsig); % gives only the estimated mixing matrix
% A and the separating matrix W.

EEGdata=Out2*Out1;
EEGdataRestore=Out2*Out1;

Out1=Out3*EEGdata;
Out1Restore=Out3*EEGdata;

%%
fasticag
%%

figure;
plotx(Out1(:,10000:20000));
figure;
plotx(Out1Restore(:,10000:20000));

figure;
plotx(EEGdata(:,10000:20000));
figure;
plotx(EEGdataRestore(:,10000:20000));


imagesc(Out2);
topoplot(Out2(1,:),EEG.chanlocs(1:128),'conv','on');
topoplot(Out2(:,1),EEG.chanlocs(1:128),'conv','on');

imagesc(Out3);
topoplot(Out3(1,:),EEG.chanlocs(1:128),'conv','on');
topoplot(Out3(:,1),EEG.chanlocs(1:128),'conv','on');


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
y1=sin(x);y2=cos(2*x); % two component
figure;plot(x,y1);hold on; plot(x,y2);legend('y1','y2'); title('2 source');

ymix1=y1+y2;ymix2=0.1*y1+2*y2;ymix3=2*y1+0.1*y2;
figure;plot(x,ymix1);hold on;plot(x,ymix2);plot(x,ymix3);legend('ymix1','ymix2','ymix3');title('mix signals');

mixedsig=[ymix1;ymix2;ymix3];

[icasig, A, W] = fastica(mixedsig); % mixing matrix A and the separating matrix W

figure; plot(x,mixedsig);legend('ymix1','ymix2','ymix3');

figure;plot(x,icasig(1,:));hold on; plot(x,icasig(2,:));legend('icasig1','icasig2');title('plot 2 ica signal');

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