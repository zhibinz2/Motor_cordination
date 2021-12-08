[Out1, Out2, Out3] = fastica(EEGdata);


mixedsig=EEGdata;


[icasig] = FASTICA (mixedsig); % the rows of icasig contain the
% estimated independent components.
%
[icasig, A, W] = FASTICA (mixedsig); % outputs the estimated separating
% matrix W and the corresponding mixing matrix A.

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
