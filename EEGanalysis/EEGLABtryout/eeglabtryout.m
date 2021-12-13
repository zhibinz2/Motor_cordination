%% import array
eegdata = rand(32, 256*100);
plot(eegdata(1:2,100:200)');
plot(eegdata(1:2,:)');

%% import events
eegdata = rand(32, 256*100); % 32 channels of random activity (100 s sampled at 256 Hz).
eegdata(33,[10:256:256*100]) = 1; % simulating a stimulus onset every second
eegdata(33,[100:256:256*100]+round(rand*128)) = 2;% simulating reaction times about 500 ms after stimulus onsets
plot(eegdata(33,:),'ro');

eegdata(33,:)=[];

%% ica
runica(eegdata);

round(EEG.etc.ic_classification.ICLabel.classifications(1:10,:)*100)

%  There are six categories of components Brain, Muscle, Eye, Heart, Line Noise, Channel Noise, and Other.