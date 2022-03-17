% Set screen flash 100HZ
% 0.01s per flash, 0.02 per frame (flip every other flash / every 20 ms)
ifi;

% Mean interval
MeanTapInterval=0.5; % second , or 500 ms
NumFramesInterval=MeanTapInterval/ifi*2; % let's say NumFramesInterval=25 frames


% 3 trials/conditions 5 blocks
% estimated about 2 min (240 taps) per trial/condition, total 30 min of tapping

%% Total time Condition 1:
TotalTimeCondition1=0.5*240 % in seconds

Showframes=[1:NumFramesInterval:NumFramesInterval*numTaps] 
% total time
Showframes(end)*20/1000/60

%% Total time Condition 2:
% Same as condition 1 except the instruction

%% Total time Condition 3:
% TotalTimeCondition2=0.5*240;
RandomIntervals = round(normrnd(NumFramesInterval,NumFramesInterval/5,[1,numTaps]))
mean(RandomIntervals) % mean interval 
plot(RandomIntervals,'ro');
% in case some interval are negative
if any(RandomIntervals<0)
    RandomIntervals(find(RandomIntervals<0))=1;
end
RandomIntervals
Showframes=cumsum(RandomIntervals)
Showframes(end)
sum(RandomIntervals)*20/1000/60 % total time
Showframes(end)*20/1000/60 % total time



%%

for n=1:max(Showframes)
    if any(Showframes(:) == n);
        Screen('flip'...

