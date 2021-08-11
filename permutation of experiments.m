% set the random number seed as the date of today in formate such as 20210809
seed=input('enter the date in format such as 20211209:'); %test % seed=20210809;

%***************************** for mix trials
rng(seed);
% number of trials
numtrials=50;
% number of blocks
numblocks=8;
% total trial number
numtotal=numtrials*numblocks; % should be 400
% num of conditions in the experiment
numconditions=16;
% how many semirandom permutation set in the experiment 
numPerm=numtotal/numconditions;
% create the whole set of random conditions for the experiment
allPerm=[];
for p=1:numPerm
    allPerm=[allPerm randperm(16)];
end
% figure;hist(allPerm,numconditions); %check to see if numbers of each condition is the same
figure;plot(allPerm); %check if sequence of the randomized elements is the same


%**************************** for learning trials
rng(seed);
% number of trials
numtrials=25;
% number of blocks
numblocks=16;
% total trial number
numtotal=numtrials*numblocks; % should be 400
% num of conditions in the experiment
numconditions=16;
% create the whole set of random conditions for the experiment
allPerm=[];
allPermSequence=randperm(numconditions);
for b=1:numblocks
    allPerm=[allPerm allPermSequence(b)*ones(1,numtrials)];
end
figure;hist(allPerm,numconditions); %check to see if numbers of each condition is the same
figure;plot(allPerm); %check if sequence of the randomized elements is the same

%% the set of conditions
conditionfuctions = {@A1_A2, @A1_A1, @A1_U2, @A1_U1, ...
@A2_A1, @A2_A2, @A2_U1, @A2_U2, ...
@U1_U2, @U1_U1, @U1_A2, @U1_A1, ...
@U2_U1, @U2_U2, @U2_A1, @U2_A2 };
% randomize the conditions and pick a condition
conditionPermutations = randperm(16); 
for cP = conditionPermutations
    [xL,yL,PosL,PosL0,PosL1,PosL2,PosL3,ThicknessL,ColorL,...
        xR,yR,PosR,PosR0,PosR1,PosR2,PosR3,ThicknessR,ColorR] ...
        = conditionfuctions{cP}(steplength,yCenter,screenXpixels);
end

