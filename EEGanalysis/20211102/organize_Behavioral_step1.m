%% load behaviral data
% clear;
% cd /home/zhibin/Documents/GitHub/Motor_cordination/behaviraldata
% load matlab20211031.mat;

%% Reorgainize behavioral data into BehavData and condition
behaviraldata.subjectnumber=seed;
subjectnumber=behaviraldata.subjectnumber;
dataBlock=behaviraldata.dataBlock;

% Just to try out with one trial and plot
dataTrialL=dataBlock.dataTrialL;  
dataTrialR=dataBlock.dataTrialR;  
dataTrialJ=dataBlock.dataTrialJ;

% Plot one trial
% xLyL=dataTrialL.xLyL; plot(xLyL(:,1),xLyL(:,2),'r'); hold on;
% xRyR=dataTrialR.xRyR; plot(xRyR(:,1),xRyR(:,2),'b');
% xJyJ=dataTrialJ.xJyJ; plot(xJyJ(:,1),xJyJ(:,2),'g'); hold off;
% set(gca, 'YDir', 'reverse');xlim([0,screenXpixels]);ylim([0,screenYpixels]);


blockNumber=dataBlock.blockNumber;
xLyL=dataTrialL.xLyL;
condition=dataTrialL.condition;

numBlock=length(dataBlock);
BehavData={};
CondiData=[]; 

for b=1:numBlock
    for t=1:length(behaviraldata.dataBlock(b).dataTrialL)
        dataxLyLtemp=behaviraldata.dataBlock(b).dataTrialL(t).xLyL;
        dataxRyRtemp=behaviraldata.dataBlock(b).dataTrialR(t).xRyR;
        dataxJyJtemp=behaviraldata.dataBlock(b).dataTrialJ(t).xJyJ;
        dataxLyLxRyRxJyJtemp=[dataxLyLtemp dataxRyRtemp dataxJyJtemp];
        dataConditiontemp=behaviraldata.dataBlock(b).dataTrialL(t).condition;
        BehavData=[BehavData; {dataxLyLxRyRxJyJtemp}];
        CondiData=[CondiData; dataConditiontemp];% same as allPerm
    end
end

% ans=BehavData{1};
% CondiData

NumTrialsRecorded=length(CondiData);

