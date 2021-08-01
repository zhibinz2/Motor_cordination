%   ********* https://www.wacom.com/en-us/support/product-support/drivers

sca;
close all;
clearvars;


ifi=0.0167; % macbook
samplingRate=1/ifi; 
trialLength=5;
DrawArcSlowLoop(samplingRate,trialLength);

%% draw out the trace
load('tabletTestData.mat')
plot(xPos,yPos);

