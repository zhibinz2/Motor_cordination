%   ********* https://www.wacom.com/en-us/support/product-support/drivers
% https://biomch-l.isbweb.org/forum/biomch-l-forums/biomch-l-1988-2010/20448-using-the-psychtoolbox-to-interface-with-a-graphics-tablet

sca;
close all;
clear;
clearvars;


try
    
ifi=0.0167; % macbook
samplingRate=1/ifi;
trialLength=5;

DrawArcSlowLoop(samplingRate,trialLength);


catch
  sca
  psychrethrow(psychlasterror);
end 

%% draw out the trace
load('tabletTestData.mat')
plot(xPos,yPos);

