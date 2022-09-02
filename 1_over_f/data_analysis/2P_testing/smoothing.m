function [smoothx] = smoothing(x,win)
% x could be time series of errors or button presses intervals in colums
% win is the number of samples in the smoothing window
smoothx=[];
for i=1:size(x,1)-win
    smoothx(i,:,:)=median(x(i:i+win-1,:,:));
end
end

