function [ErrorSmooth] = smoothError(Error,win)
% win is the number of samples in the smoothing window
ErrorSmooth=[];
for i=1:length(Error)-win
    ErrorSmooth(i)=mean(Error(i:i+win-1));
end
end

