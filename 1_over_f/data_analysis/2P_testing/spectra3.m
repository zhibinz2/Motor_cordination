function [xx] = spectra3(x,sr,win)
% Reorganize EEG into time x chan x chunks
% Then fft and get power freq x chan
% Use a window of 10s with df 0.1Hz and sum up each one Hz range

%%
% Input: x is EEG (time x chan)
%%
if win==[]; win=10; end % chunck of 10 seconds
SamplesLength=sr*win; % number of samples in the window

% maximum number of 10 second windows (# of chuncks)
% nTrials=min([floor(length(x1)/SamplesLength) floor(length(x2)/SamplesLength)]);
nTrials=floor(size(x,1)/SamplesLength);

% zcore transform
% x1=zscore(x1);x2=zscore(x2);
x=zscore(x);

% reorganize the time series into trials (time x chunks)
% xx1=[];xx2=[];
xx=[];
for i=1:nTrials
%     xx1(:,i)=x1((SamplesLength*(i-1)+1):(SamplesLength*i),:);
%     xx2(:,i)=x2((SamplesLength*(i-1)+1):(SamplesLength*i),:);
    xx(:,:,nTrials)=x((SamplesLength*(i-1)+1):(SamplesLength*i),:);
end

end

