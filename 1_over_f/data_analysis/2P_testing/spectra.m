function [P1,P2,Px,coh,phase,freqs,xx1,xx2] = spectra(x1,x2,maxfreq,sr,win);
% x1 and x2 should be of the same length in the first dimension
% This function can be used for Powerspectra or Cross spectra
% It reorgainze the timeseries into chunks and then fft
% It was used in
% /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing/P2_Syncopate_plots.m

if win=[]; win=10; end % chunck of 10 seconds
SamplesLength=sr*win; % number of samples in the window

% maximum number of 10 second windows (# of chuncks)
nTrials=min([floor(length(x1)/SamplesLength) floor(length(x2)/SamplesLength)]);

% zcore transform
x1=zscore(x1);x2=zscore(x2);

% reorganize the time series into trials (time x chan x chunks)
xx1=[];xx2=[];
for i=1:nTrials
    xx1(:,i)=x1((SamplesLength*(i-1)+1):(SamplesLength*i),:);
    xx2(:,i)=x2((SamplesLength*(i-1)+1):(SamplesLength*i),:);
end

T=SamplesLength/sr; % in second 
df=1/T;
nbins = ceil(maxfreq/df) + 1;
freqs = [0:(nbins-1)]*df;

fcoef1 = fft(ndetrend(xx1,1),[],1)./SamplesLength;
fcoef2 = fft(ndetrend(xx2,1),[],1)./SamplesLength;

P1=mean((fcoef1.*conj(fcoef1)),2);
P2=mean((fcoef2.*conj(fcoef2)),2);
Px=mean((fcoef1.*conj(fcoef2)),2);
P1=P1(1:nbins);
P2=P2(1:nbins);
Px=Px(1:nbins);
coh=(abs(Px).^2)./(P1.*P2);
phase=angle(Px);

end

