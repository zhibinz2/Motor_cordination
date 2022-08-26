function [delta, theta, alpha, beta, gamma] = sum5band(EEG,win)
%This funciton get the valuses of power sum from the 5 EEG band width
%   Detailed explanation goes here

%%
% Input: x is EEG in time x chan
% Output: xx is EEG in time x chan x chunks

%%

df=1/win;
nbins = ceil(maxfreq/df) + 1;
freqs=(0:1:nbins-1)*df;
fcoef = fft(xx,[],1);
pow= squeeze(var(fcoef(1:nbins,:,:),[],3));
eppow = abs(squeeze(mean(fcoef(1:nbins,:,:),3)));


end