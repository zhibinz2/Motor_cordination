function [P1,P2,Px,coh,phase,freqs] = spectra2(xx1,xx2,maxfreq,sr,win);
% It utilize xx1 and xx2 output from the function of spectra (time x chan x chunks)
% open /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing/spectra.m
% xx1 and xx2 are already in chunck of win length (time x chan x chunks)
% This function can be used for Powerspectra or Cross spectra
% /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing/P2_Syncopate_plots.m

SamplesLength=sr*win; % chunck of 10 seconds
T=SamplesLength/sr; % in second 
df=1/T;
nbins = ceil(maxfreq/df) + 1;
freqs = [0:(nbins-1)]*df;

% open /home/zhibin/Documents/GitHub/matlab/signal/allspectra.m
fcoef1 = fft(ndetrend(xx1,1),[],1)./SamplesLength;
fcoef2 = fft(ndetrend(xx2,1),[],1)./SamplesLength;

P1=mean((fcoef1.*conj(fcoef1)),2); % shouldn't we average across the 3rd dimension?
P2=mean((fcoef2.*conj(fcoef2)),2);
Px=mean((fcoef1.*conj(fcoef2)),2);
P1=P1(1:nbins);
P2=P2(1:nbins);
Px=Px(1:nbins);
coh=(abs(Px).^2)./(P1.*P2);
phase=angle(Px);

end