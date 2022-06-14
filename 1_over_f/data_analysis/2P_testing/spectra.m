function [fcoef1,fcoef2,cprod,freqs] = spectra(x1,x2,maxfreq,sr);
% x1 and x2 should be of the same length in the first dimension
% This function can be used for Powerspectra or Cross spectra

% cut x1 and x2 into 20s per trials
nTrials=min([floor(length(x1)/(sr*20)) floor(length(x2)/(sr*20))]);
xx1=[];xx2=[];
for i=1:nTrials
    xx1(:,i)=x1((sr*20*(i-1)+1):(sr*20*i),:);
    xx2(:,i)=x2((sr*20*(i-1)+1):(sr*20*i),:);;
end

SamplesLength=sr*20;
T=SamplesLength/sr; % in second 
df=1/T;
nbins = ceil(maxfreq/df) + 1;
freqs = [0:(nbins-1)]*df;

fcoef1 = mean(fft(ndetrend(xx1,1),[],1)/SamplesLength,2);
fcoef2 = mean(fft(ndetrend(xx2,1),[],1)/SamplesLength,2);

cprod=[];
for i=1:nbins
    cprod(i,:)=fcoef1(i,:).*conj(fcoef2(i,:))/SamplesLength;
end

end

