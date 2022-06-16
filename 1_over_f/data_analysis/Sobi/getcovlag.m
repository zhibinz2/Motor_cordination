function [covavg,goodchans,lags] = getcovlag(thedata,samplingrate, lagstep,  maxlag,badchans,dofilter);
%
%sets up covariance by lag structure for SOBI analysis.  Optionally
%applies appropriate band pass filtering based on the lags to be examined
%thereby detrending and remcving any potential aliasing
%
%
%inputs
%thedata = (time x channels);
%samplingrate = samplingrate
%lagatep = samples between lags 
%maxlag = maximum lag (in samples) 
%badchans = bad channels
%optional inputs
%dofilter = 1 to filter data to bandwidth defined by lags, 0 to not filter
%outputs (defaults to no filters)
%covavg = lags x channels x channels
%threshhold = amplitude threshhold to apply after filtering
%goodchans = good channels
if nargin < 6 
  dofilter = 0;
end;
lags = [0:lagstep:maxlag]; 
nchans = size(thedata,2);
minfreq = samplingrate/maxlag;		
maxfreq = samplingrate/(lagstep*2);
if dofilter
[N, Wn] = buttord([minfreq maxfreq]/(samplingrate/2),[0.8*minfreq ...
		    1.2*maxfreq]/(samplingrate/2),3,6);
[b,a] = butter(N,Wn);

parfor j = nchans
  thedata(:,j) = filtfilt(b,a,thedata(:,j));
end;
end;
temp = max(abs(thedata));
temp = min(abs(thedata));
goodchans = setdiff(1:nchans,badchans);
		    
%DID THE COVARIANCE MYSELF MANUALLY SEE IF CCORR works
covavg = zeros(length(goodchans),length(goodchans),length(lags));
parfor k = 1:length(lags);
      temp1 = squeeze(thedata(1:end-lags(k),goodchans));
      temp2 = squeeze(thedata(lags(k)+1:end,goodchans));

      ip1 = temp1'*temp2; 
      covavg(:,:,k) = ip1/length(temp1);
end;