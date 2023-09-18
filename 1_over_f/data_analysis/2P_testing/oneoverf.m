function [freqs,fcoef,beta,xx,yy,FitValues] = oneoverf(y,Fs);
% this function compute do fft, plot power spectrum
%   Dand fit the 1overf 
    fcoef=fft(y);
    N=length(y);
    fcoef=fcoef/N;
    halfN=floor(N/2);
    % Fs=1/2;% assuming a sampling frequency of 0.5Hz (2 second per sample)
    freqs=linspace(0,Fs/2,halfN+1); % same df
    fcoef=2*fcoef(1:halfN);
    amplitude = abs(fcoef);
    p=amplitude.^2;
    % plot(freqs(2:end),p);
    
    % xx=log10(freqs(2:end));yy=log10(p(2:length(freqs)))';
    xx=log10(freqs(2:end));yy=log10(p(1:length(p)))';
    % xx=smoothing(log10(freqs(2:end)),20);yy=smoothing(log10(p(1:length(p))),20);
%     plot(xx, yy,'bx');xlabel('Log10(f)');ylabel('Log10(power)');
    hold on;
    
    % select freqs on power scale with base of 2 [2 4 8 16 32 ...] for the fit
    freqsind=2.^[1:floor(log2(length(yy)))];
    % open nextpow2 % % open nextpow2(length(yy))
    
    % method 1: linear fit using a table
    tbl=table(xx(freqsind)',yy(freqsind)');
    tbl=table(xx(1:10)',yy(1:10)');
    tbl=table(xx',yy');
    mdl=fitlm(tbl,'linear');
    plotAdded(mdl);
%     beta=-1*table2array(mdl.Coefficients(2,1));
    % xlabel('Log(f)');ylabel('Log(power)');title('C. Spectrum');
    %ylim([-3 3]);
    %ylim([-7 0]);
    
%     % linear fit
%     tbl=table(xx(11:end)',yy(11:end)');
%     mdl=fitlm(tbl,'linear');
%     plotAdded(mdl);

    % method 2: polyfit
%     A=polyfit(xx(freqsind)',yy(freqsind)',1);
    A=polyfit(xx',yy',1);
    beta=-A(1); % the slope, or beta, the first order polynomial coefficient from polyfit (same way with Hurst Componenet)
    % plot the fit
    FitValues=polyval(A,xx');
%     
    end

