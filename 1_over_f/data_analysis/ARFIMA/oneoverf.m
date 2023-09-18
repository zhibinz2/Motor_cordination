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
    h=plot(xx, yy,'x');xlabel('Log10(f)');ylabel('Log10(power)');
    % the following line skip the name of the previous plot from the legend
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hold on;
    

    A=polyfit(xx,yy,1);
    beta=A(1); % the slope, or beta, the first order polynomial coefficient from polyfit (same way with Hurst Componenet)
    % plot the fit
    FitValues=polyval(A,xx);
    plot(xx,FitValues,'r--');
    legend({['Beta = ' num2str(round(-1*beta,2))]},'Location','southwest');
    end

