function [beta] = getEEGbeta(data,maxfreq,Fs);
% This function use EEG as input and output the beta
% beta from linear fit of freq 15-50

    tempf = fft(data,[],1);
%     for m = 1:50 % 50 frequencies
%         freqs = (m-1)*10+1:m*10; % downsample freqs to only 500 samples
%         power = sum(abs(tempf(freqs,:)),1); % sum every 10 freqs
%         pow.EEG(m,:) = power;
%     end
    

    df = Fs/size(data,1);
    nbins = ceil(maxfreq/df) + 1;
    freqs = [0:(nbins-1)]*df;
    % sum all freqs in each one Hz increment
    % NincreF=floor(nbins/maxfreq); % number of bins in each increment for method 1
    pow=[];closetIndexStart=1;
    for m = 1:maxfreq % 50 frequencies
%         % method 1
%         Selectfreqs = (m-1)*NincreF+1:m*NincreF; % downsample freqs according to increment
%         power = sum(abs(tempf(Selectfreqs,:)),1); % sum every all freqs in one increment
%         pow(m) = power;
        % method 2 (more accurate)
        [minValue,closetIndexEnd]=min(abs(freqs-m)); % closetIndex in BottonPressTime
        Selectfreqs=closetIndexStart:closetIndexEnd;
        closetIndexStart=closetIndexEnd;
        power = sum(abs(tempf(Selectfreqs,:)),1); 
        pow(m) = power;
    end
    
%     plot(1:50,pow);
%     xlim([15 50]);
%     
%     plot(15:50,pow(15:50));
%     
%     loglog(15:50,pow(15:50));
    
    
%     selectInd = 1:round(length(freqs)/(maxfreq*10)):length(freqs); % downsample
%     selectInd = selectInd(ceil((15/30)*length(selectInd)):end); % select from 15 Hz - 50 Hz
%     fcoef = fft(ndetrend(data(:,1),1),[],1)/size(data,1);
%     eppow = abs(squeeze((fcoef(1:nbins,:)))).^2; 
    % figure; plot(eppow); 
    % plot(freqs,eppow,'.');xlim([15 50]);
    % plot(log10(freqs),log10(eppow),'.');xlim([log10(15) log10(50)]);
    % plot(freqs(selectInd),eppow(selectInd),'.');
    % pow = squeeze(var(fcoef(1:nbins,:),[],3)); % figure; plot(pow)
    
    xx=log10(15:50); yy=log10(pow(15:50));
%     plot(xx',yy');
    A=polyfit(xx',yy',1);
    beta=-A(1); % the slope, or beta, the first order polynomial coefficient from polyfit (same way with Hurst Componenet)
    % plot the fit
%     FitValues=polyval(A,xx');
%     hold on;
%     plot(xx',FitValues,'r--');
%     legend({'Data',['Fit (Beta=' num2str(beta) ')']},'Location','northeast');
%     hold off;
end

