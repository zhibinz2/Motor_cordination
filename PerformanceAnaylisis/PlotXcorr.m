function PlotXcorr(testdataC3,testdataC4,Fs)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    
subplot(1,5,1);
%delta 1– 3 Hz
%fdesign.bandpass((Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, Astop2)
bandpassSpecs=fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',1.5, 2, 4, 4.5, 100, 0.02, 100,Fs);
bandpassFilt = design(bandpassSpecs, 'ellip');
% fvtool(bandpassFilt);
filtered_C3=filtfilthd(bandpassFilt,testdataC3);%plot(filtered_C3);
filtered_C4=filtfilthd(bandpassFilt,testdataC4);%plot(filtered_C4);
[r,lags]=xcorr(filtered_C3, filtered_C4);
plot(lags./2,r);xlabel('time [ms]');ylabel('cross correlation');xlim([-250 250]);ylim([-2*1e8 2*1e8]);
title('delta')

subplot(1,5,2);
%theta 4 – 7 Hz
bandpassSpecs=fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',3.5, 4, 7, 7.5, 100, 0.02, 100,Fs);
bandpassFilt = design(bandpassSpecs, 'ellip');
filtered_C3=filtfilthd(bandpassFilt,testdataC3);%plot(filtered_C3);
filtered_C4=filtfilthd(bandpassFilt,testdataC4);%plot(filtered_C4);
[r,lags]=xcorr(filtered_C3, filtered_C4);
plot(lags./2,r);xlabel('time [ms]');ylabel('cross correlation');xlim([-250 250]);ylim([-2*1e8 2*1e8]);
title('theta')

subplot(1,5,3);
%alpha 8 –13 Hz
bandpassSpecs=fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',7.5, 8, 13, 13.5, 100, 0.02, 100,Fs);
bandpassFilt = design(bandpassSpecs, 'ellip');
filtered_C3=filtfilthd(bandpassFilt,testdataC3);%plot(filtered_C3);
filtered_C4=filtfilthd(bandpassFilt,testdataC4);%plot(filtered_C4);
[r,lags]=xcorr(filtered_C3, filtered_C4);
plot(lags./2,r);xlabel('time [ms]');ylabel('cross correlation');xlim([-250 250]);ylim([-2*1e8 2*1e8]);
title('alpha')

subplot(1,5,4);
%beta 14–29 Hz
bandpassSpecs=fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',13.5, 14, 29, 29.5, 100, 0.02, 100,Fs);
bandpassFilt = design(bandpassSpecs, 'ellip');
filtered_C3=filtfilthd(bandpassFilt,testdataC3);%plot(filtered_C3);
filtered_C4=filtfilthd(bandpassFilt,testdataC4);%plot(filtered_C4);
[r,lags]=xcorr(filtered_C3, filtered_C4);
plot(lags./2,r);xlabel('time [ms]');ylabel('cross correlation');xlim([-250 250]);ylim([-2*1e8 2*1e8]);
title('beta')

subplot(1,5,5);
%gamma 30 –50 Hz
d=fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',29.5, 30, 50, 55.5, 100, 0.02, 100,Fs);
bandpassFilt = design(bandpassSpecs, 'ellip');
filtered_C3=filtfilthd(bandpassFilt,testdataC3);%plot(filtered_C3);
filtered_C4=filtfilthd(bandpassFilt,testdataC4);%plot(filtered_C4);
[r,lags]=xcorr(filtered_C3, filtered_C4);
plot(lags./2,r);xlabel('time [ms]');ylabel('cross correlation');xlim([-250 250]);ylim([-2*1e8 2*1e8]);
title('gamma')

% supertitle
sgtitle('cross-correlation between C3 and C4') 
%     suptitle('I am a super title')

end