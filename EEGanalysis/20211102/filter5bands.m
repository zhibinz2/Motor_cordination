detrend_data=double(detrend_data);

eeglab;
close all;

%delta 1– 3 Hz
%fdesign.bandpass((Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, Astop2)
bandpassSpecs=fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',1.5, 2, 4, 4.5, 100, 0.02, 100,Fs);
bandpassFilt = design(bandpassSpecs, 'ellip');
% fvtool(bandpassFilt);
filtered_delta=filtfilthd(bandpassFilt,detrend_data);

%theta 4 – 7 Hz
bandpassSpecs=fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',3.5, 4, 7, 7.5, 100, 0.02, 100,Fs);
bandpassFilt = design(bandpassSpecs, 'ellip');
filtered_theta=filtfilthd(bandpassFilt,detrend_data);

%alpha 8 –13 Hz
bandpassSpecs=fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',7.5, 8, 13, 13.5, 100, 0.02, 100,Fs);
bandpassFilt = design(bandpassSpecs, 'ellip');
filtered_alpha=filtfilthd(bandpassFilt,detrend_data);


%beta 14–29 Hz
bandpassSpecs=fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',13.5, 14, 29, 29.5, 100, 0.02, 100,Fs);
bandpassFilt = design(bandpassSpecs, 'ellip');
filtered_beta=filtfilthd(bandpassFilt,detrend_data);


%gamma 30 –50 Hz
bandpassSpecs=fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',29.5, 30, 50, 50.5, 100, 0.02, 100,Fs);
bandpassFilt = design(bandpassSpecs, 'ellip');
filtered_gamma=filtfilthd(bandpassFilt,detrend_data);

%%

subplot(1,5,1);
plot(filtered_delta(10000:12000,28));% try out
title('delta');

subplot(1,5,2);
plot(filtered_theta(10000:12000,28));
title('theta');

subplot(1,5,3);
plot(filtered_alpha(10000:12000,28));
title('alpha');

subplot(1,5,4);
plot(filtered_beta(10000:12000,28));
title('beta');

subplot(1,5,5);
plot(filtered_gamma(10000:12000,28));
title('gamma');

suptitle('5 freq band');

% plot(detrend_data);

