rng default

Fs = 1000;
t = 0:1/Fs:1-1/Fs;

% ｘ cause y
x = cos(2*pi*100*t) + sin(2*pi*200*t) + 0.5*randn(size(t));
y = 2.3*cos(2*pi*100*t) + 1.2*sin(2*pi*200*t)+ 0.5*randn(size(t));
[Cxy1,F1] = mscohere(x,y,hamming(100),80,100,Fs);
[Cxy2,F2] = mscohere(x(201:1000),y(1:800),hamming(100),80,100,Fs);
figure
subplot(121)
plot(1:1000,x);hold on;plot(1:1000,y);
legend({'x','y'});
subplot(122)
plot(F1,Cxy1);hold on;plot(F2,Cxy2);
title('Magnitude-Squared Coherence')
xlabel('Frequency (Hz)')
grid
sgtitle('x cause y')
legend({'x','y'});
% second case: one of them are periodic signal,shift by n periods
x = cos(2*pi*100*t) + sin(2*pi*200*t) + 0.5*randn(size(t));
y = cos(2*pi*100*t - pi/4) + 0.35*sin(2*pi*100*t - pi/2); 
figure
subplot(121)
plot(1:1000,x);hold on;plot(1:1000,y);
legend({'x','y'});
subplot(122)
plot(F1,Cxy1);hold on;plot(F2,Cxy2);
title('Magnitude-Squared Coherence')
legend({'x','y'});
xlabel('Frequency (Hz)')
grid
sgtitle('y is periodic, shifted by 20 periods')

% if we used cross-spectral density normalized by the power spectral
% densities as coherence as we did with EEG
csd1 =cov(fft([x' y'])/length(t));
coh1 = (abs(csd1).^2) ./ (diag(csd1) * diag(csd1)'); % coherence before shift
csd2 =cov(fft([x(201:1000)' y(1:800)'])/length(t));
coh2 = (abs(csd2).^2) ./ (diag(csd2) * diag(csd2)'); % coherence after shift
% alternative method
sf1=corrcoef(fft([x' y'])/length(t));
coh1=abs(sf1).^2;% coherence after shift
sf2=corrcoef(fft([x(201:1000)' y(1:800)'])/length(t));
coh2=abs(sf2).^2;% coherence after shift