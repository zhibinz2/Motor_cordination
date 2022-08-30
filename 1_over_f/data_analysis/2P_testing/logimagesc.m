function [logpow] = logimagesc(pow,freqs)
%UNTITLED2 Summary of this function goes here
%   pow=y1(i).pow; freqs;
% % my band class
% delta 1-3, theta 4-7, alpha 8-12, beta 13-30, gamma 30-50
deltainds=find(freqs>0 & freqs <3.1);
thetainds=find(freqs>3 & freqs <7.9);
alphainds=find(freqs>7.9 & freqs <12.1);
betainds=find(freqs>12.1 & freqs <30.1);
gammainds=find(freqs>30.1 & freqs <50);
delta = sum(pow(deltainds,:),1); % sum_of_pow x chan
theta = sum(pow(thetainds,:),1); 
alpha = sum(pow(alphainds,:),1); 
beta  = sum(pow(betainds,:),1); 
gamma = sum(pow(gammainds,:),1); 
logpow = [gamma; beta; alpha; theta; delta];
imagesc (logpow);colorbar;colormap('jet');
yticks([1:5]);yticklabels({'\gamma','\beta','\alpha','\theta','\delta'});
end