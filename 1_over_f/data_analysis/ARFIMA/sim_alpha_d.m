%% examine beta or alpha or the slope
freqs=[0.1:0.1:50];

figure;
alpha1=0.5;
P=1./(freqs.^alpha1);
plot(freqs,P,'g');% 1overf
xlabel('freqs');ylabel('1/f');
% xlim([0 10]);ylim([0 10])
hold on; 

alpha2=1;
P=1./(freqs.^alpha2);
plot(freqs,P,'b');% 1overf
xlabel('freqs');ylabel('1/f');
% xlim([0 10]);ylim([0 10])

alpha3=2;
P=1./(freqs.^alpha3);
plot(freqs,P,'r');% 1overf
xlabel('freqs');ylabel('1/f');
% xlim([0 4]);ylim([0 4])

alpha4=3;
P=1./(freqs.^alpha4);
plot(freqs,P,'m');% 1overf
xlabel('freqs');ylabel('1/f');
xlim([0 4]);ylim([0 4])

legend(['alpha1= ' num2str(alpha1)], ['alpha2= ' num2str(alpha2)],...
    ['alpha3= ' num2str(alpha3)], ['alpha4= ' num2str(alpha4)]);

%% fit the slope to find d
addpath(genpath('/home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/ARFIMA'));
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/DFA
%%% ARFIMA(0,d,0)
figure('units','normalized','outerposition',[0 0 1 0.5]);
N=600;d=1;stdx=20;
[Z] = ARFIMA_SIM(N,[],[],d,stdx);
y=Z+750;
subplot(1,4,1); plot(y,'b');
subplot(1,4,2); autocorr(y,N-1);
subplot(1,4,3); Fs=1.3; [freqs,fcoef] = oneoverf(y,Fs);xlabel('Log(f)');ylabel('Log(power)');
subplot(1,4,4); [D,Alpha1]=DFA_main(y);


addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/ARFIMA/dgp_arfima
% ARFIMA(0,0.8,0) >>>>  Results=dgp_arfima(0,[],[],600,1,0.8);
figure('units','normalized','outerposition',[0 0 1 0.5]);
T=600;SD=1;d=0.5;F=0;
Z=dgp_arfima(0,[],[],T,SD,d,F);
y=Z+750;
subplot(1,4,1); plot(y,'b');
subplot(1,4,2); autocorr(y,T-1);
subplot(1,4,3); Fs=1.3; [freqs,fcoef] = oneoverf(y,Fs);xlabel('Log(f)');ylabel('Log(power)');
subplot(1,4,4); [D,Alpha1]=DFA_main(y);
