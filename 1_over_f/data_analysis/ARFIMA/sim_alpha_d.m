freqs=[0.1:0.1:50];
figure;
alpha=0.6;
P=1./(freqs.^alpha);
plot(freqs,P,'g');% 1overf
xlabel('freqs');ylabel('1/f');
xlim([0 10]);ylim([0 10])
hold on;

alpha=0.8;
P=1./(freqs.^alpha);
plot(freqs,P,'b');% 1overf
xlabel('freqs');ylabel('1/f');
xlim([0 10]);ylim([0 10])

alpha=1;
P=1./(freqs.^alpha);
plot(freqs,P,'r');% 1overf
xlabel('freqs');ylabel('1/f');
xlim([0 4]);ylim([0 4])

legend('alpha=0.6','alpha=0.8','alpha=1');

%% fit the slope to find d