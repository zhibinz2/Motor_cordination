addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/MSE-VARFI
addpath(genpath('/home/zhibin/Documents/GitHub/granger'));
addpath(genpath('/home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing'));
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/DFA
addpath(genpath('/home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/ARFIMA'));

%%

Ns=[100 250 500]; ds=[-0.5:0.1:1]; repeat=100; stdx=20;
i=2; d=0.7;
N=Ns(i);
% [Z] = ARFIMA_SIM(N,[],[],d,stdx); % Good H for shorter length
Z=dgp_arfima(0,[],[],N,stdx,d,0); % Good H for all length

d_min=-0.5; d_max=1;
%d estimation (Whittle semiparametric estimator)
est_d=d_estimation(Z,d_min,d_max);

%d filter
[Ytmp,~]=remove_d(Z',est_d); %filtered data

%d estimation again (Whittle semiparametric estimator)
est_d2=d_estimation(Ytmp,d_min,d_max);

figure('units','normalized','outerposition',[0 0 1 0.3]);
subplot(2,1,1);plot(Z);title('simulated time series: d=0.7');
subplot(2,1,2);plot(Ytmp);title(['removed est-d of ' num2str(est_d) '; Now d= ' num2str(est_d2)]);

