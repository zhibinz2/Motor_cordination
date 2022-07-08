%% https://www.mathworks.com/help/econ/varm.html
% openExample('econ/CreateVAR0ModelExample')
Mdl = varm
Mdl.AR = {NaN}

load Data_USEconModel
cpi = DataTable.CPIAUCSL;
unrate = DataTable.UNRATE;

Mdl = varm(2,4)
Mdl.Trend = NaN

c = [1; 1; 0];
Phi1 = {[0.2 -0.1 0.5; -0.4 0.2 0; -0.1 0.2 0.3]}; 
delta = [1.5; 2; 0];
Sigma = [0.1 0.01 0.3; 0.01 0.5 0; 0.3 0 1];

Mdl = varm('Constant',c,'AR',Phi1,'Trend',delta,'Covariance',Sigma)

Mdl2 = Mdl;
Phi = [zeros(3,3) Phi1];
Mdl2.AR = Phi

load Data_USEconModel

figure;
plot(DataTable.Time,DataTable.CPIAUCSL);
title('Consumer Price Index');
ylabel('Index');
xlabel('Date');

figure;
plot(DataTable.Time,DataTable.UNRATE);
title('Unemployment Rate');
ylabel('Percent');
xlabel('Date');

rcpi = price2ret(DataTable.CPIAUCSL);
unrate = DataTable.UNRATE(2:end);

Mdl = varm(2,4)

EstMdl = estimate(Mdl,[rcpi unrate])

summarize(EstMdl)

load Data_USEconModel
cpi = DataTable.CPIAUCSL;
unrate = DataTable.UNRATE;

rcpi = price2ret(cpi);
unrate = unrate(2:end);
Y = [rcpi unrate];

Mdl = varm(2,4);
EstMdl = estimate(Mdl,Y(1:(end-10),:));

YF = forecast(EstMdl,10,Y(1:(end-10),:));

figure;
plot(DataTable.Time(end - 50:end),rcpi(end - 50:end));
hold on
plot(DataTable.Time((end - 9):end),YF(:,1))
h = gca;
fill(DataTable.Time([end - 9 end end end - 9]),h.YLim([1,1,2,2]),'k',...
    'FaceAlpha',0.1,'EdgeColor','none');
legend('True CPI growth rate','Forecasted CPI growth rate',...
    'Location','NW')
title('Quarterly CPI Growth Rate: 1947 - 2009');
ylabel('CPI growth rate');
xlabel('Year');
hold off

figure;
plot(DataTable.Time(end - 50:end),unrate(end - 50:end));
hold on
plot(DataTable.Time((end - 9):end),YF(:,2))
h = gca;
fill(DataTable.Time([end - 9 end end end - 9]),h.YLim([1,1,2,2]),'k',...
    'FaceAlpha',0.1,'EdgeColor','none');
legend('True unemployment rate','Forecasted unemployment rate',...
    'Location','NW')
title('Quarterly Unemployment Rate: 1947 - 2009');
ylabel('Unemployment rate');
xlabel('Year');
hold off

%%

addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/MSE-VARFI
addpath(genpath('/home/zhibin/Documents/GitHub/granger'));
addpath(genpath('/home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/2P_testing'));
addpath /home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/DFA
addpath(genpath('/home/zhibin/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/ARFIMA'));
%d estimation (Whittle semiparametric estimator)
est_d=d_estimation(Yo(m,:),d_min,d_max);
ed(m)=est_d;
    
%%  color scheme
red   = [1 0 0];
pink  = [1 0.65 0.75];
blue  = [0 0 1];
green = [0 1 0];
darkgreen = [0 0.5 0];
grey  = [0.5 0.5 0.5];
yellow  = [1 1 0];
deepyellow  = [1 0.8 0.2];
megenta = [1 0 1];% fill([0 1 1 0],[0 0 1 1],megenta)
cyan = [0 1 1]; % fill([0 1 1 0],[0 0 1 1],cc)
purple = [0.6 0.1 0.9];
ARFIMAcolors=[deepyellow;darkgreen;pink];

%% PLOT-1 : 2 plots with error bar for 3 different length
Ns=[100 250 500]; ds=[-0.5:0.1:1]; repeat=100; stdx=20;
figure;
tic;
for i=1:3
    [Betass,BetassErr,Hss,HssErr,Destss,DestErr] = IterARFIMA(Ns(i),ds,repeat,stdx);
    subplot(1,3,1);
        hold on;
        errorbar(ds,mean(Betass,2)/2,BetassErr/2,'.','color',ARFIMAcolors(i,:));
        xlabel('d');ylabel('Beta');xlim([-0.6 1.1]);
        title(['PSA estimate of d with standard error']);
        hold off;
        if i==3;
            d_est=ds; hold on; plot(ds,d_est,'color',purple);hold off;
            legend({'N=100','N=250','N=500','B=2*d'},'Location','southeast');
        end
    subplot(1,3,2);
        hold on;
        errorbar(ds,mean(Hss,2)-0.5,HssErr,'.','color',ARFIMAcolors(i,:))
        xlabel('d');ylabel('DFA estimate of d');xlim([-0.6 1.1]);
        title(['DFA estimate of d with standard error']);
        hold off;
        if i==3;
            d_est=ds; hold on; plot(ds,d_est,'color',purple);hold off;
            legend({'N=100','N=250','N=500','H=d+0.5'},'Location','southeast');
        end
    subplot(1,3,3); % Luca's MSE_VARFI (http://www.lucafaes.net/LMSE-MSE_VARFI.html)
        hold on;
        errorbar(ds,mean(Destss,2),DestErr,'.','color',ARFIMAcolors(i,:))
        xlabel('d');ylabel('d-est');xlim([-0.6 1.1]);
        title(['d-est with standard error']);
        hold off;
        if i==3;
            d_est=ds; hold on; plot(ds,d_est,'color',purple);hold off;
            legend({'N=100','N=250','N=500','H=d+0.5'},'Location','southeast');
        end
end
suptitle(['Beta and H with different d repeating 100 times (ARFIMA([],d,[]); length of ' num2str(Ns) ')']);
toc
