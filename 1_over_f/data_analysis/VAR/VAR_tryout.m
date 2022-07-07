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