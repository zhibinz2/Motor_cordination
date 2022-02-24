TrialTime=1/srate*[1:1:size(time_series1_trialsL_sessions,1)]-10.5;
Ylim=5;
Chanselected=[1:3];

% cd D:\360MoveData\Users\alienware\Documents\GitHub\Motor_cordination\fNIR_test
cd C:\Users\zhibi\Documents\GitHub\Motor_cordination\fNIR_test\

figure;
subplot(2,2,1);
TrialsLO2ChansL=squeeze(mean(time_series1_trialsL_sessions(:,O2HbChansL(Chanselected),:),2))'; % keep all trials, but average across all right oxy channels 
TrialsLHHChansL=squeeze(mean(time_series1_trialsL_sessions(:,HHbChansL(Chanselected),:),2))';
plotSEM(TrialTime,TrialsLO2ChansL,'r');hold on;
plotSEM(TrialTime,TrialsLHHChansL,'b');hold off;grid;ylim([-Ylim Ylim]);
title('Left cortex - Left stim');

subplot(2,2,2);
TrialsLO2ChansR=squeeze(mean(time_series1_trialsL_sessions(:,O2HbChansR(Chanselected),:),2))'; % keep all trials, but average across all right oxy channels 
TrialsLHHChansR=squeeze(mean(time_series1_trialsL_sessions(:,HHbChansR(Chanselected),:),2))';
plotSEM(TrialTime,TrialsLO2ChansR,'r');hold on;
plotSEM(TrialTime,TrialsLHHChansR,'b');hold off;grid;ylim([-Ylim Ylim]);
title('Right cortex - Left stim');

subplot(2,2,3);
TrialsRO2ChansL=squeeze(mean(time_series1_trialsR_sessions(:,O2HbChansL(Chanselected),:),2))'; % keep all trials, but average across all right oxy channels 
TrialsRHHChansL=squeeze(mean(time_series1_trialsR_sessions(:,HHbChansL(Chanselected),:),2))';
plotSEM(TrialTime,TrialsRO2ChansL,'r');hold on;
plotSEM(TrialTime,TrialsRHHChansL,'b');hold off;grid;ylim([-Ylim Ylim]);
title('Left cortex - Right stim');

subplot(2,2,4);
TrialsRO2ChansR=squeeze(mean(time_series1_trialsR_sessions(:,O2HbChansR(Chanselected),:),2))'; % keep all trials, but average across all right oxy channels 
TrialsRHHChansR=squeeze(mean(time_series1_trialsR_sessions(:,HHbChansR(Chanselected),:),2))';
plotSEM(TrialTime,TrialsRO2ChansR,'r');hold on;
plotSEM(TrialTime,TrialsRHHChansR,'b');hold off;grid;ylim([-Ylim Ylim]);
title('Right cortex - Right stim');


%% plot with 95 interval

% x = 1:100;                                          % Create Independent Variable
% y = randn(50,100);                                  % Create Dependent Variable ‘Experiments’ Data
% N = size(y,1);                                      % Number of ‘Experiments’ In Data Set
% yMean = mean(y);                                    % Mean Of All Experiments At Each Value Of ‘x’
% ySEM = std(y)/sqrt(N);                              % Compute ‘Standard Error Of The Mean’ Of All Experiments At Each Value Of ‘x’
% CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
% yCI95 = bsxfun(@times, ySEM, CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’
% figure
% plot(x, yMean,'p')                                  % Plot Mean Of All Experiments
% hold on
% 
% % plot(x, yCI95+yMean,'-r')                           % Plot 95% Confidence Intervals Of All Experiments
% patch([x, fliplr(x)], [yCI95(1,:) fliplr(yCI95(2,:))], 'b', 'EdgeColor','none', 'FaceAlpha',0.25)
% hold off
% grid

%%
TrialTime=1/srate*[1:1:size(time_series1_baseline_sessions,1)]-10.5;
Ylim=5;
Chanselected=[1:11];

% cd D:\360MoveData\Users\alienware\Documents\GitHub\Motor_cordination\fNIR_test
cd C:\Users\zhibi\Documents\GitHub\Motor_cordination\fNIR_test

figure;
subplot(2,1,1);
BaselineO2ChansL=squeeze(mean(time_series1_baseline_sessions(:,O2HbChansL(Chanselected),:),2))'; % keep all trials, but average across all right oxy channels 
BaselineHHChansL=squeeze(mean(time_series1_baseline_sessions(:,HHbChansL(Chanselected),:),2))';
plotSEM(TrialTime,BaselineO2ChansL,'r');hold on;
plotSEM(TrialTime,BaselineHHChansL,'b');hold off;grid;ylim([-Ylim Ylim]);
title('Left cortex');

subplot(2,1,2);
BaselineO2ChansR=squeeze(mean(time_series1_baseline_sessions(:,O2HbChansR(Chanselected),:),2))'; % keep all trials, but average across all right oxy channels 
BaselineHHChansR=squeeze(mean(time_series1_baseline_sessions(:,HHbChansR(Chanselected),:),2))';
plotSEM(TrialTime,BaselineO2ChansR,'r');hold on;
plotSEM(TrialTime,BaselineHHChansR,'b');hold off;grid;ylim([-Ylim Ylim]);
title('Right cortex');