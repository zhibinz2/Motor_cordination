function plot95(x,y,color);
% Plost average signal with 95 CI
%   Detailed explanation goes here
% color='b';
% x = 1:100;                                          % Create Independent Variable
% y = randn(50,100);                                  % Create Dependent Variable ‘Experiments’ Data
N = size(y,1);                                      % Number of ‘Experiments’ In Data Set
yMean = mean(y);                                    % Mean Of All Experiments At Each Value Of ‘x’
ySEM = std(y)/sqrt(N);                              % Compute ‘Standard Error Of The Mean’ Of All Experiments At Each Value Of ‘x’
CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
yCI95 = bsxfun(@times, ySEM, CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’
% figure
plot(x, yMean,color)                                  % Plot Mean Of All Experiments
hold on

% plot(x, yCI95+yMean,'-r')                           % Plot 95% Confidence Intervals Of All Experiments
patch([x, fliplr(x)], [yCI95(1,:) fliplr(yCI95(2,:))], color, 'EdgeColor','none', 'FaceAlpha',0.15)
hold off


end

