%% anova
% Example data for three groups with unequal sizes
group1 = [5, 6, 7, 8, 9];          % Group 1 data
group2 = [4, 5, 6, 7, 8, 9, 10];   % Group 2 data (larger size)
group3 = [3, 4, 1, 2];             % Group 3 data (smaller size)

% Combine the data into a single vector
data = [group1, group2, group3];

% Create a grouping variable
group = [ones(size(group1)), 2 * ones(size(group2)), 3 * ones(size(group3))];

% Perform one-way ANOVA
[p, tbl, stats] = anova1(data, group, 'off');

% Display the ANOVA table
disp('ANOVA Table:');
disp(tbl);

% Check for significant differences using multcompare
c = multcompare(stats, 'Display', 'off');

% Extract significant comparisons
significant_comparisons = c(c(:, 6) < 0.05, :);

% Calculate means and standard deviations for each group
means = [mean(group1), mean(group2), mean(group3)];
stds = [std(group1), std(group2), std(group3)];

% Plot the bar plot with error bars
figure;
bar(means, 'FaceColor', [0.7 0.7 0.7]);
hold on;
errorbar(1:3, means, stds, 'k.', 'LineWidth', 1.5);
xlabel('Groups');
ylabel('Mean Value');
title('Bar Plot with Error Bars and Significant Comparisons');
xticks(1:3);
xticklabels({'Group 1', 'Group 2', 'Group 3'});

% Mark significant comparisons with asterisks
for i = 1:size(significant_comparisons, 1)
    group1_idx = significant_comparisons(i, 1);
    group2_idx = significant_comparisons(i, 2);
    y_max = max(means(group1_idx), means(group2_idx)) + stds(group1_idx);
    line([group1_idx, group2_idx], [y_max, y_max], 'Color', 'k', 'LineWidth', 1.5);
    text(mean([group1_idx, group2_idx]), y_max + 0.1, '*', 'HorizontalAlignment', 'center', 'FontSize', 14);
end

hold off;

%% t test
% Example data for two groups with unequal sizes
group1 = [5, 6, 7, 8, 9];          % Group 1 data
group2 = [4, 5, 6, 7, 8, 9, 10];   % Group 2 data (larger size)

% Perform a two-sample t-test (assuming unequal variances)
[h, p, ci, stats] = ttest2(group1, group2, 'Vartype', 'unequal');

% Display the results
disp('Two-Sample t-Test Results:');
disp(['h = ', num2str(h)]);        % h = 1 indicates rejection of the null hypothesis
disp(['p = ', num2str(p)]);        % p-value
disp(['t-statistic = ', num2str(stats.tstat)]); % t-statistic
disp(['Degrees of freedom = ', num2str(stats.df)]); % Degrees of freedom
disp(['Confidence interval = [', num2str(ci(1)), ', ', num2str(ci(2)), ']']);

% Interpret the results
if h == 1
    disp('There is a significant difference between the two groups (reject null hypothesis).');
else
    disp('There is no significant difference between the two groups (fail to reject null hypothesis).');
end