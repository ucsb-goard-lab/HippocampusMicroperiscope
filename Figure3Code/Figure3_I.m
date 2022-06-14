%% Plot the percent of spines subtracted across 11 days
% Used for Figure 3I

% Written by Nora Wolcott 05/21/21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load data
subtraction = importdata('Figure3I_PercentSpinesSubtracted.mat');
subtraction(~subtraction) = NaN;

%% Calculate average and standard error
sub_avg = mean(subtraction,2,'omitnan')';
sub_sterror = std(subtraction,0,2,'omitnan')/sqrt(length(subtraction))';

%% Plot percent addition
figure
x = 1:10;
plot(x,subtraction, 'o', 'MarkerEdgeColor', [132, 169, 172]/255,'MarkerFaceColor', [132, 169, 172]/255)
hold on
errorbar(x,sub_avg,sub_sterror,'Color',[83, 53, 74]/255,'LineWidth', 4)
ylim([0 60])
hold on
xlabel('Day #')
ylabel('% subtracted')
title('Percent spines subtracted')
hold off