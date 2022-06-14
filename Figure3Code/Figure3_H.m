%% Plot the percent of spines added across 11 days
% Used for Figure 3H

% Written by Nora Wolcott 05/21/21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load data
addition = importdata('Figure3H_PercentSpinesAdded.mat');
addition(~addition) = NaN; % convert 0s to NaNs

%% Calculate average and standard error
add_avg = mean(addition,2,'omitnan')';
add_sterror = std(addition,0,2,'omitnan')/sqrt(length(addition))';

%% Plot percent addition
figure
x = 1:10;
plot(x,addition,'o', 'MarkerEdgeColor', [132, 169, 172]/255,'MarkerFaceColor', [132, 169, 172]/255)
hold on
errorbar(x,add_avg,add_sterror,'Color',[83, 53, 74]/255,'LineWidth', 4)
ylim([0 60])
hold on
xlabel('Day #')
ylabel('% added')
title('Percent spines added')
hold off