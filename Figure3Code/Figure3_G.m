%% Plot the spine survival fraction
% Used for Figure 3G

% Written by Nora Wolcott 08/03/21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load Data
survival = importdata('Figure3G_SpineSurvivalFraction.mat');
survival(~survival) = NaN; % convert 0s to NaNs

%% Calculate average and standard error
sur_avg = mean(survival,2,'omitnan')';
sur_sterror = std(survival,0,2,'omitnan')/sqrt(length(survival))';

%% Plot survival fraction
figure
x = 1:11;
plot(x,survival, 'o', 'MarkerEdgeColor', [132, 169, 172]/255,'MarkerFaceColor',[132, 169, 172]/255)
hold on
errorbar(x,sur_avg,sur_sterror,'Color',[83, 53, 74]/255,'LineWidth',4)
ylim([0 100])
hold on
xlabel('Day #')
ylabel('Percent spines survived')
title('Survival Fraction')
hold off