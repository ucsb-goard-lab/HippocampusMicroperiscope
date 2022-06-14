%% Create violin plots of average spine density and turnover by class
% Calculate ANOVA for data

% Used for Figure 3E and 3F

% Written by Nora Wolcott 07/27/21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load data
spine_density = importdata('Figure3E_SpineDensity.mat');
spine_turnover = importdata('Figure3F_SpineTurnover.mat');

%% Calculate ANOVA
group = {'Stubby','Thin','Mushroom','Filopodia'};
[p,tbl1,stats] = anova1(spine_density,group,'off');
[p2,tbl2,stats2] = anova1(spine_turnover,group,'off');
% save('spineDensity_ANOVA.mat','tbl1')
% save('spineTurnover_ANOVA.mat','tbl2')

%% Plot 
% Plot spine density per 10um (Figure 3E)
figure
violinplot(spine_density)
title('Average # spines per 10um')
xticklabels({'Stubby','Thin','Mushroom','Filopodia'})
ylim([0 3.5])

% Plot % spine turnover (Figure 3F)
figure
violinplot(spine_turnover)
title('Percent spine turnover')
xticklabels({'Stubby','Thin','Mushroom','Filopodia'})
ylim([0 20])