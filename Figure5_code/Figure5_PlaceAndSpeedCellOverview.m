%% PlaceAndSpeedCellOverview
%-------------------------------------------------------------------------%
% This script generates all the figures for the place and speed cell 
% analysis for each of the hippocampal subregions (Figs. 5C-H). 
%
% Requires the package Violin Plot
% (https://www.mathworks.com/matlabcentral/fileexchange/45134-violin-plot).
%
% Written by WTR 04/29/2022 // Last updated by WTR 04/29/2022
%-------------------------------------------------------------------------%
%% Loading data
CA1_LapActivity = load('CA1_lap_by_lap_activity_DFF.mat'); 
CA1_PCs = load('CA1_place_cells.mat'); 
CA1_SCs = load('CA1_speed_cells.mat'); 
CA3_LapActivity = load('CA3_lap_by_lap_activity_DFF.mat'); 
CA3_PCs = load('CA3_place_cells.mat'); 
CA3_SCs = load('CA3_speed_cells.mat');
DG_LapActivity = load('DG_lap_by_lap_activity_DFF.mat'); 
DG_PCs = load('DG_place_cells.mat'); 
DG_SCs = load('DG_speed_cells.mat');

load('all_regions_place_field_width.mat'); 
load('all_regions_spatial_info_place_cells.mat'); 

%% Example place cells
nPos = size(CA1_LapActivity.L, 2);
trackLength = 19.5*pi; 
cmPerBin = trackLength / nPos; 

CA1pcExs = [108, 123, 139, 141];

figure
for ii = 1:4
    A = CA1_LapActivity.L(:, :, CA1pcExs(ii)); 
    laps = find(nansum(A, 2) ~= 0); 
    A = A(laps, :);  
    meanResp = nanmean(A) * 100; 
    semResp = nanstd(A) ./ sqrt(sum(~isnan(A))) * 100;
    
    subplot(2, 2, ii)
    x = [1:nPos fliplr(1:nPos)]*cmPerBin;
    y = [meanResp + semResp, fliplr(meanResp - semResp)];
    fill(x,y,[0 167/255 157/255],'edgecolor',[0 167/255 157/255])
    hold on
    plot([1:nPos]*cmPerBin,meanResp,'k','linewidth',3)
    xlabel('Position (cm)')
    ylabel('DF/F')
    title(['CA1 neuron #' num2str(CA1pcExs(ii)) ' place field'])
    set(gcf,'color',[1 1 1])
end

CA3pcExs = [84, 60, 24, 27];

figure
for ii = 1:4
    A = CA3_LapActivity.L(:, :, CA3pcExs(ii)); 
    laps = find(nansum(A, 2) ~= 0); 
    A = A(laps, :);  
    meanResp = nanmean(A) * 100; 
    semResp = nanstd(A) ./ sqrt(sum(~isnan(A))) * 100;
    
    subplot(2, 2, ii)
    x = [1:nPos fliplr(1:nPos)]*cmPerBin;
    y = [meanResp + semResp, fliplr(meanResp - semResp)];
    fill(x,y,[239/255 65/255 54/255],'edgecolor',[239/255 65/255 54/255])
    hold on
    plot([1:nPos]*cmPerBin,meanResp,'k','linewidth',3)
    xlabel('Position (cm)')
    ylabel('DF/F')
    title(['CA3 neuron #' num2str(CA1pcExs(ii)) ' place field'])
    set(gcf,'color',[1 1 1])
end

DGpcExs = [11, 83, 111, 31];

figure
for ii = 1:4
    A = DG_LapActivity.L(:, :, DGpcExs(ii)); 
    laps = find(nansum(A, 2) ~= 0); 
    A = A(laps, :);  
    meanResp = nanmean(A) * 100; 
    semResp = nanstd(A) ./ sqrt(sum(~isnan(A))) * 100;
    
    subplot(2, 2, ii)
    x = [1:nPos fliplr(1:nPos)]*cmPerBin;
    y = [meanResp + semResp, fliplr(meanResp - semResp)];
    fill(x,y,[127/255 63/255 152/255],'edgecolor',[127/255 63/255 152/255])
    hold on
    plot([1:nPos]*cmPerBin,meanResp,'k','linewidth',3)
    xlabel('Position (cm)')
    ylabel('DF/F')
    title(['DG neuron #' num2str(CA1pcExs(ii)) ' place field'])
    set(gcf,'color',[1 1 1])
end
   
%% Example speed cells
CA1scExs = [185, 202]; 

figure
for ii = 1:2
    meanResp = CA1_SCs.activity_binned(CA1scExs(ii), :); 
    semResp = CA1_SCs.sde_activity_binned(CA1scExs(ii), :);

    subplot(1, 2, ii)
    x = [CA1_SCs.speed_binned(CA1scExs(ii), :), fliplr(CA1_SCs.speed_binned(CA1scExs(ii), :))]; 
    y = [meanResp + semResp, fliplr(meanResp - semResp)];
    
    fill(x,y,[0 167/255 157/255],'edgecolor',[0 167/255 157/255])
    hold on
    plot(CA1_SCs.speed_binned(CA1scExs(ii), :),meanResp,'k','linewidth',3)
    xlabel('Speed (mm/s)')
    ylabel('DF/F')
    title(['CA1 neuron #' num2str(CA1scExs(ii)) ' speed tuning'])
    set(gcf,'color',[1 1 1])
end

CA3scExs = [37, 56]; 

figure
for ii = 1:2
    meanResp = CA3_SCs.activity_binned(CA3scExs(ii), :); 
    semResp = CA3_SCs.sde_activity_binned(CA3scExs(ii), :);

    subplot(1, 2, ii)
    x = [CA3_SCs.speed_binned(CA3scExs(ii), :), fliplr(CA3_SCs.speed_binned(CA3scExs(ii), :))]; 
    y = [meanResp + semResp, fliplr(meanResp - semResp)];
    
    fill(x,y,[239/255 65/255 54/255],'edgecolor',[239/255 65/255 54/255])
    hold on
    plot(CA3_SCs.speed_binned(CA3scExs(ii), :),meanResp,'k','linewidth',3)
    xlabel('Speed (mm/s)')
    ylabel('DF/F')
    title(['CA3 neuron #' num2str(CA3scExs(ii)) ' speed tuning'])
    set(gcf,'color',[1 1 1])
end
    
DGscExs = [102, 239]; 

figure
for ii = 1:2
    meanResp = DG_SCs.activity_binned(DGscExs(ii), :); 
    semResp = DG_SCs.sde_activity_binned(DGscExs(ii), :);

    subplot(1, 2, ii)
    x = [DG_SCs.speed_binned(DGscExs(ii), :), fliplr(DG_SCs.speed_binned(DGscExs(ii), :))]; 
    y = [meanResp + semResp, fliplr(meanResp - semResp)];
    
    fill(x,y,[127/255 63/255 152/255],'edgecolor',[127/255 63/255 152/255])
    hold on
    plot(DG_SCs.speed_binned(DGscExs(ii), :),meanResp,'k','linewidth',3)
    xlabel('Speed (mm/s)')
    ylabel('DF/F')
    title(['DG neuron #' num2str(DGscExs(ii)) ' speed tuning'])
    set(gcf,'color',[1 1 1])
end

%% Cross-validated responses 
figure

CA1pcsIdx = find(CA1_PCs.place_cell_vec == 1); 
CA1CrossValResp = nan(length(CA1pcsIdx), nPos);
CA1PrefPos = zeros(1, length(CA1pcsIdx)); 

for ii = 1:length(CA1pcsIdx)
    A = CA1_LapActivity.L(:, :, CA1pcsIdx(ii)); 
    laps = find(nansum(A, 2) ~= 0); 
    A = A(laps, :); 
    nLaps = length(laps); 
    halfLaps = round(nLaps/2);
    randvec = randperm(nLaps);
    half1_idx = randvec(1:halfLaps);
    half2_idx = randvec(halfLaps+1:end);
    [~, CA1PrefPos(ii)] = max(nanmean(A(half1_idx, :))); 
    CA1CrossValResp(ii, :) = nanmean(A(half2_idx, :)); 
    CA1CrossValResp(ii, :) = (CA1CrossValResp(ii, :) - min(CA1CrossValResp(ii, :))) / (max(CA1CrossValResp(ii, :)) - min(CA1CrossValResp(ii, :))); 
    CA1CrossValResp(ii, :) = smooth(CA1CrossValResp(ii, :)); 
end

[~,sort_idx] = sort(CA1PrefPos,'ascend');
CA1SortedCrossValResp = CA1CrossValResp(sort_idx,:);

subplot(1, 3, 1)
imagesc(CA1SortedCrossValResp)
title('CA1')
xlabel('Position');
ylabel('Place cells')

CA3pcsIdx = find(CA3_PCs.place_cell_vec == 1); 
CA3CrossValResp = nan(length(CA3pcsIdx), size(CA3_LapActivity.L, 2));
CA3PrefPos = zeros(1, length(CA3pcsIdx)); 

for ii = 1:length(CA3pcsIdx)
    A = CA3_LapActivity.L(:, :, CA3pcsIdx(ii)); 
    laps = find(nansum(A, 2) ~= 0); 
    A = A(laps, :); 
    nLaps = length(laps); 
    halfLaps = round(nLaps/2);
    randvec = randperm(nLaps);
    half1_idx = randvec(1:halfLaps);
    half2_idx = randvec(halfLaps+1:end);
    [~, CA3PrefPos(ii)] = max(nanmean(A(half1_idx, :))); 
    CA3CrossValResp(ii, :) = nanmean(A(half2_idx, :)); 
    CA3CrossValResp(ii, :) = (CA3CrossValResp(ii, :) - min(CA3CrossValResp(ii, :))) / (max(CA3CrossValResp(ii, :)) - min(CA3CrossValResp(ii, :))); 
    CA3CrossValResp(ii, :) = smooth(CA3CrossValResp(ii, :)); 
end

[~,sort_idx] = sort(CA3PrefPos,'ascend');
CA3SortedCrossValResp = CA3CrossValResp(sort_idx,:);

subplot(1, 3, 2)
imagesc(CA3SortedCrossValResp)
title('CA3')
   
DGpcsIdx = find(DG_PCs.place_cell_vec == 1); 
DGCrossValResp = nan(length(DGpcsIdx), size(DG_LapActivity.L, 2));
DGPrefPos = zeros(1, length(DGpcsIdx)); 

for ii = 1:length(DGpcsIdx)
    A = DG_LapActivity.L(:, :, DGpcsIdx(ii)); 
    laps = find(nansum(A, 2) ~= 0); 
    A = A(laps, :); 
    nLaps = length(laps); 
    halfLaps = round(nLaps/2);
    randvec = randperm(nLaps);
    half1_idx = randvec(1:halfLaps);
    half2_idx = randvec(halfLaps+1:end);
    [~, DGPrefPos(ii)] = max(nanmean(A(half1_idx, :))); 
    DGCrossValResp(ii, :) = nanmean(A(half2_idx, :)); 
    DGCrossValResp(ii, :) = (DGCrossValResp(ii, :) - min(DGCrossValResp(ii, :))) / (max(DGCrossValResp(ii, :)) - min(DGCrossValResp(ii, :))); 
    DGCrossValResp(ii, :) = smooth(DGCrossValResp(ii, :)); 
end

[~,sort_idx] = sort(DGPrefPos,'ascend');
DGSortedCrossValResp = DGCrossValResp(sort_idx,:);

subplot(1, 3, 3)
imagesc(DGSortedCrossValResp)
title('DG')

%% Place field width and spatial info distributions
figure
subplot(1, 2, 1)
violinplot(F, {'CA1', 'CA3', 'DG'}); 
title('Place field width'); 
ylabel('Place field width (cm)'); 

subplot(1, 2, 2)
violinplot(S, {'CA1', 'CA3', 'DG'}); 
title('Spatial info'); 
ylabel('Spatial info (bits/inf.spike)')

%% Pie chart functional type breakdown
figure

CA1conjunct = sum(CA1_SCs.speed_cell_vec .* CA1_PCs.place_cell_vec); 
CA1pcs = sum(CA1_PCs.place_cell_vec) - CA1conjunct; 
CA1scs = sum(CA1_SCs.speed_cell_vec) - CA1conjunct; 
CA1none = length(CA1_PCs.place_cell_vec) - (CA1pcs + CA1scs + CA1conjunct); 
subplot(1, 3, 1)
pie([CA1pcs, CA1scs, CA1conjunct, CA1none]); 
title('CA1') 
legend('place cells', 'speed cells', 'conjunctive', 'neither'); 

CA3conjunct = sum(CA3_SCs.speed_cell_vec .* CA3_PCs.place_cell_vec); 
CA3pcs = sum(CA3_PCs.place_cell_vec) - CA3conjunct; 
CA3scs = sum(CA3_SCs.speed_cell_vec) - CA3conjunct; 
CA3none = length(CA3_PCs.place_cell_vec) - (CA3pcs + CA3scs + CA3conjunct); 
subplot(1, 3, 2)
pie([CA3pcs, CA3scs, CA3conjunct, CA3none]); 
title('CA3') 

DGconjunct = sum(DG_SCs.speed_cell_vec .* DG_PCs.place_cell_vec); 
DGpcs = sum(DG_PCs.place_cell_vec) - DGconjunct; 
DGscs = sum(DG_SCs.speed_cell_vec) - DGconjunct; 
DGnone = length(DG_PCs.place_cell_vec) - (DGpcs + DGscs + DGconjunct); 
subplot(1, 3, 3)
pie([DGpcs, DGscs, DGconjunct, DGnone]); 
title('DG') 



