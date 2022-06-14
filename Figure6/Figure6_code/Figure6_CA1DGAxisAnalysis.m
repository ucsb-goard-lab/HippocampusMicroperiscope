function [] = Figure6_CA1DGAxisAnalysis()
%-------------------------------------------------------------------------%
%   This function analyzes place cell width and spatial information as a
%   function of position along the CA1-DG axis. 
%
%   Last updated by WTR 05/31/2022
%-------------------------------------------------------------------------%
%% Load data
% clear; close all

load('CA1DG_neuron_pos.mat')
load('CA1DG_DFF_place_cells.mat'); 

% Place field width (uncomment if you want to run PF width analysis)
% pos = neuron_pos(place_cell_vec == 1); 
% y = FWHM_vec(place_cell_vec == 1);

% Spatial information (uncomment if you want to run spatial info analysis)
pos = neuron_pos; 
y = spatial_info;

%% Apply a sliding window to spatial information/place field width 
% Parameters
win_size = 150;
sigma = 75;
pos_min = -500;                 % minimum position to analyze (in microns)
pos_max = 500;                  % maximum position to analyze (in microns)
test_pos = pos_min:pos_max;
iterations = 1000;              % iterations for random shuffling
info_by_pos = zeros(1,length([pos_min:pos_max]));
info_by_pos_sem = nan(1,length([pos_min:pos_max]));
info_by_pos_rand = nan(iterations,length([pos_min:pos_max]));

% Make Gaussian window
X = [-win_size:win_size];
G = normpdf(X,0,sigma);
G = G/max(G);

% Determine spatial_info for each point (weighted by distance from center)
for i = pos_min:pos_max
    y_curr = 0;
    weight_curr = 0;
    i_y = []; 
    for j = 1:length(pos)
        if (pos(j) >= i-win_size+1) &&  (pos(j) <= i+win_size-1)
            if ~isnan(y(j))
                pos_diff = round(i-pos(j)+win_size);
                weight = G(pos_diff);            
                y_curr = y_curr+y(j)*weight;
                i_y = [i_y, y(j)];
                weight_curr = weight_curr+weight;
            end
        end
    end
    
    info_by_pos(i-pos_min+1) = y_curr/weight_curr;
    info_by_pos_sem(i-pos_min+1) = std(i_y)/sqrt(length(i_y));

end

% Randomize and determine spatial info/place field width for each point 
for iter = 1:iterations
    rand_vec = randperm(length(pos));
    pos_rand = pos(rand_vec);
    for i = pos_min:pos_max
        y_curr = 0;
        weight_curr = 0;
        for j = 1:length(pos)
            if (pos_rand(j) >= i-win_size+1) &&  (pos_rand(j) <= i+win_size-1)
                if ~isnan(y(j))
                    pos_diff = round(i-pos_rand(j)+win_size);
                    weight = G(pos_diff);
                    y_curr = y_curr+y(j)*weight;
                    weight_curr = weight_curr+weight;
                end
            end
        end
        info_by_pos_rand(iter,i-pos_min+1) = y_curr/weight_curr;
    end
end
info_by_pos_rand_mean = nanmean(info_by_pos_rand); % mean of random spatial info
info_by_pos_rand_sem = nanstd(info_by_pos_rand);   % bootstrap SEM

%% Comparing the real values to the shuffled distribution
alpha = 0.05;         % alpha value to assess signifigance to shuffled distribution
k_min = 3;            % minimum number of consecutive positions that are significant 
sig_vec = zeros(1, length(test_pos)); 

for i = 1:length(test_pos)
    if info_by_pos_rand_mean(i) < info_by_pos(i)
        if sum(info_by_pos(i) < info_by_pos_rand(:, i)) < (iterations * alpha) 
            sig_vec(i) = 1;
        end
    else
        if sum(info_by_pos(i) > info_by_pos_rand(:, i)) < (iterations * alpha) 
            sig_vec(i) = 1;
        end
    end
end

%% Plotting
figure
rand_x = [test_pos fliplr(test_pos)];
rand_y = [info_by_pos_rand_mean+info_by_pos_rand_sem fliplr(info_by_pos_rand_mean-info_by_pos_rand_sem)];
fill(rand_x,rand_y,[0.85 0.85 0.85])
hold on
plot(test_pos,info_by_pos_rand_mean,'color',[0.5 0.5 0.5],'linewidth',2)

data_x = [test_pos fliplr(test_pos)];
data_y = [info_by_pos+info_by_pos_sem fliplr(info_by_pos-info_by_pos_sem)];
fill(data_x,data_y,[0.5 0.75 1])
plot(test_pos,info_by_pos,'color',[0 0.25 0.75],'linewidth',2)

sig_bar_y = max(info_by_pos+info_by_pos_sem);
CC = bwconncomp(sig_vec); 
for i = 1:CC.NumObjects
    if length(CC.PixelIdxList{i}) >= k_min
        plot(test_pos(CC.PixelIdxList{i}), sig_bar_y * ones(1, length(CC.PixelIdxList{i})), 'g-');
    end
end

xlabel('Position along CA1-to-DG Axis')
ylabel('Place field width'); 
set(gcf,'color',[1 1 1])
set(gcf,'position',[450 400 750 550])
shg
