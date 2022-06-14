function [] = Figure6_CA1DGAxisStats()
%-------------------------------------------------------------------------%
%   This function analyzes computes the F statistic for the place field
%   width/spatial info as a function of distance along the CA1-DG axis. 
%
%   Last updated by WTR 05/31/2022
%-------------------------------------------------------------------------%
%% Load data
clear
close all

load('CA1DG_neuron_pos.mat')
load('CA1DG_DFF_place_cells.mat'); 

% Place field width (uncomment if you want to run PF width analysis)
% pos = neuron_pos(place_cell_vec == 1); 
% y = FWHM_vec(place_cell_vec == 1);

% Spatial information (uncomment if you want to run spatial info analysis)
pos = neuron_pos; 
y = spatial_info;

%% Remove any cells with NaN spatial info, clip cells at ends if desired (end cells can cause overfitting)
lower_bound_dist = -600;
upper_bound_dist = 600;
exclude_idx1 = isnan(y);
exclude_idx2 = pos<lower_bound_dist;
exclude_idx3 = pos>upper_bound_dist;
include_idx = ~or(or(exclude_idx1,exclude_idx2),exclude_idx3);
pos = pos(include_idx);
y = y(include_idx);

%% Fit average data with spline (full model)
smooth_param = 1e-8;    % Use lower values for fewer params: 1e-8 was used for spatial info analysis and 5e-8 was used for place field width analysis
[curve, ~, output] = fit(pos',y','smoothingspline','SmoothingParam',smooth_param);
num_param = ceil(output.numparam);

%% Calculate sum of squares error for full model (spline) fit
y_model = curve(pos);
full_model_sse = sum((y'-y_model).^2);

%% Calculate sum of squares error for reduced model
y_model = mean(y)*ones(length(pos),1);
reduc_model_sse = sum((y'-y_model).^2);

%% Plot
subplot(1,2,1)
plot(curve, pos',y'); hold on 
xlabel('Position on DG-to-CA1 axis')
ylabel('Spatial Info (bit/event)')
title(['Smoothing spline model of data, number of parameters = ' ...
    num2str(num_param) ', sse = ' num2str(full_model_sse)])

subplot(1,2,2)
scatter(pos',y','b.')
hold on
plot(pos, mean(y)*ones(1,length(pos)),'r')
xlabel('Position on DG-to-CA1 axis')
ylabel('Spatial Info (bit/event)')
title(['Null hypothesis fit, number of parameters = 1, sse = ' num2str(reduc_model_sse)])

set(gcf,'Position',[300 450 1200 500])
set(gcf,'color',[1 1 1])
pause
close

%% General Linear F-test
dF_full = length(pos)-num_param; % degrees of freedom for full model (n - params)
dF_reduc = length(pos)-1;        % degrees of freedom for reduced model (n - 1)

F_stat = ((reduc_model_sse - full_model_sse) / (dF_reduc- dF_full)) / (full_model_sse / dF_full);
dF_numerator = dF_reduc-dF_full;
df_denominator = dF_full;
p_value = 1-fcdf(F_stat,dF_numerator,df_denominator);
disp(['F(' num2str(dF_numerator) ',' num2str(df_denominator) ') = ' num2str(F_stat) ', p = ' num2str(p_value)])

