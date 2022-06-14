%% PlaceCellClassifier
%-------------------------------------------------------------------------%
%   This script classifies cells as place cells, given that the following
%   criteria are met:
%       1) The cell's activity is reliable across laps
%       2) The cell's activity is well described by a Gaussian 
%
%   Written by MJG 06/18/2021 // Last updated by WTR 04/22/2022
%-------------------------------------------------------------------------%
%% Load data
clear; close all
load DG_lap_by_lap_activity_DFF.mat

%% User parameters
save_file_name = 'DG_DFF_place_cells';
alpha_val = 0.01;             % alpha value for reliability testing
iterations = 500;             % number of iterations for reliability testing
track_length = 19.5*pi;       % calculate form path radius
FWHM_max = track_length / 2;  % maximum place field width
FWHM_min = 2.50;              % minimum place field width 
gof_thresh = 0.375;           % minimum r^2 value on gaussian
amp_thresh = 0.50;            % minmum ratio of Gaussian amplitude to baseline
max_dist2center = 5;          % maximum distance a trial peak can be from the overall peak
min_peak_ratio = 0.50;        % minimum percent of overall peak a lap peak has to have
cohen_thresh = 0.5;           % Cohen's D threshold

save_flag = 1;                % save data?
plot_fits = 0;                % plot all reliable cells and fits
plot_resp = 0;                % plot responses of significant cells
shuffle_data = 0;             % shuffle actual data - only for testing! (default = 0)
 
c = [127/255 63/255 152/255]; % colors for the place cells lap-by-lap activity        
cmap = cat(1, ones(10, 3), ...
    [linspace(1, c(1), 100); ...
    linspace(1, c(2), 100); linspace(1, c(3), 100)]');

%% Calculated parameters
lap_activity = L;
numSamp = size(lap_activity,2);
numCell = size(lap_activity,3);
cm_per_bin = track_length/72;

%% Initialize vectors
reliable_cell_vec = zeros(1,numCell);
nan_cell_vec = zeros(1,numCell);
baseline_vec = zeros(1, numCell); 
amp_vec = zeros(1,numCell);
FWHM_vec = zeros(1,numCell);
gof_vec = zeros(1,numCell);
place_cell_vec = zeros(1,numCell);
preferred_pos = zeros(1,numCell);
cross_validated_resp = zeros(numCell,numSamp); 

fit_func = fittype('a0+a1*exp(-((x-b1)/c1)^2)'); % define function: gaussian w/ baseline
initial_param = [0 0 36 20];                     % initial parameters: [a0 a1 b1 c1]

%% Test lap-wise reliability of each cell
for i = 1:numCell
    i
    % initialize vectors
    bt_CC_data = zeros(1,iterations);
    bt_CC_rand = zeros(1,iterations);
    
    for j = 1:iterations
        
        % Remove empty trials (and keep track of trials with NaN values)
        full_trials = nanmean(lap_activity(:,:,i),2);
        numLaps = sum(isfinite(full_trials));
        halfLaps = round(numLaps/2);
        finite_trials = mean(lap_activity(:,:,i),2);
        
        % Extract neuron data and initialize random matrix
        activity_data = lap_activity(1:numLaps,:,i);
        activity_rand = zeros(size(activity_data));
               
        for k = 1:numLaps           
            % interpolate any NaN values
            if isnan(finite_trials(k))
                trial_interp = activity_data(k,:);
                nanx = isnan(trial_interp);   
                t = 1:numel(trial_interp);
                trial_interp(nanx) = interp1(t(~nanx), trial_interp(~nanx), t(nanx));
                if ~isempty(intersect(find(nanx == 1), 1))
                    trial_interp(1) = (trial_interp(2) + trial_interp(end))/2; 
                elseif ~isempty(intersect(find(nanx == 1), numSamp))
                    trial_interp(end) = (trial_interp(1) + trial_interp(end - 1))/2; 
                end
                
                activity_data(k,:) = trial_interp;
            end
            
       
            % circularly shuffle random matrices
            activity_rand(k,:) = circshift(activity_data(k,:),randi(numSamp));
            if shuffle_data == 1                  % shuffle actual data if desired (only for testing)
                activity_data(k,:) = circshift(activity_data(k,:),randi(numSamp));
            end            
        end
               
        % Randomly split laps into two halves
        randvec = randperm(numLaps);
        half1_idx = randvec(1:halfLaps);
        half2_idx = randvec(halfLaps+1:end);
        
        % Calculate CC from two halves of data
        lap_half1 = activity_data(half1_idx,:);
        lap_half2 = activity_data(half2_idx,:);

        CC = corrcoef(nanmean(lap_half1),nanmean(lap_half2));
        bt_CC_data(j) = CC(2);
    
        % Calculate CC from two halves of randomized data
        lap_half1 = activity_rand(half1_idx,:);
        lap_half2 = activity_rand(half2_idx,:);
        CC = corrcoef(nanmean(lap_half1),nanmean(lap_half2));
        bt_CC_rand(j) = CC(2);
        
    end

    % Test actual CC distribution against shuffled distribution
    try
        P = ranksum(bt_CC_data, bt_CC_rand, 'tail', 'right');
        nan_cell_vec(i) = 0;
        
        x1 = bt_CC_data; 
        x2 = bt_CC_rand; 
        n1 = numel(x1);
        n2 = numel(x2);
        mean_x1 = nanmean(x1);
        mean_x2 = nanmean(x2);
        var_x1  = nanvar(x1);
        var_x2  = nanvar(x2);
        meanDiff = (mean_x1 - mean_x2);
        sv1 = ((n1-1)*var_x1);
        sv2 = ((n2-1)*var_x2);
        numer =  sv1 + sv2;
        denom = (n1 + n2 - 2);
        pooledSD =  sqrt(numer / denom); % pooled Standard Deviation
        s = pooledSD;             % re-name
        d =  meanDiff / s;        % Cohen's d (for independent samples)
        
    catch
        disp(['Error: Unable to calculate reliability for neuron #' num2str(i)])
        nan_cell_vec(i) = 1;
        P = 1;
    end
    
    if P < alpha_val && d > cohen_thresh
        reliable_cell_vec(i) = 1;
    end
    
    %% Test Gaussian fit for each cell
    % Extract mean response
    mean_resp = nanmean(activity_data)*100;
    
    % Shift max response to center
    [~,max_idx] = max(smooth(mean_resp));
    centered_resp = circshift(mean_resp,round(numSamp/2)-max_idx);
    
    % Fit gaussian and save full width at half-maximum and goodness-of-fit
    fit_func = fittype('a0+a1*exp(-((x-b1)/c1)^2)');        % define function: gaussian w/ baseline
    initial_param = [0 0 36 20];                            % initial parameters: [a0 a1 b1 c1]
    [fitobj,gof] = fit([1:numSamp]',centered_resp',fit_func,'StartPoint',initial_param);
    baseline_vec(i) = fitobj.a0; 
    amp_vec(i) = fitobj.a1;
    FWHM_vec(i) = 2*sqrt(log(2)) * fitobj.c1 * cm_per_bin;
    gof_vec(i) = gof.rsquare;
    
    % Plot fit if desired
    if plot_fits==1 && reliable_cell_vec(i)==1
        figure(1)
        plot(fitobj,[1:numSamp],centered_resp);
        title(['Neuron #' num2str(i) ', FWHM = ' num2str(FWHM_vec(i),'%.1f')...
            ' cm, r^2 = ' num2str(gof_vec(i),'%.2f')])
        figure(2)
        imagesc(circshift(activity_data, [1, round(numSamp/2)-max_idx]))
        colormap(cmap)
        pause
    end
    
    % Determine if current cell is place cell (then plot response if desired)
    if reliable_cell_vec(i)==1 && amp_vec(i)>0 && FWHM_vec(i)<FWHM_max && FWHM_vec(i)>FWHM_min && gof_vec(i)>gof_thresh && (amp_vec(i)/abs(baseline_vec(i)))>amp_thresh
        
        % indicate whether a place cell
        place_cell_vec(i) = 1;
        
        % calculate cross_validated response
        [~,preferred_pos(i)] = max(nanmean(activity_data(half1_idx,:)));
        cross_validated_resp(i,:) = nanmean(activity_data(half2_idx,:));
        
        % plot place cell response if desired
        if plot_resp==1
            
            % calculate sem (fix NaNs at end if necessary)
            sem_resp = std(activity_data)/sqrt(size(activity_data,1))*100;
            
            if isnan(sem_resp(1))
                counter = 2; 
                while isnan(sem_resp(counter))
                    counter = counter + 1;
                end
                sem_resp(1) = sem_resp(counter); 
            end
            
            for kk = 2:length(sem_resp)
                if isnan(sem_resp(kk))
                    sem_resp(kk) = sem_resp(kk - 1);
                end
            end

            % plot mean +/- sem
            x_vec = [1:numSamp fliplr(1:numSamp)]*cm_per_bin;
            y_vec = [mean_resp+sem_resp fliplr(mean_resp-sem_resp)];
            hold off
            fill(x_vec,y_vec,[0 0.75 1],'edgecolor',[0 0.75 1])
            hold on
            plot([1:numSamp]*cm_per_bin,mean_resp,'k','linewidth',3)
            xlabel('Position (cm)')
            ylabel('DF/F')
            title(['Neuron #' num2str(i) ' place field'])
            set(gcf,'color',[1 1 1])
            hold off
            pause
            
        end
    end
    
end

% Place cells
PCs = find(place_cell_vec == 1); 

%% Plot cross-validated place cell response
pref_positions = preferred_pos(PCs);
cross_val_resp_mat = cross_validated_resp(PCs,:);
[sorted_pos,sort_idx] = sort(pref_positions,'ascend');
sorted_cross_val_resp = cross_val_resp_mat(sort_idx,:);
for i = 1:size(sorted_cross_val_resp,1)
    sorted_cross_val_resp(i,:) = smooth(sorted_cross_val_resp(i,:));
    sorted_cross_val_resp(i,:) = sorted_cross_val_resp(i,:)-min(sorted_cross_val_resp(i,:));
    sorted_cross_val_resp(i,:) = sorted_cross_val_resp(i,:)/max(sorted_cross_val_resp(i,:));
end

imagesc([1:numSamp]*cm_per_bin,1:sum(place_cell_vec),sorted_cross_val_resp)
xlabel('Position (cm)')
ylabel('Neuron (sorted)')
title('Place cells plotted by position')
set(gcf,'color',[1 1 1])
if save_flag
    saveas(gcf,save_file_name)
end

place_cell_pctl = sum(place_cell_vec)/numCell*100;
disp(['Percentage of place cells = ' num2str(place_cell_pctl,'%.1f') '%'])

%% Saving 
if save_flag==1
    save(save_file_name,'reliable_cell_vec','FWHM_vec','gof_vec','place_cell_vec', 'spatial_info');
end

    






