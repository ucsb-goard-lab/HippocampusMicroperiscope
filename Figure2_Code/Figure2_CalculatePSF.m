%%% Figure2_CalculatePSF
%%%
%%% Code used for: Fig.2D-G, 2I, 2K of Redman et al (eLife, 2022).
%%%
%%% Description: Measure the point spread function of fluorescent 
%%% microspheres. Fluorescent microspheres (beads) should be smaller than 
%%% theoretical resolution of microscope (typically 0.1-0.2 um).
%%%
%%% Written: MG 200816, last update: 220511
%%%
%%% Inputs: 
%%%     None - user will be prompted to load file.
%%%     File should be first TIF file from Z stack of microspheres
%%%     (may need to change file indexing to make compatible with new data)
%%% Outputs:
%%%     xy_plane: XY image through center of each bead (X*Y*cell matrix)
%%%     xz_plane: XZ image through center of each bead (X*cell*Z matrix)
%%%     xy_curves/xz_curves: curves through Y centroid for X and Z
%%%     xy_curves_norm/xz_curves_norm: as above, but normalized 0-1
%%%     x_FWHM/z_FWHM: calculated full width at half max for X and Z curves
%%%     x_fits/z_fits: quality of gaussian fits

%% Initialize workspace
clear
close all

%% Imaging parameters (User input)
num_z = 101;        % number of image planes (default: 101)
int_z = 0.5;        % interval between planes (default: 0.5 um)
op_zoom = 16;       % optical zoom (default: 16x, max for Prairieview)
im_sz_um = 829;     % image size in microns at 1x (default: 829 um for 16x/0.8NA)
im_sz_px = 760;     % number of pixels (default: 760 px)
lambda = 920;       % wavelength (default: 920 nm)
NA_obj = 0.8;       % NA of the objective (default: Nikon 16x = 0.8)
ref_index = 1.33;   % refractive index (default: water = 1.33)

%% User parameters
max_bead_mult = 3;   % multiples of theoretical lateral bead area (squared) for max (default: 2.5-3.0)
min_neighbor = 25;   % closest neaerest neighbor distance (default: 25)
xy_window = 75;      % number of pixels around centroid to include in PSF (default: 75)
z_window = 15;       % number of slices around centroid to include in PSF (default: 15)
gof_thresh = 0.5;    % minimum r^2 value on gaussian (default: 0.5)

%% Calculate theoretical maximum resolution (Zipfel et al., Nat Biotech, 2003)
lat_res = (2*sqrt(log(2)) * 0.325*lambda) / (1000*sqrt(2)*NA_obj^0.91); % lat res
axial_res = (2*sqrt(log(2)) * ((0.532*lambda) / (1000*sqrt(2))))...
    * (1/(ref_index-sqrt(ref_index^2-NA_obj^2))); % axial res
disp('Theoretical diffraction-limited resolution:')
disp(['Lateral: ' num2str(lat_res) ' um'])
disp(['Axial: ' num2str(axial_res) ' um'])

%% calculate min and max bead size (to eliminate artifacts or clumped beads)
um_per_pix = (im_sz_um/im_sz_px)/op_zoom;          % micron per pixel
bead_diam = lat_res/um_per_pix;                    % smallest diam in pixels
area_per_bead = pi*(bead_diam/2)^2;                % minimum theoretical bead area
minBeadSize = round(area_per_bead);                % min bead size (exclude smaller beads)
maxBeadSize = round(max_bead_mult^2*minBeadSize);  % max bead size (avoid multibead clusters)

%% Load files (may need to change indexing if using with other data)
 disp('Load first file')
[filename,~] = uigetfile('.tif');
for i = 1:num_z
    filename(end-10:end-8) = sprintf('%03d',i);
    imageStack(:,:,i) = imread(filename,1);
end 

%% Detect and segment potential beads
% calcualte average projection
avg_projection = mean(imageStack,3);
avg_proj_norm = (avg_projection-min(min(avg_projection)));
avg_proj_norm = avg_proj_norm/max(max(avg_projection));

% Segmentation instructions
disp('Image Segmenter opening, follow the directions below:')
disp('1) Under "Create mask tab, select "Find Circles"')
disp('2) Set Min Diameter = 4, Max Diameter = 50, Sensitivity = 0.9')
disp('3) Click "Find Circles" (green arrow) and wait for all circles to populate')
disp('4) Click "Create mask" (check mark)')
disp('5) Click "Export" (check mark) and then "OK" to save the masks')

% Segment
imageSegmenter(avg_proj_norm);
f = findall(0,'name','Segmentation');
while exist('BW','var')==0
    pause(0.1)
end
imageSegmenter close
segmentedROIs = bwlabel(BW);
numROIs = max(max(segmentedROIs));

% Find centroids 
props = regionprops(segmentedROIs,avg_projection,'Centroid','WeightedCentroid');
centroids = vertcat(props.Centroid);

% display
imshow(avg_proj_norm)
hold on
for i = 1:numROIs
    plot(centroids(i,1),centroids(i,2),'r.','MarkerSize',15)
end
set(gcf,'color',[1 1 1])
title('Check centroids: large ROIs might be innacurate (they will be removed)')
shg
pause
close

% show bead size histogram
beadSize = zeros(1,numROIs);
for i = 1:numROIs
    curr_mask = (segmentedROIs==i);
    beadSize(i) = sum(sum(curr_mask));
end
[N,X] = hist(beadSize,100);
bar(X,N)
hold on
plot(ones(1,50)*minBeadSize,linspace(0,ceil(max(N)/10)*10,50),'r:','linewidth',2)
plot(ones(1,50)*maxBeadSize,linspace(0,ceil(max(N)/10)*10,50),'r:','linewidth',2)
set(gcf,'color',[1 1 1])
title('Bead size histogram')
shg
pause
close

%% Find well-isolated beads 
% 1) Exclude ROIs on edge if image
% 2) Exclude ROIs that are too close to other ROIs
% 3) Exclude ROIs whose area is >mult_beads x difraction-limited bead size (non-singlet)

% initialize
selectFlag = zeros(numROIs,3);
selectBeads = zeros(im_sz_px);
beadCount = 0;

% find distances between centroids
dist_vec = pdist(centroids,'euclidean');
dist_mat = squareform(dist_vec);

% calculate nearest neighbor distance for each centroid
nearest_neighbor = zeros(1,numROIs);
for i = 1:numROIs
    nearest_neighbor(i) = min(dist_mat([1:i-1 i+1:numROIs],i));
end

% check criteria for each bead
for i = 1:numROIs
    
    % check that ROI is far enough from edge
    ROI_x = round(centroids(i,1));
    ROI_y = round(centroids(i,2));
    Z_profile = squeeze(max(imageStack(ROI_y,ROI_x,:),[],2));
    [~,max_idx] = max(Z_profile);
    if (ROI_x-min_neighbor)>0 &&...
            (ROI_x+min_neighbor)<im_sz_px &&...
            (ROI_y-min_neighbor)>0 &&...
            (ROI_y+min_neighbor)<im_sz_px &&...
            (max_idx-z_window)>1 &&...
            (max_idx+z_window)<num_z
        selectFlag(i,1) = 1;
    end
    
    % check that ROI is far enough from other ROIs
    if nearest_neighbor(i)>=min_neighbor
        selectFlag(i,2) = 1;
    end
    
    % check that ROI corresponds to a single bead
    if beadSize(i)>=minBeadSize && beadSize(i)<=maxBeadSize
        selectFlag(i,3) = 1;
    end
    
    % keep beads that satisfy all parameters
    if sum(selectFlag(i,:))==3 
        beadCount = beadCount+1;
        selectBeads(segmentedROIs==i) = beadCount;
        selectCentroids(beadCount,:) = centroids(i,:);
    end
    
end

% display
imshow(avg_proj_norm)
hold on
for i = 1:beadCount
    plot(selectCentroids(i,1),selectCentroids(i,2),'r.','MarkerSize',15)
end
set(gcf,'color',[1 1 1])
title('Check centroids: If misaligned - adjust Min Diameter, Max Diameter, and Sensitivity paraemters.')
shg
pause
close

% Calculate average point spread function across beads
xy_plane = zeros(2*xy_window+1,2*xy_window+1,beadCount);
xz_plane = zeros(2*xy_window+1,beadCount,2*z_window+1);

% Extract xy and xz planes and curves from all good beads
for i = 1:beadCount
    x_coor = round(selectCentroids(i,1));
    y_coor = round(selectCentroids(i,2));
        
    x_pixels = [x_coor-xy_window+im_sz_px:x_coor+xy_window+im_sz_px];
    y_pixels = [y_coor-xy_window+im_sz_px:y_coor+xy_window+im_sz_px];
    
    [~,max_idx] = max(imageStack(y_coor,x_coor,:));
    xy_plane_rep = repmat(imageStack(:,:,max_idx),3);
    xy_plane(:,:,i) = xy_plane_rep(y_pixels,x_pixels);
    xy_curves(:,i) = mean([xy_plane(:,xy_window+1,i) xy_plane(xy_window+1,:,i)'],2);
    xz_plane_rep = repmat(imageStack(:,:,:),[1 3 1]);
    xz_plane(:,i,:) = xz_plane_rep(y_coor,x_pixels,[max_idx-z_window:max_idx+z_window]);
    xz_curves(:,i) = xz_plane(xy_window+1,i,:);
end

% Normalize curves
x_axis = linspace(-xy_window*um_per_pix,xy_window*um_per_pix,2*xy_window+1);
z_axis = linspace(-z_window*int_z,z_window*int_z,2*z_window+1);
xy_curves_norm = zeros(size(xy_curves));
xz_curves_norm = zeros(size(xz_curves));
for i = 1:beadCount
    [KSD,Xi] = ksdensity(reshape(xy_plane(:,:,i),[1 size(xy_plane,1)*size(xy_plane,2)]));
    [~,maxIdx]= max(KSD);
    xy_baseline = Xi(maxIdx);
    xy_curves_norm(:,i) = xy_curves(:,i)-xy_baseline;
    xy_curves_norm(:,i) = xy_curves_norm(:,i)/max(xy_curves_norm(:,i));

    [KSD,Xi] = ksdensity(reshape(xz_plane(:,i,:),[1 size(xz_plane,1)*size(xz_plane,3)]));
    [~,maxIdx]= max(KSD);
    xz_baseline = Xi(maxIdx);
    xz_curves_norm(:,i) = xz_curves(:,i)-xz_baseline;
    xz_curves_norm(:,i) = xz_curves_norm(:,i)/max(xz_curves_norm(:,i));
end

% Fit with gaussians and calculate FWHM
goodFit_idx = [];
for i = 1:beadCount
    curr_curve = xy_curves_norm(:,i);
    fit_func = fittype('exp(-((x-b1)/c1)^2)');     % define function: gaussian normlaized to amp = 1
    initial_param = [0 lat_res];                   % initial parameters: [b1 c1]
    [fitobj,gof] = fit(x_axis',curr_curve,fit_func,'StartPoint',initial_param);
    x_FWHM(i) = 2*sqrt(log(2))*fitobj.c1;
    x_fits(:,i) = fitobj(x_axis);
    x_gof(i) = gof.rsquare;

    curr_curve = xz_curves_norm(:,i);
    fit_func = fittype('exp(-((x-b1)/c1)^2)');     % define function: gaussian normlaized to amp = 1
    initial_param = [0 axial_res];                 % initial parameters: [b1 c1]
    [fitobj,gof] = fit(z_axis',curr_curve,fit_func,'StartPoint',initial_param);
    z_FWHM(i) = 2*sqrt(log(2))*fitobj.c1;
    z_fits(:,i) = fitobj(z_axis);
    z_gof(i) = gof.rsquare;
    
    if x_gof(i) > gof_thresh && z_gof(i) > gof_thresh
        goodFit_idx = [goodFit_idx i];
    end
end

% save beads with good gaussian fits
xy_plane = xy_plane(:,:,goodFit_idx);
xz_plane = xz_plane(:,goodFit_idx,:);
xy_curves = xy_curves(:,goodFit_idx);
xz_curves = xz_curves(:,goodFit_idx);
xy_curves_norm = xy_curves_norm(:,goodFit_idx);
xz_curves_norm = xz_curves_norm(:,goodFit_idx);
x_FWHM = x_FWHM(goodFit_idx);
z_FWHM = z_FWHM(goodFit_idx);
x_fits = x_fits(:,goodFit_idx);
z_fits = z_fits(:,goodFit_idx);

% Calculate means
mean_xy = mean(xy_plane,3);
xz_mat = mean(xz_plane,2);
mean_xz = reshape(xz_mat,[2*xy_window+1 2*z_window+1])';
x_curve_avg = mean(xy_curves_norm,2);
z_curve_avg = mean(xz_curves_norm,2);

%% Plot

% set axes
x_axis = linspace(-xy_window*um_per_pix,xy_window*um_per_pix,2*xy_window+1);
z_axis = linspace(-z_window*int_z,z_window*int_z,2*z_window+1);
xy_ticklabels = round(linspace(-xy_window*um_per_pix,xy_window*um_per_pix,9)*10)/10;
xy_ticks = linspace(1,(2*xy_window+1),numel(xy_ticklabels));
z_ticklabels = round(linspace(-z_window*int_z,z_window*int_z,9)*10)/10;
z_ticks = linspace(1,(2*z_window+1),numel(z_ticklabels));

subplot(2,2,1)
imagesc(mean_xy)
set(gca,'XTick',xy_ticks,'XTickLabel',xy_ticklabels)
set(gca,'YTick',xy_ticks,'YTickLabel',xy_ticklabels)
title('Average bead profile (X-Y)')
xlabel('Lateral distance from centroid (um)')
ylabel('Lateral distance from centroid (um)')
axis square

subplot(2,2,2)
imagesc(mean_xz)
set(gca,'XTick',xy_ticks,'XTickLabel',xy_ticklabels)
set(gca,'YTick',z_ticks,'YTickLabel',z_ticklabels)
title('Average bead profile (X-Z)')
xlabel('Lateral distance from centroid (um)')
ylabel('Axial distance from centroid (um)')
axis square

subplot(2,2,3)
plot(x_axis,x_curve_avg,'linewidth',2)
title(['Lateral Resolution (FWHM = ' num2str(median(x_FWHM),'%.1f') ' um)'])
axis([x_axis(1) x_axis(end) 0 1])
xlabel('Distance from centroid (um)')
ylabel('Normalized intensity')
axis square

subplot(2,2,4)
plot(z_curve_avg,z_axis,'linewidth',2)
axis([0 1 z_axis(1) z_axis(end)])  
set(gca,'YDir','reverse')
title(['Axial Resolution (FWHM = ' num2str(median(z_FWHM),'%.1f') ' um)'])
xlabel('Normalized intensity')
ylabel('Distance from centroid (um)')
axis square

% set figure size/color/position and display
set(gcf,'color',[1 1 1])
set(gcf,'Position',[500 300 800 700])
shg

%% Store variabels in data structure
PSF_data.xy_plane = xy_plane;
PSF_data.xz_plane = xz_plane;
PSF_data.xy_curves = xy_curves;
PSF_data.xz_curves = xz_curves;
PSF_data.xy_curves_norm = xy_curves_norm; 
PSF_data.xz_curves_norm = xz_curves_norm;
PSF_data.x_FWHM = x_FWHM;
PSF_data.z_FWHM = z_FWHM;
PSF_data.x_fits = x_fits; 
PSF_data.z_fits = z_fits;
clearvars -except PSF_data

