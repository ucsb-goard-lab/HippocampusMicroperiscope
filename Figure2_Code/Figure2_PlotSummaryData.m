%%% Figure2_PlotSummaryData
%%%
%%% Collate data and plot Fig.2D-G
%%% Written MG 220511
%%%
%%% Inputs: None
%%% Outputs: None (Data is plotted in Figure)

%% User input (files used for analysis)
filenames{1} = {'beadProfiles_v1CA1_vol1','beadProfiles_v1CA1_vol2'};
filenames{2} = {'beadProfiles_v2HPC_vol1','beadProfiles_v2HPC_vol2','beadProfiles_v2HPC_vol3'};
filenames{3} = {'beadProfiles_window_vol1','beadProfiles_window_vol2','beadProfiles_window_vol3'};
opticstype = {'v1_C_A_1','v2_H_P_C','Window'};

%% User parameters
xy_window = 75;      % number of pixels around centroid to include in PSF
z_window = 15;       % number of slices around centroid to include in PSF
int_z = 0.5;         % interval between planes (default: 1 um)
op_zoom = 16;        % optical zoom (default: 16x, max for Prairieview)
im_sz_um = 829;      % image size in microns at 1x (default: 829 um for 16x/0.8NA)
im_sz_px = 760;      % number of pixels (default: 760 px)
um_per_pix = (im_sz_um/im_sz_px)/op_zoom;  

%% Go to directory (might need to be altered depending on directory structure)
cd ..
cd('Figure2_Data')
cd('ProcessedData')

%% Collate data
for i = 1:length(filenames)
    curr_filename = filenames{i};
    xy_plane_collated = [];
    xz_plane_collated = [];
    xy_curves_collated = [];
    xz_curves_collated = [];
    x_FWHM_collated = [];
    z_FWHM_collated = [];
    for j = 1:length(curr_filename)
        load(curr_filename{j})
        xy_plane_collated = cat(3,xy_plane_collated,xy_plane);
        xz_plane_collated = cat(2,xz_plane_collated,xz_plane);
        xy_curves_collated = cat(2,xy_curves_collated,xy_curves_norm);
        xz_curves_collated = cat(2,xz_curves_collated,xz_curves_norm);
        x_FWHM_collated = [x_FWHM_collated x_FWHM];
        z_FWHM_collated = [z_FWHM_collated z_FWHM];
    end
    xy_plane_mean{i} = mean(xy_plane_collated,3);
    xz_plane_mean{i} = reshape(mean(xz_plane_collated,2),[2*xy_window+1 2*z_window+1])';
    xy_curves_mean{i} = mean(xy_curves_collated,2);
    xz_curves_mean{i} = mean(xz_curves_collated,2);
    x_FWHM_mean{i} = median(x_FWHM_collated);
    z_FWHM_mean{i} = median(z_FWHM_collated);
end
    
%% Plot summary data
x_axis = linspace(-xy_window*um_per_pix,xy_window*um_per_pix,2*xy_window+1);
z_axis = linspace(-z_window*int_z,z_window*int_z,2*z_window+1);
xy_ticklabels = round(linspace(-xy_window*um_per_pix,xy_window*um_per_pix,9)*10)/10;
xy_ticks = linspace(1,(2*xy_window+1),numel(xy_ticklabels));
z_ticklabels = round(linspace(-z_window*int_z,z_window*int_z,9)*10)/10;
z_ticks = linspace(1,(2*z_window+1),numel(z_ticklabels));

for i = 1:length(filenames)
    
    
    figure(i)
    subplot(2,2,1)
    imagesc(xy_plane_mean{i})
    set(gca,'XTick',xy_ticks,'XTickLabel',xy_ticklabels)
    set(gca,'YTick',xy_ticks,'YTickLabel',xy_ticklabels)
    title([opticstype{i} ' average bead profile (X-Y)'])
    xlabel('Lateral distance from centroid (um)')
    ylabel('Lateral distance from centroid (um)')
    axis square
    
    subplot(2,2,2)
    imagesc(xz_plane_mean{i})
    set(gca,'XTick',xy_ticks,'XTickLabel',xy_ticklabels)
    set(gca,'YTick',z_ticks,'YTickLabel',z_ticklabels)
    title([opticstype{i} ' average bead profile (X-Z)'])
    xlabel('Lateral distance from centroid (um)')
    ylabel('Axial distance from centroid (um)')
    axis square
    
    subplot(2,2,3)
    plot(x_axis,xy_curves_mean{i},'linewidth',2)
    title(['Lateral Resolution (FWHM = ' num2str(x_FWHM_mean{i},'%.1f') ' um)'])
    axis([x_axis(1) x_axis(end) 0 1])
    xlabel('Distance from centroid (um)')
    ylabel('Normalized intensity')
    axis square
    
    subplot(2,2,4)
    plot(xz_curves_mean{i},z_axis,'linewidth',2)
    axis([0 1 z_axis(1) z_axis(end)])
    set(gca,'YDir','reverse')
    title(['Axial Resolution (FWHM = ' num2str(z_FWHM_mean{i},'%.1f') ' um)'])
    xlabel('Normalized intensity')
    ylabel('Distance from centroid (um)')
    axis square
    
    set(gcf,'color',[1 1 1])
    set(gcf,'Position',[500 300 800 700])
    shg
end
