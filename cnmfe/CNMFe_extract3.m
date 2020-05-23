function [metadata] = CNMF_extract3(nam,varargin)



%% -------------CNMF-angeloIV: modified from 'demo_large_data_1p'(CNMF_e-master/demos)-------------------

%% clear the workspace and select data
%clear; clc; close all;
%% clear the workspace and select data
close all;

neuron = Sources2D();
use_prev = 0; % use previous extraction?
% nam = './Dsampled.tif';
if exist('nam') ==1;
nam = neuron.select_data(nam);         %if nam is [], then select data interactively
else
    nam = neuron.select_data;         %if nam is [], then select data interactively
end
%% parameters
% -------------------------    COMPUTATION    -------------------------  %
pars_envs = struct('memory_size_to_use', 20, ...  % GB, memory space you allow to use in MATLAB
    'memory_size_per_patch', 20, ...              % GB, space for loading data within one patch
    'patch_dims', [200, 200]);                    % px, patch size [200, 200] usually

% -------------------------      SPATIAL      -------------------------  %

metadata.cnmfe.gSig = 3;           % pixel, gaussian width of a gaussian kernel for filtering the data. 0 means no filtering
metadata.cnmfe.gSiz = 4*metadata.cnmfe.gSig+1;    % pixel, approximate neuron diameter
metadata.cnmfe.ssub = 1;           % spatial downsampling factor


% determine the search locations by selecting a round area
updateA_search_method = 'ellipse';
updateA_dist = 2;      %was 3, this parameter (sensitive!) regulates the expansion of the spatial footprints during update
updateA_bSiz = neuron.options.dist;

spatial_constraints = struct('connected', true, 'circular', false);  % you can include following constraints: 'circular'
spatial_algorithm = 'hals_thresh';

% -------------------------      TEMPORAL     -------------------------  %
try
metadata.cnmfe.Fs = 30./metadata.temporal_downsample;             % frame rate
catch
    disp(' no metadata, assuming 6fps...')
    metadata.cnmfe.Fs = 6;             % frame rate
end
metadata.cnmfe.tsub = 1;           % temporal downsampling factor
deconv_flag = true; % run deconvolution or not

%--constrained foopsi deconvolution
deconv_options = struct('type', 'ar2', ...      % model of the calcium traces. {'ar1', 'ar2'}
    'deconv_method','constrained_foopsi',...    % activity deconvolution method
    'optimize_pars', true, ...                  % optimize AR coefficients
    'temporal_iter',2,...                       % number of block-coordinate descent steps
    'fudge_factor',0.98);                       % bias correction for AR coefficients

nk = 3;                     % detrending the slow fluctuation. usually 1 is fine (no detrending)or try some integers smaller than total_frame/(Fs*30)
detrend_method = 'spline';  % compute the local minimum as an estimation of trend.

% -------------------------     BACKGROUND    -------------------------  %
bg_model = 'ring';                              % model of the background {'ring', 'svd'(default), 'nmf'}
nb = 1;                                         % number of background sources for each patch (only be used in SVD and NMF model)
bg_neuron_factor = 2;                           %this was 1.4
ring_radius = round(bg_neuron_factor * metadata.cnmfe.gSig);   % when the ring model used, it is the radius of the ring used in the background model.
num_neighbors = [];                             % number of neighbors for each neuron
metadata.cnmfe.bg_ssub = 2;                                    % downsample background for a faster speed

% -------------------------  INITIALIZATION   -------------------------  %
K = [];                         % maximum number of neurons per patch. when K=[], take as many as possible.
metadata.cnmfe.min_corr = 0.85;                % minimum local correlation for a seeding pixel
metadata.cnmfe.min_pnr = 40;                   % minimum peak-to-noise ratio for a seeding pixel (low values are discouraged)
metadata.cnmfe.min_pixel = 6*metadata.cnmfe.gSig^2;           % minimum number of nonzero pixels for each neuron (was gSig^2)
bd = 8;                         % number of rows/columns to be ignored in the boundary (mainly for motion corrected data)
frame_range = [];               % when [], uses all frames
save_initialization = false;    % save the initialization procedure as a video.
use_parallel = true;            % use parallel computation for parallel computing
show_init = true;               % show initialization results
center_psf = true;              % set the value as true when the background fluctuation is large (usually 1p data)

% -------------------------      MERGING      -------------------------  %
show_merge = false;                    % if true, manually verify the merging step
metadata.cnmfe.merge_thr = 0.8;                       % temporal correlation threshold for merging neurons;
method_dist = 'mean';                  % method for computing neuron distances {'mean', 'max'}
dmin = 0.5*metadata.cnmfe.gSiz;                       % minimum distances between two neurons. it is used together with merge_thr (it was 1*gSiz)
dmin_only = metadata.cnmfe.gSig;                      % merge neurons if their distances are smaller than dmin_only.
metadata.cnmfe.merge_thr_spatial = [0.9, 0.6, -inf];  % merge components with highly correlated spatial shapes and temporal shapes

% -------------------------  Residual   -------------------------  %
metadata.cnmfe.min_corr_res = 0.6;
metadata.cnmfe.min_pnr_res = 10;
seed_method_res = 'manual';  % method for initializing neurons from the residual
update_sn = true;

% -------------------------  FINAL RESULTS   -------------------------  %
save_demixed = true;    % save the demixed file or not
kt = 3;                 % frame intervals

% -------------------------    UPDATE ALL    -------------------------  %
% Be aware that what is declared here won't change unless you call again
% the function neuron.updateParams()
neuron.updateParams('gSig', metadata.cnmfe.gSig, ...       % -------- spatial --------
    'gSiz', metadata.cnmfe.gSiz, ...
    'ring_radius', ring_radius, ...
    'ssub', metadata.cnmfe.ssub, ...
    'search_method', updateA_search_method, ...
    'bSiz', updateA_bSiz, ...
    'dist', updateA_dist, ...
    'spatial_constraints', spatial_constraints, ...
    'spatial_algorithm', spatial_algorithm, ...
    'tsub', tsub, ...                       % -------- temporal --------
    'deconv_flag', deconv_flag, ...
    'deconv_options', deconv_options, ...
    'nk', nk, ...
    'detrend_method', detrend_method, ...
    'background_model', bg_model, ...       % -------- background --------
    'nb', nb, ...
    'ring_radius', ring_radius, ...
    'num_neighbors', num_neighbors, ...
    'bg_ssub', metadata.cnmfe.bg_ssub, ...
    'merge_thr', metadata.cnmfe.merge_thr, ...             % -------- merging ---------
    'dmin', dmin, ...
    'method_dist', method_dist, ...
    'min_corr', metadata.cnmfe.min_corr, ...               % ----- initialization -----
    'min_pnr', metadata.cnmfe.min_pnr, ...
    'min_pixel', metadata.cnmfe.min_pixel, ...
    'bd', bd, ...
    'noise_method','logmexp',...
    'center_psf', center_psf);
neuron.Fs = Fs;

%% distribute data and be ready to run source extraction
neuron.getReady(pars_envs);

%% initialize neurons from the video data within a selected temporal range
[center, Cn, PNR] = neuron.initComponents_parallel(K, frame_range, save_initialization, use_parallel);
neuron.compactSpatial();
if show_init
    figure();
    ax_init= axes();
    imagesc(Cn, [0, 1]); colormap gray;
    hold on;
    plot(center(:, 2), center(:, 1), '.r', 'markersize', 10);
end
%% Show correlation and PNR images
% show correlation image, peak-to-noise ratio and pointwise product of correlation image and peak-to-noise ratio
figure('position', [10, 500, 1776, 400]);
subplot(131);   imagesc(Cn, [0, 1]); colorbar;
axis equal off tight;   title('correlation image');
subplot(132);   imagesc(PNR,[0,max(PNR(:))*0.98]); colorbar;
axis equal off tight;   title('peak-to-noise ratio');
subplot(133);   imagesc(Cn.*PNR, [0,max(PNR(:))*0.98]); colorbar;
axis equal off tight;   title('Cn*PNR');

% save PNR and correlation images as tiff
options.overwrite = true;
movie_name = cell2mat(extractBetween(nam,pwd,".tif"));  movie_name = movie_name(2:end);
saveastiff(uint16(PNR),['PNR_',movie_name,'.tif'],options);
saveastiff(uint8(Cn*255),['Cn_',movie_name,'.tif'],options);

%% estimate the background components
neuron.update_background_parallel(use_parallel);
neuron_0 = neuron.copy();

%%  merge neurons and update spatial/temporal components
neuron.merge_neurons_dist_corr(show_merge);
neuron.merge_high_corr(show_merge, metadata.cnmfe.merge_thr_spatial);
neuron_1 = neuron.copy();

%% udpate spatial&temporal components, delete false positives and merge neurons
% update spatial
if update_sn
    neuron.update_spatial_parallel(use_parallel, true);
    udpate_sn = false;
else
    neuron.update_spatial_parallel(use_parallel);
end
% merge neurons based on correlations
neuron.merge_high_corr(show_merge, metadata.cnmfe.merge_thr_spatial);

for m=1:2
    % update temporal
    neuron.update_temporal_parallel(use_parallel);

    % delete bad neurons
    neuron.remove_false_positives();

    % merge neurons based on temporal correlation + distances
    neuron.merge_neurons_dist_corr(show_merge);
end
neuron_2 = neuron.copy();
neuron.updateParams('merge_thr', 0.7, 'dmin', 0.3*metadata.cnmfe.gSiz);
neuron.options.spatial_algorithm = 'nnls';

%% run more iterations
neuron.update_background_parallel(use_parallel);
neuron.update_spatial_parallel(use_parallel);
neuron.update_temporal_parallel(use_parallel);

K = size(neuron.A,2);
tags = neuron.tag_neurons_parallel();  % find neurons with fewer nonzero pixels than min_pixel and silent calcium transients
neuron.remove_false_positives();
neuron.merge_neurons_dist_corr(show_merge);
neuron.merge_high_corr(show_merge, metadata.cnmfe.merge_thr_spatial);

if K~=size(neuron.A,2)
    neuron.update_spatial_parallel(use_parallel);
    neuron.update_temporal_parallel(use_parallel);
    neuron.remove_false_positives();
end
neuron_3 = neuron.copy();

%% Increase the quality of fitting
%might need to modify constrained_oasisAR1
better_fit = true;
if better_fit
    options = neuron.options.deconv_options;
    for i=1:size(neuron.C,1)
        [c, s, option] = deconvolveCa(neuron.C_raw(i,:), options, 'ar1');
        pars = [option.pars, 0];
        [c2, s2, option] = deconvolveCa(neuron.C_raw(i,:), options, 'pars', pars);
        neuron.C(i,:) = c2;
    end
end
neuron_4 = neuron.copy();

%% Last merging
%Merge neurons
neuron.updateParams('merge_thr', 0.60,'dmin', 0.7*metadata.cnmfe.gSiz);
neuron.merge_neurons_dist_corr(false);  %(..., temporal correlation, min distance)

neuron.compactSpatial();
neuron.update_background_parallel(use_parallel);
neuron.update_spatial_parallel(use_parallel);
neuron.update_temporal_parallel(use_parallel);

%Increase the quality of fitting
options = neuron.options.deconv_options;
for i=1:size(neuron.C,1)
    [c, s, option] = deconvolveCa(neuron.C_raw(i,:), options, 'ar1');
    pars = [option.pars, 0];
    [c2, s2, option] = deconvolveCa(neuron.C_raw(i,:), options, 'pars', pars);
    neuron.C(i,:) = c2;
end

%% Save workspace
neuron.orderROIs('snr');
%cnmfe_path = neuron.save_workspace();
movie_name = cell2mat(extractBetween(nam,pwd,".tif"));  movie_name = movie_name(2:end);
filename = ['CNMFe_',datestr(now, 'yymmdd_HHMM'),'_raw_',movie_name,'.mat'];
save(filename);

%% Show current situation
figure();   set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
subplot(2,3,1); plot_contours(neuron_0.A,PNR,0.9,true);
subplot(2,3,2); plot_contours(neuron_1.A,PNR,0.9,true);
subplot(2,3,3); plot_contours(neuron_2.A,PNR,0.9,true);
subplot(2,3,4); plot_contours(neuron_3.A,PNR,0.9,true);
subplot(2,3,5); plot_contours(neuron_4.A,PNR,0.9,true);
subplot(2,3,6); plot_contours(neuron.A,PNR,0.9,true);

%% -------------Post extraction visualization----------------
figure();   set(gcf, 'units','normalized','outerposition',[0 0.4 0.5 0.6]);
plot_contours(neuron.A,PNR,0.9,true);   %plot_contours(Aor,Cn,thr,display_numbers,max_number,Coor, ln_wd, ln_col)
figure();   set(gcf, 'units','normalized','outerposition',[0.1 0.05 0.3 0.4]);
r = corrcoef(neuron.C');   imagesc(r); colorbar;

%% Check proper demixing (if needed)
% neuron_i = 2;
% neuron_j = 12;
%
% figure();
% plot(neuron.C(neuron_i,:));   hold on;    plot(neuron.C(neuron_j,:));   hold off;
%
% figure()
% for i = 1:10
% imagesc(neuron.reshape(neuron.A(:, neuron_i), 2));    pause(0.3);
% imagesc(neuron.reshape(neuron.A(:, neuron_j), 2));    pause(0.3);
% end

%% -------------Save Data----------------


results = neuron.obj2struct();
results.metadata = metadata;
save('results.mat','results');

ImBat_ROIoverlay(results);
print('ROI_Images/overlay','-dpng')
