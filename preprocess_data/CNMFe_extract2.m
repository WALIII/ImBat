function [metadata] = CNMFe_extract2(nam,varargin)
%functionalized from CNMF-angeloII: modified from 'demo_large_data_1p'(CNMF_e-master/demos)


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

%% select data and map it to the RAM
cnmfe_choose_data;

%% parameters

% % Pass  Manual inputs
vin=varargin;
for i=1:length(vin)
    if isequal(vin{i},'roi') % manually inputing a sort order
        ROI_flag=vin{i+1};
    elseif isequal(vin{i},'place');
        analysis_flag = vin{i+1};
      elseif isequal(vin{i},'metadata');  % pass along metadata.cnmfe.file if need be...
          metadata = vin{i+1};
          disp('Loading in metadata.cnmfe...')
    elseif isequal(vin{i},'rextract');
        reExtract=vin{i+1};
    end
end



% -------------------------    COMPUTATION    -------------------------  %
pars_envs = struct('memory_size_to_use', 60,...  % GB, memory space you allow to use in MATLAB
    'memory_size_per_patch', 40, ...              % GB, space for loading data within one patch
    'patch_dims', [200, 200]);                    % px, patch size [200, 200] usually

% -------------------------      SPATIAL      -------------------------  %
metadata.cnmfe.gSig = 4;           % pixel, gaussian width of a gaussian kernel for filtering the data. 0 means no filtering
metadata.cnmfe.gSiz = 4*metadata.cnmfe.gSig+1;    % pixel, approximate neuron diameter

with_dendrites = false;   % with dendrites or not, indicated with high-resolution/2p-movies
if with_dendrites
    % determine the search locations by dilating the current neuron shapes
    updateA_search_method = 'dilate';  %#ok<UNRCH>
    updateA_bSiz = 5;
    updateA_dist = neuron.options.dist;
else
    % determine the search locations by selecting a round area
    updateA_search_method = 'ellipse'; %#ok<UNRCH>
    updateA_dist = 3;      %this was 3
    updateA_bSiz = neuron.options.dist;
end
spatial_constraints = struct('connected', true, 'circular', false);  % you can include following constraints: 'circular'
spatial_algorithm = 'hals_thresh';

% -------------------------      TEMPORAL     -------------------------  %
try
metadata.cnmfe.Fs = 30./metadata.temporal_downsample;             % frame rate
catch
    disp(' no metadata, assuming 6fps...')
    metadata.cnmfe.Fs = 6;             % frame rate
end
deconv_flag = true; % run deconvolution or not

% %--foopsi deconvolution
% deconv_options = struct('type', 'ar2', ... % model of the calcium traces. {'ar1', 'ar2'}
%     'method', 'foopsi', ... % method for running deconvolution {'foopsi', 'constrained_foopsi(problematic A.F.)', 'thresholded'}
%     'smin', -5, ...         % minimum spike size. When the value is negative, the actual threshold is abs(smin)*noise level
%     'optimize_pars', true, ...  % optimize AR coefficients
%     'optimize_b', true, ...% optimize the baseline);
%     'max_tau', [100]);    % maximum decay time (unit: frame);

%--constrained foopsi deconvolution
deconv_options = struct('type', 'ar2', ...      % model of the calcium traces. {'ar1', 'ar2'}
    'deconv_method','constrained_foopsi',...    % activity deconvolution method
    'optimize_pars', true, ...                  % optimize AR coefficients
    'temporal_iter',2,...                       % number of block-coordinate descent steps
    'fudge_factor',0.98);                          % bias correction for AR coefficients
metadata.cnmfe.deconv_options = deconv_options

nk = 3;                     % detrending the slow fluctuation. usually 1 is fine (no detrending)or try some integers smaller than total_frame/(metadata.cnmfe.Fs*30)
detrend_method = 'spline';  % compute the local minimum as an estimation of trend.

% -------------------------     BACKGROUND    -------------------------  %
bg_model = 'ring';                              % model of the background {'ring', 'svd'(default), 'nmf'}
nb = 1;                                         % number of background sources for each patch (only be used in SVD and NMF model)
bg_neuron_factor = 2;                           %this was 1.4
ring_radius = round(bg_neuron_factor * metadata.cnmfe.gSig);   % when the ring model used, it is the radius of the ring used in the background model.
num_neighbors = [];                             % number of neighbors for each neuron
bg_ssub = 2;                                    % downsample background for a faster speed

% -------------------------  INITIALIZATION   -------------------------  %
K = [];                         % maximum number of neurons per patch. when K=[], take as many as possible.
%metadata.cnmfe.min_corr = 0.9;                 % minimum local correlation for a seeding pixel
%metadata.cnmfe.min_pnr = 50;                  % minimum peak-to-noise ratio for a seeding pixel (low values are discouraged)
metadata.cnmfe.min_pixel = metadata.cnmfe.gSig^2;             % minimum number of nonzero pixels for each neuron (was gSig^2)
bd = 6;                         % number of rows/columns to be ignored in the boundary (mainly for motion corrected data)
frame_range = [];               % when [], uses all frames
save_initialization = false;    % save the initialization procedure as a video.
use_parallel = true;            % use parallel computation for parallel computing
show_init = true;               % show initialization results
choose_params = false;          % manually choose parameters
center_psf = true;              % set the value as true when the background fluctuation is large (usually 1p data)

% -------------------------      MERGING      -------------------------  %
show_merge = false;                      % if true, manually verify the merging step
metadata.cnmfe.merge_thr = 0.65;                         % temporal correlation threshold for merging neurons;
method_dist = 'mean';                    % method for computing neuron distances {'mean', 'max'}
dmin = metadata.cnmfe.gSig*2;                             % minimum distances between two neurons. it is used together with merge_thr
dmin_only = metadata.cnmfe.gSig;                        % merge neurons if their distances are smaller than dmin_only.
metadata.cnmfe.merge_thr_spatial = [0.9, 0.6, -inf];    % merge components with highly correlated spatial shapes and temporal shapes

% -------------------------  Residual   -------------------------  %
metadata.cnmfe.min_corr_res = 0.6;
metadata.cnmfe.min_pnr_res = 10;
seed_method_res = 'manual';  % method for initializing neurons from the residual
update_sn = true;

% ----------------------  WITH MANUAL INTERVENTION  --------------------  %
with_manual_intervention = false;

% -------------------------  FINAL RESULTS   -------------------------  %
save_demixed = true;    % save the demixed file or not
kt = 3;                 % frame intervals





% -------------------------    UPDATE ALL    -------------------------  %
% Be aware that what is declared here won't change unless you call again
% the function neuron.updateParams()
neuron.updateParams('gSig', metadata.cnmfe.gSig, ...       % -------- spatial --------
    'gSiz', metadata.cnmfe.gSig, ...
    'ring_radius', ring_radius, ...
    'ssub', metadata.cnmfe.ssub, ...
    'search_method', updateA_search_method, ...
    'bSiz', updateA_bSiz, ...
    'dist', updateA_bSiz, ...
    'spatial_constraints', spatial_constraints, ...
    'spatial_algorithm', spatial_algorithm, ...
    'tsub', metadata.cnmfe.tsub, ...                       % -------- temporal --------
    'deconv_flag', deconv_flag, ...
    'deconv_options', deconv_options, ...
    'nk', nk, ...
    'detrend_method', detrend_method, ...
    'background_model', bg_model, ...       % -------- background --------
    'save_intermediate  ' ,true, ... % save intermediate results or not
    'nb', nb, ...
    'ring_radius', ring_radius, ...
    'num_neighbors', num_neighbors, ...
    'bg_ssub', bg_ssub, ...
    'merge_thr', metadata.cnmfe.merge_thr, ...             % -------- merging ---------
    'dmin', dmin, ...
    'method_dist', method_dist, ...
    'min_corr', metadata.cnmfe.min_corr, ...               % ----- initialization -----
    'min_pnr', metadata.cnmfe.min_pnr, ...
    'min_pixel', metadata.cnmfe.min_pixel, ...
    'bd', bd, ...
    'noise_method','logmexp',...
    'center_psf', center_psf);
neuron.Fs = metadata.cnmfe.Fs;

%% distribute data and be ready to run source extraction
neuron.getReady(pars_envs);

%% initialize neurons from the video data within a selected temporal range
if choose_params
    % change parameters for optimized initialization
    [metadata.cnmfe.gSig, metadata.cnmfe.gSig, ring_radius, metadata.cnmfe.min_corr, metadata.cnmfe.min_pnr] = neuron.set_parameters();
end

[center, Cn, PNR] = neuron.initComponents_parallel(K, frame_range, save_initialization, use_parallel, use_prev);
neuron.compactSpatial();
if show_init
    figure();
    ax_init= axes();
    imagesc(Cn, [0, 1]); colormap gray;
    hold on;
    plot(center(:, 2), center(:, 1), '.r', 'markersize', 10);
end
%% Show correlation and PNR images

% show correlation image
figure('position', [10, 500, 1776, 400]);
subplot(131);   imagesc(Cn, [0, 1]); colorbar;
axis equal off tight;   title('correlation image');

% show peak-to-noise ratio
subplot(132);   imagesc(PNR,[0,max(PNR(:))*0.98]); colorbar;
axis equal off tight;   title('peak-to-noise ratio');

% show pointwise product of correlation image and peak-to-noise ratio
subplot(133);   imagesc(Cn.*PNR, [0,max(PNR(:))*0.98]); colorbar;
axis equal off tight;   title('Cn*PNR');

% save PNR and correlation images as tiff
options.overwrite = true;
mkdir('ROI_Images');
saveastiff(uint16(PNR),'ROI_Images/PNR.tif',options);
saveastiff(uint8(Cn*255),'ROI_Images/Cn.tif',options);

%% estimate the background components
neuron.update_background_parallel(use_parallel);
neuron_init = neuron.copy();

%%  merge neurons and update spatial/temporal components
neuron.merge_neurons_dist_corr(show_merge);
neuron.merge_high_corr(show_merge, metadata.cnmfe.merge_thr_spatial);

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

%% add a manual intervention and run the whole procedure for a second time
neuron.options.spatial_algorithm = 'nnls';
if with_manual_intervention
    show_merge = true;
    neuron.orderROIs('snr');   % order neurons in different ways {'snr', 'decay_time', 'mean', 'circularity'}
    neuron.viewNeurons([], neuron.C_raw);

    % merge closeby neurons
    neuron.merge_close_neighbors(true, dmin_only);

    % delete neurons
    tags = neuron.tag_neurons_parallel();  % find neurons with fewer nonzero pixels than min_pixel and silent calcium transients
    ids = find(tags>0);
    if ~isempty(ids)
        neuron.viewNeurons(ids, neuron.C_raw);
    end
end

% %%
% deconv_options.type = 'ar2';
% neuron.updateParams('deconv_options', deconv_options);

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

%Increase the quality of fitting and store in neuron_1
%might need to modify constrained_oasisAR1
better_fit = true;
if better_fit
    neuron_1 = neuron.copy();
    options = neuron_1.options.deconv_options;
    for i=1:size(neuron_1.C,1)
        [c, s, option] = deconvolveCa(neuron_1.C_raw(i,:), options, 'ar1');
        pars = [option.pars, 0];
        [c2, s2, option] = deconvolveCa(neuron_1.C_raw(i,:), options, 'pars', pars);
        neuron_1.C(i,:) = c2;
    end
end

%% save the workspace for future analysis
neuron.orderROIs('snr');
%cnmfe_path = neuron.save_workspace();

%% show neuron contours
Coor = neuron.show_contours(0.6);

%% Save workspace
movie_name = cell2mat(extractBetween(nam,pwd,".tif"));  movie_name = movie_name(2:end);
filename = ['CNMFe_',datestr(now, 'yymmdd_HHMM'),'_raw_',movie_name,'.mat'];
save(filename);

% %% create a video for displaying the results
% amp_ac = 140;
% range_ac = 5+[0, amp_ac];
% multi_factor = 10;
% range_Y = 1300+[0, amp_ac*multi_factor];
% avi_filename = neuron.show_demixed_video(save_demixed, kt, [], amp_ac, range_ac, range_Y, multi_factor);

%% -------------Post extraction visualization----------------

% Look at the contours of neurons and their correlation
figure();   set(gcf, 'units','normalized','outerposition',[0 0.4 0.5 0.6]);
plot_contours(neuron.A,Cn,0.75,true,size(neuron.C_raw,1));   %plot_contours(Aor,Cn,thr,display_numbers,max_number,Coor, ln_wd, ln_col)
figure();   set(gcf, 'units','normalized','outerposition',[0.1 0.05 0.3 0.4]);
r = corrcoef(neuron.C_raw');   imagesc(r); colorbar;

%% ------------Run Sanity checks------------------

sanity_check = false;
residual_check = true;
if sanity_check

    %Merge neurons (run this more times-if needed-and play with thresholds)
    neuron.updateParams('merge_thr', 0.6,'dmin', 1.5*metadata.cnmfe.gSig);
    neuron.merge_neurons_dist_corr(true);  %(..., temporal correlation, min distance)

    %Pick neurons from the residual
    if residual_check
        [center_res, Cn_res, PNR_res] =neuron.initComponents_residual_parallel([], save_initialization, use_parallel, metadata.cnmfe.min_corr_res, metadata.cnmfe.min_pnr_res, seed_method_res);
        if show_init    %This adds new initialization dots (green) to the initial figure with red dots
            axes(ax_init);
            plot(center_res(:, 2), center_res(:, 1), '.g', 'markersize', 10);
        end
        neuron_init_res = neuron.copy();
    end

    %Updates (P-B,P-S,P-T,P-S,P-T)
    neuron.update_background_parallel(use_parallel);
    neuron.update_spatial_parallel(use_parallel);       neuron.update_temporal_parallel(use_parallel);

    %Evaluate neurons after looking at their spatial and temporal footprints
    ids = [];       neuron.viewNeurons(ids, neuron.C_raw);      neuron.orderROIs('snr');

    %Increase the quality of fitting
    %might need to modify constrained_oasisAR1
    %please do not update_temporal again after this step
    options = neuron.options.deconv_options;
    for i=1:size(neuron.C,1)
        [c, s, option] = deconvolveCa(neuron.C_raw(i,:), options, 'ar1');
        pars = [option.pars, 0];
        [c2, s2, option] = deconvolveCa(neuron.C_raw(i,:), options, 'pars', pars);
        neuron.C(i,:) = c2;
    end

    %Updates (P-B,P-S,P-T,P-S,P-T)
    neuron.update_background_parallel(use_parallel);
    neuron.update_spatial_parallel(use_parallel);

    %Evaluate neurons after looking at their spatial and temporal footprints
    ids = [];       neuron.viewNeurons(ids, neuron.C_raw);      neuron.orderROIs('snr');

    %save the workspace after sanity check
    filename = ['CNMFe_',datestr(now, 'yymmdd_HHMM'),'_refined_',movie_name,'.mat'];
    save(filename);
end

%% Extract relevant quantities

A = neuron.reshape(neuron.A, 2); A = permute(A, [3 1 2]);   %Spatial footprints: cell# x pixels x pixels
%fr_area = sum(A,[2 3])./sum(A,'all');                       %Fractional area for each Spatial footprint
fr_area = squeeze(sum(sum(A,2),3))./sum(A(:));
FC_raw = neuron.C_raw.*fr_area;                             %Raw fluorescence (normalized)
FC = neuron.C.*fr_area;                                     %Denoised fluorescence (normalized)
noise_F = (FC - FC_raw);                                    %Noise temporale traces
thr = median(FC_raw,2)+7*std(noise_F,[],2);                 %Thresholds for Ca-transient detection
time = [0:1/metadata.cnmfe.Fs:(size(FC_raw,2)-1)/metadata.cnmfe.Fs];                      %Time trace

%% Plot fluorescence traces and spatial footprints (again)
until_neuron = size(FC_raw,1);              %number of neurons to be inspected size(C_raw,1) if you want to see all neurons
figure();   set(gcf, 'units','normalized','outerposition',[0.5 0 0.5 1]);
mx = max(max(FC_raw))/3;
for i = 1:until_neuron
    plot(time,FC_raw(i,:)-((i-1)*mx));    hold on;
    ylim([-until_neuron*mx mx*3]); xlabel('Time(s)');

end
hold off;

figure();   set(gcf, 'units','normalized','outerposition',[0 0.4 0.5 0.6]);
plot_contours(neuron.A,Cn,0.75,true,until_neuron);   %plot_contours(Aor,Cn,thr,display_numbers,max_number,Coor, ln_wd, ln_col)
figure();   set(gcf, 'units','normalized','outerposition',[0.1 0.05 0.3 0.4]);
r = corrcoef(FC_raw');   imagesc(r); colorbar;

%% Show spatial footprints and thresholds on activity traces
mkdir('ROI_Images/images');
figure();   set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
for i = 1:size(neuron.C_raw,1)
    subplot(3,2,1); neuron.image(neuron.A(:,i)); colorbar;  title(['ROI' num2str(i)]);
    subplot(3,2,2); surf(reshape(A(i,:,:),size(A,2),size(A,3)));    zlim([0 max(max(max(A)))]);
    subplot(3,2,[3 4]); plot(FC_raw(i,:));        hold on;    refline(0,thr(i));  hold off;     %ylim([-0.1 1]);
    subplot(3,2,[5 6]); plot(noise_F(i,:));       hold on;    refline(0,thr(i));  hold off;     %ylim([-0.01 0.1]);
  fname_sv = ['ROI_Images/images/ROI_',num2str(i)];
    print(fname_sv,'-dpng')
    %  w = waitforbuttonpress;
  pause(0.1)
end

close all
% Save data:

%% save results
results = neuron.obj2struct();
results.metadata = metadata;
save('results.mat','results');

ImBat_ROIoverlay(results);
print('ROI_Images/overlay','-dpng')



