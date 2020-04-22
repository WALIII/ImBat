function [metadata] = ImBat_defaults;
% creates 'default' metadata file for troubleshooting pipelines:


% Default Movie Paramaters:
metadata.temporal_downsample = 5; % temporal downsampleing
metadata.spatial_downsample = 0.4; % spatial downsampling
metadata.median_filter_kernal = 3; % median filtering
metadata.artifact_reject = 1; % median filtering
metadata.initial_median_filter_kernal = 11;
% Default CNMFe Paramaters:
metadata.cnmfe.min_corr = 0.9;     % minimum local correlation for a seeding pixel
metadata.cnmfe.min_pnr = 50;       % minimum peak-to-noise ratio for a seeding pixel
% defualt alignment params:
metadata.moco.itter = 5;
metadata.moco.bin_width = 200;

metadata.processed_FN = 'temp';
mkdir(metadata.processed_FN);