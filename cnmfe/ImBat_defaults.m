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
metadata.cnmfe.gSig = 4;           % pixel, gaussian width of a gaussian kernel for filtering the data. 0 means no filtering
metadata.cnmfe.gSiz = 4*metadata.cnmfe.gSig+1;    % pixel, approximate neuron diameter
metadata.cnmfe.ssub = 1;   % 
center_psf = 1;              % set the value as true when the background fluctuation is large (usually 1p data)


%defualt alignment params:
metadata.moco.itter = 5;
metadata.moco.bin_width = 200;

metadata.processed_FN = 'temp';
mkdir(metadata.processed_FN);