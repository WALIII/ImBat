function ImBat_START(varargin)
% ImBat_START

% Main wrapper for all ImBat functions. Will Motion correct, and basic ROI
% extraction using CNMFe

% WAL3
% d05/16/2019
% d02/24/2020


% Core Functions:

% ImBat_TemporalDownSample.m  temporal downsampleing...

% FS_AV_Parse_batch:
%           Will extract .mov files to .tif and aligned .mat metadata

% ImBat_processVideos:
%           contains 'ImBat_denoise.m' For de-noising video ( artifact rejection)
%                    'CNMFe_align.m' For motion correction ( rigid/non-rigid)

% CNMFe_extract: CNMF extraction


% Default extraction Params
ROI_flag = 1; % run ROI extraction
ROI_flag_reset = 1;
reROI_extract = 0; %rerun ROI extraction

Analysis_flag = 1; % run basic ROI analysis...
extract  = 1; % extract the basic ROI timeseries
reExtract = 1; % re-extract, in the event that things have been extracted already.
extract_track = 1;
Mov_extract_flag = 0; % run .mov extraction ( Mac only...);
Bat_Cluster =1; % extract the bat_cluster tracking
alignment = 1;


% Default Movie Paramaters:
metadata.temporal_downsample = 1; % temporal downsampleing
metadata.spatial_downsample = 0.4; % spatial downsampling
metadata.median_filter_kernal = 3; % median filtering
metadata.artifact_reject = 1; % median filtering
metadata.initial_median_filter_kernal = 11;

% Motion ocrrection:
metadata.moco.itter = 1;
metadata.moco.bin_width = 1000;

% Default CNMFe Paramaters:
metadata.cnmfe.min_corr = 0.9;     % minimum local correlation for a seeding pixel
metadata.cnmfe.min_pnr = 4.0;       % minimum peak-to-noise ratio for a seeding pixel
metadata.cnmfe.tsub = 5;           % temporal downsampling factor
metadata.cnmfe.ssub = 1;           % temporal downsampling factor
metadata.cnmfe.gSig = 4;           % pixel, gaussian width of a gaussian kernel for filtering the data. 0 means no filtering
metadata.cnmfe.gSiz = 4*metadata.cnmfe.gSig+1;    % pixel, approximate neuron diameter
metadata.cnmfe.ssub = 1;   % 

% Mergting thesholds:
%metadata.cnmfe.min_pnr = 10;       % minimum peak-to-noise ratio for a seeding pixel



processed_FN = ['processed_',datestr(now,'yyyy_mm_dd__hhMM')];

% Housekeeping:
metadata.Extraction_date = datestr(now);
metadata.processed_FN = processed_FN;


% Manual inputs
vin=varargin;
for i=1:length(vin)
    if isequal(vin{i},'roi') % manually inputing a sort order
        ROI_flag=vin{i+1};
        ROI_flag_reset = ROI_flag;
    elseif isequal(vin{i},'place');
        analysis_flag = vin{i+1};
    elseif isequal(vin{i},'metadata');  % pass along metadata file if need be...
        metadata = vin{i+1};
    elseif isequal(vin{i},'rextract');
        reExtract=vin{i+1};
    elseif isequal(vin{i},'realign'); % Re-alignment-dont re-process everything
        % Default extraction Params
        ROI_flag = 0; % run ROI extraction
        ROI_flag_reset = 0;
        reROI_extract = 0; %rerun ROI extraction
        extract  = 0; % extract the basic ROI timeseries
        extract_track = 1;
        alignment = 1;
        Mov_extract_flag = 0; % run .mov extraction ( Mac only...);
        Bat_Cluster =1; % extract the bat_cluster tracking
        
    end
end


if Bat_Cluster ==1;
    metadata.tracking_file_type =  'Bat_Cluster';
else
    metadata.tracking_file_type = '_track';
end

if Mov_extract_flag  ==1;
    disp(' %===========[ EXTRACTING  .mov  FILES  ] ===========% ');
    FS_AV_Parse_batch;
end




% Get all folders in directory
files = dir(pwd);
files(ismember( {files.name}, {'.', '..','Processed'})) = [];  %remove . and .. and Processed

% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags);
% Print folder names to command window.
for k = 1 : length(subFolders)
    fprintf('Sub folder #%d = %s\n', k, subFolders(k).name);
end


%% Extraction, Motion Correction, ROI extraction


% For each folder ( extraction, and motion correction step)
for i = 1:length(subFolders);
    disp(['entering folder', char(subFolders(i).name)])
    
    cd([subFolders(i).folder,'/',subFolders(i).name]);
    
    % index every subfolder...
    flight_files = dir(pwd);
    flight_files(ismember( {flight_files.name}, {'.', '..','Processed'})) = [];  %remove . and .. and Processed
    flight_dirFlags = [flight_files.isdir];% Get a logical vector that tells which is a directory.
    flight_subFolders = flight_files(flight_dirFlags);% Extract only those that are directories.
    
    
    for ii = 1:length(flight_subFolders);
        cd([subFolders(i).folder,'/',subFolders(i).name]);
        
        % load tracking data
        if extract_track ==1;
            track_fname = flight_subFolders(ii).name;
            track_fname = extractBefore( track_fname,'_extraction');
            if strcmp(metadata.tracking_file_type, 'Bat_Cluster')
                try
                    track_fname = [track_fname,'-Bat_Cluster_track'];
                    load(track_fname);
                catch % TO DO: eliminate redundancy ( turn into 'exist' not try/catch...)
                    track_fname = flight_subFolders(ii).name;
                    track_fname = extractBefore( track_fname,'_extraction');
                    track_fname = [track_fname,'_track'];
                    load(track_fname);
                    metadata.tracking_file_type = 'track';
                end
            else
                track_fname = [track_fname,'_track'];
                load(track_fname);
                metadata.tracking_file_type = 'track';
            end
        end
        
        
        %%====[ Motion Correction ]======%%
        
        cd([flight_subFolders(ii).folder,'/',flight_subFolders(ii).name])% index into the flight_subfolder
        if extract ==1;
            % Run processing script
            mkdir(processed_FN);
            ImBat_processVideos('metadata',metadata);
            disp('processing!!');
        else
            
            % get most recent processed folder:
            temp_files = dir(pwd);
            temp_files(ismember( {temp_files.name}, {'.', '..','Processed'})) = [];  %remove . and .. and Processed
            temp_dirFlags = [temp_files.isdir];% Get a logical vector that tells which is a directory.
            temp_subFolders = temp_files(temp_dirFlags);% Extract only those that are directories.
            processed_FN = temp_subFolders(size(temp_subFolders,1)).name;
        end
        
        
        cd(processed_FN) % move to processed folder...
        
        % Check if roi extraction folder exists
        if exist([flight_subFolders(ii).folder,'/',flight_subFolders(ii).name,'/',processed_FN,'/','Motion_corrected_Data_DS_neurons'])>0;
            disp('ROIs already extracted..');
            
            if reROI_extract ==1
                disp('Re-Extracting ROIs...');
                ROI_flag = 1 ;
                
            else
                disp('Moving to the next folder...');
            end
        end
        
        
        %%====[ CNMF-e ROI Extraction ]======%%
        if ROI_flag ==1;
            disp('extracting ROIs...')
            nam = './Motion_corrected_Data.mat'
            
            CNMFe_extract2(nam,'metadata',metadata);
            
            
            ROI_flag = ROI_flag_reset;
        end 
            %%====[ Aligning Time Stamps ]======%%
            % to do: add logfile
            if alignment ==1;
                disp('Aligning Timestamps...');
                load('AV_data.mat');
                [out] = ImBat_alignTimeStamps(audio,video,AnalogSignals,Markers);
            end
        
        clear video audio Markers AnalogSignals out % Clear vars from RAM
    end
    
end
