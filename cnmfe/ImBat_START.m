function ImBat_START(varargin)
% ImBat_START

% Main wrapper for all ImBat functions. Will Motion correct, and basic ROI
% extraction using CNMFe

% WAL3
% d05/16/2019


ROI_flag = 1; % run ROI extraction
<<<<<<< Updated upstream
reROI_extract = 0; %rerun ROI extraction
Analysis_flag = 0; % run basic ROI analysis...
extract  = 1;
reExtract =1;
=======
reROI_extract = 1; %rerun ROI extraction
Analysis_flag = 1; % run basic ROI analysis...
extract  = 1; % extract the basic ROI timeseries
reExtract = 1; % re-extract, in the event that things have been extracted already.
ROI_flag_reset = 1;
extract_track = 0;
Mov_extract_flag = 0; % run .mov extraction ( Mac only...);

% Default Movie Paramaters:
metadata.temporal_downsample = 5; % temporal downsampleing
metadata.spatial_downsample = 0.4; % spatial downsampling
metadata.median_filter_kernal = 3; % median filtering
metadata.artifact_reject = 0; % median filtering
metadata.initial_median_filter_kernal = 11;
% Default CNMFe Paramaters:

processed_FN = ['processed_',datestr(now,'yyyy_mm_dd__hhMM')];


>>>>>>> Stashed changes

% Manual inputs
vin=varargin;
for i=1:length(vin)
  if isequal(vin{i},'ROI') % manually inputing a sort order
    ROI_flag=vin{i+1};
  elseif isequal(vin{i},'Place')
      Analysis_flag = vin{i+1};
  elseif isequal(vin{i},'reExtract');
      reExtract=vin{i+1};
end
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
% Get a logical vector that tells which is a directory.
flight_dirFlags = [flight_files.isdir];
% Extract only those that are directories.
flight_subFolders = flight_files(flight_dirFlags);

for ii = 1:length(flight_subFolders);
cd([subFolders(i).folder,'/',subFolders(i).name]);

% Check if folder exists 
if exist([flight_subFolders(ii).folder,'/',flight_subFolders(ii).name,'/','processed','/','Motion_corrected_Data_DS.mat'])>0;
    disp('Folder already extracted..');
    
    if reExtract ==1
        disp('Re-Extracting...');
        
    else
        disp('Moving to the next folder...');
        extract = 0 ;
    end
end
    

    % load tracking data 
    track_fname = flight_subFolders(ii).name;
    track_fname = extractBefore( track_fname,'_extraction');
    track_fname = [track_fname,'_track.mat'];
    load(track_fname);
    

% gindex into the flight_subfolder
cd([flight_subFolders(ii).folder,'/',flight_subFolders(ii).name])
if extract ==1;    
% Run processing script 
mkdir('processed');
ImBat_processVideos('downsample',0.5);
disp('processing!!');
end

cd('processed') % move to processed folder...

% Check if roi extraction folder exists 
if exist([flight_subFolders(ii).folder,'/',flight_subFolders(ii).name,'/','processed','/','Motio_corrected_Data_DS_neurons'])>0;
    disp('ROIs already extracted..');
    
    if reROI_extract ==1
        disp('Re-Extracting ROIs...');
        
    else
        disp('Moving to the next folder...');
        ROI_flag = 0 ;
    end
end


if ROI_flag ==1;
    disp('extracting ROIs...')
    nam = './Motion_corrected_Data_DS.mat'
    CNMFe_extract(nam);
end


disp('Aligning Timestamps...');
load('AV_data.mat');
[out] = ImBat_alignTimeStamps(audio,video,AnalogSignals,Markers);
% Clear vars from RAM
clear video audio Markers AnalogSignals out
extract =1;
end

end




% Check if there is a mat file

    % If no, run FS_AV_parse_batch
    
    % If yes, look into all the folders and check for processd files
        %If process files exist, break loop
        % if no, run ImBat_Process
 

%% ROI extraction

    % for each day folder:
        % enter each sub-folder, grab imaging data run ROI extraction
        
 
        
%% Align Timestamps:

    
    
    
    
    