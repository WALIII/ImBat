function ImBat_START(varargin)
% ImBat_START

% Main wrapper for all ImBat functions. Will Motion correct, and basic ROI
% extraction using CNMFe

% WAL3
% d05/16/2019


ROI_flag = 0; % run ROI extraction
Analysis_flag = 0; % run basic ROI analysis...
extract  = 1;
reExtract =0;

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
if exist([flight_subFolders(ii).folder,'/',flight_subFolders(ii).name,'/','processed','/','Alignment.mat'])>0;
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

    
    
    
    
    