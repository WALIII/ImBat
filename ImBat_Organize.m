function ImBat_Organize


DIR = pwd;
% Get all folders in directtory
files = dir(DIR);
files(ismember( {files.name}, {'.', '..','Processed'})) = [];  %remove . and .. and Processed

% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags);


for i = 1: size(subFolders,1)
    
%
folder_name = [subFolders(i).name,'/tracking/Generated_C3D_files'];  
destination_name = [subFolders(i).name,'/scope/mat'];

% Check if extracted
if exist('processed')>1
    disp(' Tracking already extracted')
else
% cd into subfolder, get tracking data
  disp('Extracting Tracking Data...')
extract_tracking_data(folder_name);
end

% put in mat folder
cd(processed);

% Move mat folders to new location...
mov_listing=dir(fullfile(pwd,'*.mat'));
matfilenames={mov_listing(:).name};
% Move files...
for ii = 1:size(matfilenames,4);
copyfile(matfilenames(ii),destination_name);
end
end


% Make a Processed Folder with all extracted data...




