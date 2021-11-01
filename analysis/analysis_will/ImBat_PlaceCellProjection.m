
% NOTE: data can be alos found:
% Load data from Processed/AllDataHeat/matlab.m

% Load data

files = dir(pwd);
files(ismember( {files.name}, {'.', '..','Processed'})) = [];  %remove . and .. and Processed

% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags);
% Print folder names to command window.

figure();
counter = 1;

for ii = 1 : length(subFolders)
        load([subFolders(ii).name, '/CellReg_files/ROI_Data/Saved_Data/Aligned_Data.mat'])
[output{ii}] = ImBat_2D_overdays(CombinedROI,flightPaths,FlightAlignedROI);
clear CombinedROI flightPaths FlightAlignedROI
end

% Load data from Processed/AllDataHeat/matlab.m

for i = 1:8;
    G(:,:,i) = output{i}.M_Unclustered;
end
figure(); imagesc(mean(G,3));