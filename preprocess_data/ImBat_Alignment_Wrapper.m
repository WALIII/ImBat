function ImBat_Alignment_Wrapper
% Re-run alignmnet on all files for robustness

% Temp function to re-do alignment.. 


% d10/09/2020
% WAL3




HomeDir = cd;
% Get all folders in directory
files = dir(pwd);
files(ismember( {files.name}, {'.', '..','Processed','plots','new3'})) = [];  %remove . and .. and Processed

% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags);
% Print folder names to command window.
for k = 1 : length(subFolders)
    fprintf('Sub folder #%d = %s\n', k, subFolders(k).name);
end

for i = 1:length(subFolders);
    % index into folders,
    
    cd([subFolders(i).folder,'/',subFolders(i).name]);
    

% Now, the meat:
ImBat_START('realign');
end
cd(HomeDir);
