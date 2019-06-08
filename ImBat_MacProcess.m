
% For all subfolders in folder


files = dir(pwd);
DIR = pwd

files(ismember( {files.name}, {'.', '..','DATA'})) = [];  %remove . and .. and Processed


% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags);
% Print folder names to command window.
for k = 1 : length(subFolders)
	fprintf('Sub folder #%d = %s\n', k, subFolders(k).name);
end

for i = 1:length(subFolders);
    try
    cd([subFolders(i).name,'scope');
    FS_AV_parse_batch
    catch
       disp([' no scope subfolder for ',subFolders(i).name]);
    end
    cd(DIR)
end

% run FS_AV_parse_batch

% index into the next folder

