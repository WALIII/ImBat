function ImBat_ReformatExtracted
% ImBat_ReformatExtracted.m

% hopefully, one time use function to reformat the processed data folders

% WAL3
% d12/12/2020

files = dir(pwd);
files(ismember( {files.name}, {'.', '..','Processed','error','plots','new3'})) = [];  %remove . and .. and Processed

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
    
    % check if any mov has been extracted,if not,
    if isfolder('extracted')==1
        
    else % move files into extracted folder
        disp('old format detected, reformatting...');
        mkdir('extracted');
        A= [subFolders(i).folder,'/',subFolders(i).name]
        B= [A,'/extracted']
        f_raw=dir(fullfile(A));
        for i = 3:size(f_raw,1)-1
            movefile(fullfile(f_raw(i).folder,f_raw(i).name),B);
        end
    end
end
