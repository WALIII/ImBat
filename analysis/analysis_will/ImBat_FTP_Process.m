function ImBat_FTP_Process(local_dir)
% ImBat_FTP_Process.m

% Server to VM function. Use FTP to transfer files locally for processing,
% then upload the finished folder when complete

% WAL3
% 12/11/2020

local_dir = '/Users/ARGO/Desktop/TEST';
Folders_per_session = 1; %how many folders to load in per session




%log on to server:
ftpobj = ftp('169.229.54.11','liberti','batsFly!19');
cd(ftpobj, 'server_home/users/tobias/flight/data_processed/ready2analyze_Gbats/Gi200408/extracted/Gio_200408_rest1-1_extraction');

% list contents of folder, take only relevant folders
filelist = dir(ftpobj);

% Only up/download in one at a time:
itter2process = find(mod(size(filelist,1),Folders_per_session) == 0);
itter2process = % add last one.. size(filelist,1)

for i = 1:size(filelist,1)
    
    % Pick Folder
    folder2use = filelist(i).name
    
    % Only up/download in one at a time:
    mget(ftpobj,folder2use,local_dir)
    
    % Process folder
    if ismember(i,itter2process)
        disp('processing data...')
        %     ImBat_BatchExtract
        
        % find all folders,
        files = dir(local_dir);
        files(ismember( {files.name}, {'.', '..'})) = [];  %remove . and .. and Processed
        
        for ii = 1:size(files,1)
            % upload the corrected folder
            mput(ftpobj,'Screen Shot 2020-12-05 at 6.58.10 PM.png');
        end
    end
    
    
    
    
end

close(ftpobj)