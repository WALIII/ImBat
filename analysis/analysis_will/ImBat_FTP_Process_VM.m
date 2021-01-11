function ImBat_FTP_Process_VM
% ImBat_FTP_Process_VM.m

% Server--> VM data transfer. Use FTP to transfer files locally for processing,
% then upload the finished folder when complete. Send text on error/completetion .

% WAL3
% 01/09/2020

% local dir should be current directory.. % example:
local_dir = cd;
% Or, hardcode it:
% example: local_dir = '/Users/ARGO/Desktop/TEST';
% example:  local_dir ='C:\Users\Tobias\Documents\DATA_2'; % for rosetta

server_dir = 'server_home/users/tobias/flight/data_processed/ready2analyze_Gbats/Processed';

Folders_per_session = 1; %how many folders to load in per session ( fewer on VM)
[ret, computer_name] = system('hostname'); % what computer are we using
computer_name = convertCharsToStrings(computer_name);
computer_name = strtrim(computer_name);

%log on to server:
ftpobj = ftp('169.229.54.11','liberti','batsFly!19','System','WINDOWS');
cd(ftpobj, server_dir);


% list contents of folder, take only relevant folders
filelist = dir(ftpobj);
dirFlags = [filelist.isdir];
filelist= filelist(dirFlags);

% Only up/download in one at a time:
itter2process = find(mod(1:size(filelist,1),Folders_per_session) == 0);
itter2process = [itter2process, size(filelist,1)];

try
    
    for i = 1:size(filelist,1)
        % Pick Folder, get .mat files
         folder2use = [filelist(i).name,'/extracted/*.mat'];        

        % Only up/download in one at a time. 
            % 1. Get tack tiles
       manifest = mget(ftpobj,folder2use,local_dir);
       manifest = dir([filelist(i).name,'/extracted/*.mat']);
              manifest = {manifest.name};

     % 2. Get tifs and mats, based on track files
     for ii = 1:length(manifest);
         folder2use2 = manifest{ii};
       %  idx2use = findstr(folder2use2,'\');
        % if length(idx2use)==0
        %      idx2use = findstr(folder2use2,'/');
        % end
         folder2use2 = folder2use2(1:end-10);
         folder2use2 = [filelist(i).name,'/extracted/',folder2use2, '_extraction'];
   
         % get mat and .tif files
       %  ind2rplce = findstr(folder2use2,'\');
         %folder2use2(ind2rplce) = '/';
           mget(ftpobj,[folder2use2,'/*.mat'],local_dir); % mat files
           mget(ftpobj,[folder2use2,'/*.tif'],local_dir);  % tif files
           clear folder2use2 idx2use
     end
     clear folder2use
        
        % Process folder
        if ismember(i,itter2process)
            disp('processing data...')
            
            % Process Data:
            ImBat_BatchExtract
            
            % find all session folders ( 2 = session)
            files2 = dir(local_dir);
            files2(ismember( {files2.name}, {'.', '..','Processed','error'})) = [];  %remove . and .. and Processed
            dirFlags2 = [files2.isdir];
            subFolders2 = files2(dirFlags2);  % Extract only those that are directories.
            
            % upload just the the corrected folder
            for ii = 1:size(subFolders2)% for all sessions..
                % find all extractions in the extracted folder ( i.e. fly, rest, etc)
                files3 = dir([subFolders2(ii).name,'/extracted']);
                files3(ismember( {files3.name}, {'.', '..','Processed','error'})) = [];  %remove . and .. and Processed
                dirFlags3 = [files3.isdir];            % Get a logical vector that tells which is a directory.
                subFolders3 = files3(dirFlags3);            % Extract only those that are directories.
                
                for iii = 1:size(subFolders3)%
                    % FIND just the most recent extracttion in each folder ( fly-1/processed_date, rest-1/processed_date, etc)
                    files4 = dir([subFolders3(iii).folder,'/', subFolders3(iii).name]);
                    files4(ismember( {files4.name}, {'.', '..','Processed','error'})) = [];  %remove . and .. and Processed
                    dirFlags4 = [files4.isdir];            % Get a logical vector that tells which is a directory.
                    subFolders4 = files4(dirFlags4);            % Extract only those that are directories.
                    % FIND the last entry- this is the most recet folder:
                    P2use_name = subFolders4(size(subFolders4 ,1)).name;
                    P2use_folder = subFolders4(size(subFolders4 ,1)).folder;
                    % This is the processed folder path to upload:
                    to_upload = [P2use_folder,'/',P2use_name];
                    % upload this subfolder via FTP:
                    % TO DO: cd into the right place on the server first...,
                    temp_path = [subFolders2(ii).name,'/extracted/',subFolders3(iii).name];
                    cd(ftpobj, temp_path);
                    mput(ftpobj,to_upload);
                    cd(ftpobj, '../../../');       

                end
                        fclose('all')
                        rmdir(subFolders2(ii).name,'s')
                   
            end
        end
    end
    
catch
    % ERROR in processing- send a text
    
    message2send = ['Error processing: ',filelist(i).name];
    send_text_message2('617-529-0762','Verizon',computer_name,message2send)
    
end

% Processing complete- send a text

send_text_message2('617-529-0762','Verizon',computer_name,'Processing Complete!')
close(ftpobj)