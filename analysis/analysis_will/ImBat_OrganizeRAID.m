function [log] = ImBat_OrganizeRAID()
% organize the data on the RAID to make sure that each directory has the
% most recent extraction on it.

% TO DO: update Locations:


try
    comp_str = '/Volumes/';
    cd(comp_str);
    
catch
    disp('select path to RAID');
selpath = uigetdir;
end_str = strfind(selpath,'server_home');
comp_str = selpath(1:end_str-1);
end

counter = 1;
% 1. ZBats: /Volumes/server_home/users/tobias/flight/data_processed/ready2analyze_Zbats
Locations2Use{1} = [comp_str, 'server_home/users/tobias/flight/data_processed/ready2analyze_Zbats'];

% 2. GBats: /Volumes/server_home/users/tobias/flight/data_processed/ready2analyze_Gbats
Locations2Use{2} = [comp_str, 'server_home/users/tobias/flight/data_processed/ready2analyze_Gbats'];


% 3. /Volumes/server_home/users/tobias/flight/data_processed/topQualityData/analysis_done
Locations2Use{3} = [comp_str, 'server_home/users/tobias/flight/data_processed/topQualityData/analysis_done/Gal'];
Locations2Use{4} = [comp_str, 'server_home/users/tobias/flight/data_processed/topQualityData/analysis_done/Ge'];
Locations2Use{5} = [comp_str, 'server_home/users/tobias/flight/data_processed/topQualityData/analysis_done/Ge_2'];
Locations2Use{6} = [comp_str, 'server_home/users/tobias/flight/data_processed/topQualityData/analysis_done/Gi'];
Locations2Use{7} = [comp_str, 'server_home/users/tobias/flight/data_processed/topQualityData/analysis_done/Go'];
Locations2Use{8} = [comp_str, 'server_home/users/tobias/flight/data_processed/topQualityData/analysis_done/z2'];
Locations2Use{9} = [comp_str, 'server_home/users/tobias/flight/data_processed/topQualityData/analysis_done/za'];
Locations2Use{10} = [comp_str, 'server_home/users/tobias/flight/data_processed/topQualityData/analysis_done/zu'];



% 4. Tobias's subfolders:
% Locations2Use{6} = [comp_str, 'server_home/users/tobias/flight/data_processed/topQualityData/analysis_done/new3/z2'];
% Locations2Use{7} = [comp_str, 'server_home/users/tobias/flight/data_processed/topQualityData/analysis_done/new3/Gen'];
% Locations2Use{8} = [comp_str, 'server_home/users/tobias/flight/data_processed/topQualityData/analysis_done/new3/Gal'];
%



% Run through Locations, and get the right folders:
for i = 1:length(Locations2Use);
    for ii = 1:length(Locations2Use);
        if i == ii;
            disp('do not compare the same folder, skipping...');
        else
            % get folder manifest for each location:
            %% Folder Manifest:
            % Folder 1;
            files = dir(Locations2Use{i});
            files(ismember( {files.name}, {'.', '..','Processed'})) = [];  %remove . and .. and Processed
            % Get a logical vector that tells which is a directory.
            dirFlags = [files.isdir];
            % Extract only those that are directories.
            subFolders = files(dirFlags);
            % Get only data Folders:
            subFolders(find(cellfun(@length, {subFolders.name})<8)) = [];
            FoldMan_1 = subFolders;
            
            % Folder 2
            files2 = dir(Locations2Use{ii});
            files2(ismember( {files2.name}, {'.', '..','Processed'})) = [];  %remove . and .. and Processed
            % Get a logical vector that tells which is a directory.
            dirFlags2 = [files2.isdir];
            % Extract only those that are directories.
            subFolders2 = files2(dirFlags2);
            % Get only data Folders:
            subFolders2(find(cellfun(@length, {subFolders2.name})<8)) = [];
            FoldMan_2 = subFolders2;
            
            %% Compare Manifests:
            for iii = 1:length(FoldMan_1)
                idx = ismember( {FoldMan_2.name}, {FoldMan_1(iii).name});
                a = find(idx==1);
                if a>1;
                    disp('!');
                    % check for files to transfer:
                    %                     Dir1 = FoldMan_2(a).name;
                    %                     Dir2 = FoldMan_1(iii).name;
                    % check all processing subfolders:
                    files3 = dir([FoldMan_2(a).folder,'/',FoldMan_2(a).name,'/extracted/']);
                    files3(ismember( {files3.name}, {'.', '..','Processed'})) = [];
                    dirFlags3 = [files3.isdir];
                    subFolders3 = files3(dirFlags3);
                    
                    files3b = dir([FoldMan_1(iii).folder,'/',FoldMan_1(iii).name,'/extracted/']);
                    files3b(ismember( {files3b.name}, {'.', '..','Processed'})) = [];
                    dirFlags3b = [files3b.isdir];
                    subFolders3b = files3b(dirFlags3b);
                    
                    for iv = 1:length(subFolders3);
                        % Compare Folders
                        
                        files4 = dir([subFolders3(iv).folder,'/',subFolders3(iv).name]);
                        files4(ismember( {files4.name}, {'.', '..','Processed'})) = [];
                        dirFlags4 = [files4.isdir];
                        subFolders4 = files4(dirFlags4);
                        X = subFolders4(contains({subFolders4.name},'processed_'));
                        X_toCompare_01 = X(length(X));
                        
                        files4b = dir([subFolders3b(iv).folder,'/',subFolders3b(iv).name]);
                        files4b(ismember( {files4b.name}, {'.', '..','Processed'})) = [];
                        dirFlags4b = [files4b.isdir];
                        subFolders4b = files4b(dirFlags4b);
                        Xb = subFolders4b(contains({subFolders4b.name},'processed_'));
                        X_toCompare_02 = Xb(length(Xb));
                        
                        % NOW copy the bigger folder
                        if X_toCompare_01.datenum>X_toCompare_02.datenum;
                            disp(['coping', X_toCompare_01.folder, ' : ', X_toCompare_01.name ' to: ']);
                            disp([' `------>', X_toCompare_02.folder]);
                            copyfile([X_toCompare_01.folder,'/',X_toCompare_01.name],[X_toCompare_02.folder,'/',X_toCompare_01.name]);
                            
                        elseif X_toCompare_02.datenum>X_toCompare_01.datenum;
                            disp(['coping', X_toCompare_02.folder, ' : ', X_toCompare_02.name ' to: ']);
                            disp([' `------>', X_toCompare_01.folder]);           
                            copyfile([X_toCompare_02.folder,'/',X_toCompare_02.name],[X_toCompare_01.folder,'/',X_toCompare_02.name]);
                            % test:
                            % copyfile([X_toCompare_02.folder,'/',X_toCompare_02.name], ['/Users/ARGO/Desktop','/',X_toCompare_02.name);
                            
                            %elseif
                            % report date num for log
                            
                        end
                    end
                    
                    
                end
                
            end
        end
    end
end
    
