function ImBat_BatchExtract
% extract c3D files and ,mov files, then run ImBat_Start ( within a single
% script)

%

% 
% Remove competing code from path:
try
 addpath(genpath('C:\Users\WAL3\Documents\MATLAB\CNMF_E'));
rmpath(genpath('C:\Users\WAL3\Documents\MATLAB\CaImAn-MATLAB'));
catch
    disp(' not using Wills computer..');
end
parseFiles = 1; % MAC ONLY
HomeDir = cd;
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


for i = 1:length(subFolders);
    % index into folders,
    
    cd([subFolders(i).folder,'/',subFolders(i).name]);
    
    % check if any mov has been extracted,if not,
    if isfolder('extracted')==1
        disp('.mov files already extracted...');
        % move all data into 'extracted folder
        
    elseif parseFiles ==1; % mac only 
            %   extract .mov files:
         %   FS_AV_Parse_batch(pwd,'mat_dir','/extracted')
            %   extract c3d files:
             mkdir('extracted');

            processed_dir = [pwd,'/extracted/'];
            
            c3dList = dir([pwd filesep '*.c3d']);
            
            %run through list of c3d files in that directory, convert to mat, and save
            %to processed directory
            for i = 1:length(c3dList)
                fileName = extractBefore(c3dList(i).name,'-Bat_Cluster.c3d'); %Bat
                %fileName = extractBefore(c3dList(i).name,'-Bat_Cluster.c3d'); %Bat
                dateSesh = datestr(datetime(c3dList(i).date),'yymmdd');
                batName = extractBefore(fileName,['_' dateSesh]);
                %sessionNum = fileName(end);
                
                
                [Markers,VideoFrameRate,AnalogSignals,AnalogFrameRate,Event,ParameterGroup,CameraInfo,ResidualError]=readC3D_analog([cd filesep c3dList(i).name]); %convert file
                %plot ttl impulses to check they are linear and not missing ttl
                event_ttls = AnalogSignals(:,2);
                [R,LT,UT,LL,UL] = risetime(event_ttls,VideoFrameRate);
                
                %save new mat file in both original directory and copied directory for processing
                save([processed_dir fileName '_track' '.mat'],'AnalogFrameRate','AnalogSignals','Markers','VideoFrameRate');
                %save([copy_dir fileName '_track' '.mat'],'AnalogFrameRate','AnalogSignals','Markers','VideoFrameRate');
            end
            
        else % move files into extracted folder ( also see: ImBat_ReformatExtracted)
            mkdir('extracted');
            A= [subFolders(i).folder,'/',subFolders(i).name]
            B= [A,'/extracted']
            f_raw=dir(fullfile(A));
            for i = 3:size(f_raw,1)-1
                movefile(fullfile(f_raw(i).folder,f_raw(i).name),B);
            end
        end
    
        
        % Run extraction
        close all;
        ImBat_START;
        
    end
    
    cd(HomeDir);
    
    
    
    
