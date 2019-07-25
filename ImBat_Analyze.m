function ImBat_Analyze

% Main wrapper for all analyzing bat calcium imaging and flight data.
% Goals:
% 1. Sanity checks: max projection w/ ROI overlay, flight + traces, flights in 3d
% 2. Quantify behavior: #total flights, #rewarded flights, variability of
% flights (covariance matrix), stability of flights over days (Karthik code)
% 3. Cells aligned to behavior:
%    a. Place cell plots: heat maps, plot3, schnitzer plots
%    b. Preflight activity: 2 sec, 20 sec before
%    c. Postflight activity: 2 sec, 20 sec, rewarded vs nonrewarded

% TAS
% d07/24/2019

% Manual inputs
analysisFlag = 1;
reAnalyze = 1;
sanityCheckFlag = 1;

% Get all folders in directory
files = dir(pwd);
files(ismember( {files.name}, {'.', '..','.DS_Store','Thumbs.db'})) = [];  %remove . and .. and DS_Store

% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags);
% Print folder names to command window.
for k = 1 : length(subFolders)
    fprintf('Date folder #%d = %s\n', k, subFolders(k).name);
end

%% Perform analysis on each folder
for i = 1:length(subFolders)
    disp(['entering folder ', char(subFolders(i).name)])
    cd([subFolders(i).folder,'/',subFolders(i).name]);
    
    % index every subfolder...
    trackFiles = dir('*track.mat');
    imageFolders = dir('*_extraction');
    %print image data folder names
    for kk = 1 : length(imageFolders)
        fprintf('Image data folder #%d = %s\n', kk, imageFolders(kk).name);
        %check that track data matches image data
        if strcmp(imageFolders(kk).name(1:end-10),trackFiles(kk).name(1:end-9)) == 0
            fprintf('Tracking and image data do not match');
            analysisFlag = 0;
        end
        
        % Check if folder exists
        if exist([imageFolders(kk).folder,'/',imageFolders(kk).name,'/','analysis'])>0;
            disp('Folder already analyzed..');
            if reAnalyze ==1
                disp('Re-Analyzing...');
            else
                disp('Moving to the next folder...');
                analysisFlag = 0 ;
            end
        end
        
        % load tracking and cell data
        if analysisFlag == 1
            trackData = load([trackFiles(kk).folder,'/',trackFiles(kk).name]);
            cellData = load([imageFolders(kk).folder,'/',imageFolders(kk).name,'/','processed','/','Motion_corrected_Data_DS_results.mat']);
            videoData = load([imageFolders(kk).folder,'/',imageFolders(kk).name,'/','processed','/','Motion_corrected_Data_DS.mat']);
            
            % make new analysis directory for .fig and .tif files
            cd([imageFolders(kk).folder,'/',imageFolders(kk).name]);
            mkdir('analysis');
            disp('Analyzing!!');
        end
        
        %perform basic sanity checks on imaging and flight data
        if sanityCheckFlag == 1
            % max projection from each session without ROI overlay
            [Ymax, Y] = ImBat_Dff(videoData.Y);
            % max projection with ROI overlay
            
        end
        
    end
    
    
    
end