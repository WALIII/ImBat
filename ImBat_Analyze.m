function ImBat_Analyze

global batName dateSesh sessionType topROI

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

%topROI is top% of cells you want to look at
topROI = 30;
% Manual inputs
analysisFlag = 1;
reAnalyze = 1;
%roi plot flags
plotROIFlag = 1;
centroidFlag = 1;
roiHeatFlag = 1;
binaryMaskFlag = 1;
%flight plot flags
plotFlightsFlag = 1;
flightPathsAllFlag = 1;
flightPathsFeederFlag = 1;

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
            alignment = load([imageFolders(kk).folder,'/',imageFolders(kk).name,'/','processed','/','Alignment.mat']);
            fileName = extractBefore(imageFolders(kk).name,'_extraction'); %get filename for titles
            dateSesh = imageFolders(kk).folder(end-5:end);
            batName = extractBefore(fileName,['_' dateSesh]);
            sessionType = extractAfter(fileName,[dateSesh '_']);
            % make new analysis directory for .fig and .tif files
            cd([imageFolders(kk).folder,'/',imageFolders(kk).name]);
            mkdir('analysis');
            disp('Analyzing!!');
        end
        
        %%perform basic sanity checks on imaging and flight data
        if plotROIFlag == 1
            % max projection from each session without ROI overlay
            [Ymax, Y, maxFig] = ImBat_Dff(videoData.Y);
            hold on
            %save fig and tif of max projection
            savefig(maxFig,[imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/' fileName '_maxProject.fig']);
            saveas(maxFig, [imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/' fileName '_maxProject.tif']);
            % max projection with ROI overlay
            [ROIoverlay] = ImBat_ROIoverlay(cellData.results,videoData.Ysiz,centroidFlag,binaryMaskFlag,roiHeatFlag);
            if centroidFlag == 1 || binaryMaskFlag == 1
                %save fig and tif of max projection
                savefig(maxFig,[imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/' fileName '_maxProjectROI.fig']);
                saveas(maxFig, [imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/' fileName '_maxProjectROI.tif']);
            end
            if roiHeatFlag == 1
                savefig(ROIoverlay,[imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/' fileName '_maxProjectROIheatMap.fig']);
                saveas(ROIoverlay, [imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/' fileName '_maxProjectROIheatMap.tif']);
            end
            hold off
        end
        
        if plotFlightsFlag == 1
            %plot all flights in 3D
            if flightPathsAllFlag == 1
                [flightPathsAll,flightPathsStartStop, trajectories_continuous, flight_starts, flight_ends, fstartxyz, fendxyz, batSpeed]...
                    = ImBat_plotFlights(trackData);
                save([imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/' fileName '_flightPaths.mat'],...
                    'trajectories_continuous', 'flight_starts','flight_ends', 'fstartxyz', 'fendxyz', 'batSpeed');
                savefig(flightPathsAll,[imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/' fileName '_flightPathsAll.fig']);
                saveas(flightPathsAll, [imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/' fileName '_flightPathsAll.tif']);
                savefig(flightPathsStartStop,[imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/' fileName '_flightPathsAllStartStop.fig']);
                saveas(flightPathsStartStop, [imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/' fileName '_flightPathsAllStartStop.tif']);
            end
            %plot flights to/from feeder in 3D
            if flightPathsFeederFlag == 1
                [flightPathsToFeeder, flightPathsFromFeeder, flightToFeederEnd, flightToFeederStart, flightFromFeederEnd, flightFromFeederStart]...
                    = ImBat_plotFlightsToFeeder(trackData);
                savefig(flightPathsToFeeder,[imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/' fileName '_flightPathsToFeeder.fig']);
                saveas(flightPathsToFeeder, [imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/' fileName '_flightPathsToFeeder.tif']);
                savefig(flightPathsFromFeeder,[imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/' fileName '_flightPathsFromFeeder.fig']);
                saveas(flightPathsFromFeeder, [imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/' fileName '_flightPathsFromFeeder.tif']);
                save([imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/' fileName '_flightsToFromFeeders.mat'],...
                    'flightToFeederEnd', 'flightToFeederStart', 'flightFromFeederEnd', 'flightFromFeederStart');
            end
            %plot flights over traces
            flightPathsAll = load([imageFolders(kk).folder,'/',imageFolders(kk).name,'/','analysis','/',fileName '_flightPaths.mat']);
            [flightVsVelocity,smoothAvgSpiking,smoothVelocity] = ImBat_plotFlightsVsCells(cellData,alignment,flightPathsAll);
            save([imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/' fileName '_flightPaths.mat'],'smoothVelocity','smoothAvgSpiking','-append');
            savefig(flightVsVelocity,[imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/' fileName '_flightVsVelocity.fig']);
            saveas(flightVsVelocity, [imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/' fileName '_flightVsVelocity.tif']);
        end
        
        
        
    end
    
    
    
end