function ImBat_Analyze

global topROI

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
% d04/19/2020

%topROI is top% of cells you want to look at
topROI = 60;
% Manual inputs
AngeloData = 1;
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
clustManualFlag = 0;
flightPathsFeederFlag = 0;
plotFlightvsCellsFlag = 1;
%place cells plot flags
plotPlaceCellsFlag = 1;
%snake/schnitz plot flags
plotSnakesFlag = 1;
%align max projections of specific flights across trajectories
plotROI3dayFlag = 0;

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
    %     if AngeloData == 1
    %         trackFiles = dir('*track.mat');
    %         cnmfeFiles = dir('CNMFe*');
    %         % load tracking and cell data
    %         if analysisFlag == 1
    %             trackData = load([trackFiles(1).folder,'/',trackFiles(1).name]);
    %             cellData = load([cnmfeFiles(1).folder,'/',cnmfeFiles(1).name]);
    %             videoData = load('Motion_corrected_Data.mat');
    %             alignment = load('Alignment.mat');
    %             alignment.out.video_timesDS = alignment.out.video_times(1:5:end);
    %             fileName = extractBefore(trackFiles(1).name,'_track'); %get filename for titles
    %             dateSesh = trackFiles(1).name(5:10);
    %             batName = trackFiles(1).name(1:3);
    %             sessionType = extractAfter(fileName,[dateSesh '_']);
    %             % make new analysis directory for .fig and .tif files
    %             cd([imageFolders(kk).folder,'/',imageFolders(kk).name]);
    %             mkdir('analysis');
    %             disp('Analyzing!!');
    %         end
    %     else
    extractedFolderFlag = 0;
    % Check if this folder is in the newer extracted version with an extracted folder or older version of cnmfe/extraction
    subFolderFiles = dir(pwd);
    subFolderFiles(ismember( {subFolderFiles.name}, {'.', '..','.DS_Store','Thumbs.db'})) = [];  %remove . and .. and DS_Store
    % Get a logical vector that tells which is a directory.
    subDirFlags = [subFolderFiles.isdir];
    % Extract only those that are directories.
    subDirFolders = subFolderFiles(subDirFlags);
    % Print folder names to command window.
    for k = 1 : length(subDirFolders)
        fprintf('Sub folder #%d = %s\n', k, subDirFolders(k).name);
        if strcmp(subDirFolders(k).name,'extracted') == 1
            cd([subDirFolders(k).folder,'/',subDirFolders(k).name]);
            extractedFolderFlag = 1;
            break
        end
    end
    
    % index every subfolder...
    trackFiles = dir('*track.mat');
    trackFiles(ismember( {trackFiles.name}, {'_track.mat'})) = [];  %'_track.mat'
    imageFolders = dir('*_extraction');
    %print image data folder names
    for kk = 1 : length(imageFolders)
        fprintf('Image data folder #%d = %s\n', kk, imageFolders(kk).name);
        %check that track data matches image data
        %         if strcmp(imageFolders(kk).name(1:end-10),trackFiles(kk).name(1:end-9)) == 0
        %             fprintf('Tracking and image data do not match');
        %             analysisFlag = 0;
        %         end
        
        % Check if folder exists
        if exist([imageFolders(kk).folder,'/',imageFolders(kk).name,'/','analysis*'])>0;
            disp('Folder already analyzed..');
            if reAnalyze ==1
                disp('Re-Analyzing...');
            else
                disp('Moving to the next folder...');
                analysisFlag = 0 ;
            end
        end
        
        processedFolders = dir([imageFolders(kk).folder,'/',imageFolders(kk).name,'/','processed*']);
        processedNewest = sort({processedFolders(end).name});
        processedNewest = char(processedNewest);
        
        
        
        % load tracking and cell data
        if analysisFlag == 1
            trackData = load([trackFiles(kk).folder,'/',trackFiles(kk).name]);
            cellData = load([imageFolders(kk).folder,'/',imageFolders(kk).name,'/',processedNewest,'/','results.mat']);%'Motion_corrected_Data_DS_results.mat']);
            videoData = load([imageFolders(kk).folder,'/',imageFolders(kk).name,'/',processedNewest,'/','Motion_corrected_Data_DS.mat']);
            alignment = load([imageFolders(kk).folder,'/',imageFolders(kk).name,'/',processedNewest,'/','Alignment.mat']);
            video_timesDS = alignment.out.video_times(1:5:end);
            alignment.out.video_timesDS = video_timesDS;
            fileName = extractBefore(imageFolders(kk).name,'_extraction'); %get filename for titles
            if extractedFolderFlag == 1
                dateSesh = imageFolders(kk).folder(end-15:end-10);
            else
                dateSesh = imageFolders(kk).folder(end-5:end);
            end
            batName = extractBefore(fileName,['_' dateSesh]);
            sessionType = extractAfter(fileName,[dateSesh '_']);
            save([imageFolders(kk).folder '/' imageFolders(kk).name '/' processedNewest '/Alignment.mat'],'video_timesDS','-append');
            % make new analysis directory for .fig and .tif files
            cd([imageFolders(kk).folder,'/',imageFolders(kk).name]);
            analysis_Folder = ['analysis_',datestr(now,'yyyy_mm_dd__hhMM')];
            mkdir(analysis_Folder);
            disp('Analyzing!!');
        end
        %end
        %         if strcmp(extractBefore(sessionType,'-'),'fly')
        %             plotFlightsFlag = 1;
        %             flightPathsAllFlag = 1;
        %             flightPathsFeederFlag = 1;
        %             plotPlaceCellsFlag = 1;
        %         else
        %             plotFlightsFlag = 0;
        %             flightPathsAllFlag = 0;
        %             flightPathsFeederFlag = 0;
        %             plotPlaceCellsFlag = 0;
        %         end
        %%perform basic sanity checks on imaging and flight data
        if plotROIFlag == 1
            mkdir([analysis_Folder filesep 'ROI']);
            % max projection from each session without ROI overlay
            [Ymax, Y, maxFig] = ImBat_Dff(videoData.Y,'batname',batName,'datesesh',dateSesh,'sessiontype',sessionType);
            hold on
            %save fig and tif of max projection
            %set(findall(maxFig,'-property','FontSize'),'FontSize',20);
            savefig(maxFig,[imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/ROI/' fileName '_maxProject.fig']);
            saveas(maxFig, [imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/ROI/' fileName '_maxProject.tif']);
            % max projection with ROI overlay
            [ROIoverlay,centroidMax] = ImBat_ROIoverlay(cellData.results,'centroid',centroidFlag,'binarymask',binaryMaskFlag,'roiheat',roiHeatFlag,'batname',batName,'datesesh',dateSesh,'sessiontype',sessionType);
            if centroidFlag == 1 || binaryMaskFlag == 1
                %save fig and tif of max projection
                %set(findall(centroidMax,'-property','FontSize'),'FontSize',20);
                savefig(centroidMax,[imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/ROI/' fileName '_maxProjectROI.fig']);
                saveas(centroidMax, [imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/ROI/' fileName '_maxProjectROI.tif']);
            end
            if roiHeatFlag == 1
                %set(findall(ROIoverlay,'-property','FontSize'),'FontSize',20);
                savefig(ROIoverlay,[imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/ROI/' fileName '_maxProjectROIheatMap.fig']);
                saveas(ROIoverlay, [imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/ROI/' fileName '_maxProjectROIheatMap.tif']);
            end
            hold off
        end
        
        %plot flights in 3d and their clusters
        if plotFlightsFlag == 1 && strcmp(extractBefore(sessionType,'-'),'fly')
            mkdir([analysis_Folder filesep 'flights'])
            %plot all flights in 3D
            if flightPathsAllFlag == 1 && strcmp(extractBefore(sessionType,'-'),'fly')
                %[flightPaths] = ImBat_plotFlights(trackData,'batname',batName,'datesesh',dateSesh,'sessiontype',sessionType,'clustmanualflag',clustManualFlag);
                [flightPaths] = ImBat_flightsAng(trackData,'batname',batName,'datesesh',dateSesh,'sessiontype',sessionType);
                save([imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/' fileName '_flightPaths.mat'],'flightPaths');
                %set(findall(flightPaths.flightPathsAll,'-property','FontSize'),'FontSize',20);
                savefig(flightPaths.flightPathsAll,[imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/flights/' fileName '_flightPathsAll.fig']);
                saveas(flightPaths.flightPathsAll, [imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/flights/' fileName '_flightPathsAll.tif']);
                saveas(flightPaths.flightPathsAll, [imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/flights/' fileName '_flightPathsAll.svg']);
                %set(findall(flightPaths.flightPathsStartStop,'-property','FontSize'),'FontSize',20);
                savefig(flightPaths.flightPathsStartStop,[imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/flights/' fileName '_flightPathsAllStartStop.fig']);
                saveas(flightPaths.flightPathsStartStop, [imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/flights/' fileName '_flightPathsAllStartStop.tif']);
                saveas(flightPaths.flightPathsStartStop, [imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/flights/' fileName '_flightPathsAllStartStop.svg']);
                %set(findall(flightPaths.flightPathsClusterEach,'-property','FontSize'),'FontSize',20);
                savefig(flightPaths.flightPathsClusterEach,[imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/flights/' fileName '_flightPathsClusterEach.fig']);
                saveas(flightPaths.flightPathsClusterEach, [imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/flights/' fileName '_flightPathsClusterEach.tif']);
                saveas(flightPaths.flightPathsClusterEach, [imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/flights/' fileName '_flightPathsClusterEach.svg']);
                savefig(flightPaths.flightTimeline,[imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/flights/' fileName '_flightTimeline.fig']);
                saveas(flightPaths.flightTimeline, [imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/flights/' fileName '_flightTimeline.tif']);
                saveas(flightPaths.flightTimeline, [imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/flights/' fileName '_flightTimeline.svg']);
                savefig(flightPaths.clusterDistance,[imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/flights/' fileName '_clusterDistance.fig']);
                saveas(flightPaths.clusterDistance, [imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/flights/' fileName '_clusterDistance.tif']);
                saveas(flightPaths.clusterDistance, [imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/flights/' fileName '_clusterDistance.svg']);
                %set(findall(flightPaths.flightPathsClusterAll,'-property','FontSize'),'FontSize',20);
                %savefig(flightPaths.flightPathsClusterAll,[imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/flights/' fileName '_flightPathsClusterAll.fig']);
                %saveas(flightPaths.flightPathsClusterAll, [imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/flights/' fileName '_flightPathsClusterAll.tif']);
            end
            %plot flights to/from feeder in 3D
            if flightPathsFeederFlag == 1 && strcmp(extractBefore(sessionType,'-'),'fly')
                [flightPathsToFeeder, flightPathsFromFeeder,flightPathsClusterToFeederEach, flightPathsClusterToFeederAll, flightFeedersStartStop] = ImBat_plotFlightsToFeeder(trackData,'batname',batName,'datesesh',dateSesh,'sessiontype',sessionType);
                %set(findall(flightPathsToFeeder,'-property','FontSize'),'FontSize',20);
                savefig(flightPathsToFeeder,[imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/flights/' fileName '_flightPathsToFeeder.fig']);
                saveas(flightPathsToFeeder, [imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/flights/' fileName '_flightPathsToFeeder.tif']);
                %set(findall(flightPathsFromFeeder,'-property','FontSize'),'FontSize',20);
                savefig(flightPathsFromFeeder,[imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/flights/' fileName '_flightPathsFromFeeder.fig']);
                saveas(flightPathsFromFeeder, [imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/flights/' fileName '_flightPathsFromFeeder.tif']);
                %set(findall(flightPathsClusterToFeederEach,'-property','FontSize'),'FontSize',20);
                savefig(flightPathsClusterToFeederEach,[imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/flights/' fileName '_flightPathsClusterToFeederEach.fig']);
                saveas(flightPathsClusterToFeederEach, [imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/flights/' fileName '_flightPathsClusterToFeederEach.tif']);
                %set(findall(flightPathsClusterToFeederAll,'-property','FontSize'),'FontSize',20);
                savefig(flightPathsClusterToFeederAll,[imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/flights/' fileName '_flightPathsClusterToFeederAll.fig']);
                saveas(flightPathsClusterToFeederAll, [imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/flights/' fileName '_flightPathsClusterToFeederAll.tif']);
                save([imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/' fileName '_flightsToFromFeeders.mat'],...
                    'flightFeedersStartStop');
            end
            % plot flights over traces
            if plotFlightvsCellsFlag == 1
                load([imageFolders(kk).folder,'/',imageFolders(kk).name,'/',analysis_Folder,'/',fileName '_flightPaths.mat']);
                [flightVsVelocity,smoothAvgSpiking,smoothVelocity] = ImBat_plotFlightsVsCells(cellData,alignment,flightPaths,'batname',batName,'datesesh',dateSesh,'sessiontype',sessionType);
                save([imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/' fileName '_flightPaths.mat'],'smoothVelocity','smoothAvgSpiking','-append');
                %set(findall(flightVsVelocity,'-property','FontSize'),'FontSize',20);
                savefig(flightVsVelocity,[imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/flights/' fileName '_flightVsVelocity.fig']);
                saveas(flightVsVelocity, [imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/flights/' fileName '_flightVsVelocity.tif']);
            end
        end
        
        %number of total flights and rewarded flights
        %numTotalFlights = length(flightPaths.flight_starts_idx);
        %numRewardedFlights = length(flightFeedersStartStop.flightToFeederStart);
        
        %variability of flights
        %covariance matrix
        
        %place cells
        %spike activity over flight trajectories
        if plotPlaceCellsFlag == 1 && strcmp(extractBefore(sessionType,'-'),'fly')
            mkdir([imageFolders(kk).folder,'/',imageFolders(kk).name,'/' analysis_Folder '/placeCells']);
            cd([imageFolders(kk).folder,'/',imageFolders(kk).name,'/' analysis_Folder '/placeCells']);
            ImBat_PlaceCells_Tobias(flightPaths, cellData, alignment,'batname',batName,'datesesh',dateSesh,'sessiontype',sessionType)
            %savefig(flightPathsToFeeder,[imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/flights/' fileName '_flightPathsToFeeder.fig']);
            %saveas(flightPathsToFeeder, [imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/flights/' fileName '_flightPathsToFeeder.tif']);
            
        end
        
        %snake plots: raw fluorescent traces for each cell sorted by timing of
        %peak activity for the smoothed zscored mean across all trials in a cluster
        if plotSnakesFlag == 1 && strcmp(extractBefore(sessionType,'-'),'fly')
            cd([imageFolders(kk).folder ,'/',imageFolders(kk).name,'/' analysis_Folder])
            mkdir('snakePlots')
            cd([imageFolders(kk).folder,'/',imageFolders(kk).name,'/' analysis_Folder '/snakePlots'])
            %load([subFolders(i).folder,'/',subFolders(i).name,'/',imageFolders(kk).name,'/' analysis_Folder '/',batName,'_',dateSesh,'_',sessionType,'_flightPaths.mat']);
            [snakeTrace_cRaw,snakeTrace_c,snakeTrace_s] = ImBat_snakeData(cellData,flightPaths,alignment)
            [snakeTrace] = ImBat_plotSnake(snakeTrace_cRaw)
            saveas(snakeTrace.snakePlot_clust, [imageFolders(kk).folder filesep imageFolders(kk).name filesep analysis_Folder filesep 'snakePlots' filesep fileName '_snakePlots_clust.svg']);
            saveas(snakeTrace.snakePlot_clustOddEven, [imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/snakePlots/' fileName '_snakePlots_clustOddEven.svg']);
            saveas(snakeTrace.snakePlot_clustBy1, [imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/snakePlots/' fileName '_snakePlots_clustBy1.svg']);
            saveas(snakeTrace.snakePlot_prefEachPrePostFlight, [imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/snakePlots/' fileName '_snakePlot_prefEachPrePostFlight.svg']);
            saveas(snakeTrace.snakePlot_clustPrePostFlight, [imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/snakePlots/' fileName '_snakePlot_clustPrePostFlight.svg']);
            saveas(snakeTrace.snakePlot_clust, [imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/snakePlots/' fileName '_snakePlots_clustAll.tif']);
            saveas(snakeTrace.snakePlot_clustOddEven, [imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/snakePlots/' fileName '_snakePlots_clustOddEven.tif']);
            saveas(snakeTrace.snakePlot_clustBy1, [imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/snakePlots/' fileName '_snakePlots_clustBy1.tif']);
            saveas(snakeTrace.snakePlot_prefEachPrePostFlight, [imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/snakePlots/' fileName '_snakePlot_prefEachPrePostFlight.tif']);
            saveas(snakeTrace.snakePlot_clustPrePostFlight, [imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/snakePlots/' fileName '_snakePlot_clustPrePostFlight.tif']);
            saveas(snakeTrace.snakePlot_prefAll, [imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/snakePlots/' fileName '_snakePlot_prefAll.tif']);
            saveas(snakeTrace.snakePlot_prefAll, [imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/snakePlots/' fileName '_snakePlot_prefAll.svg']);
            saveas(snakeTrace.snakePlot_clustBy1PrePostFlight, [imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/snakePlots/' fileName '_snakePlot_clustBy1PrePostFlight.tif']);
            saveas(snakeTrace.snakePlot_clustBy1PrePostFlight, [imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/snakePlots/' fileName '_snakePlot_clustBy1PrePostFlight.svg']);
            
            save([imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/' fileName '_snakePlotData.mat'],'snakeTrace');
            savefig(snakeTrace.snakePlot_clust, [imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/snakePlots/' fileName '_snakePlots_clustAll.fig']);
            savefig(snakeTrace.snakePlot_clustOddEven, [imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/snakePlots/' fileName '_snakePlots_clustOddEven.fig']);
            savefig(snakeTrace.snakePlot_clustBy1, [imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/snakePlots/' fileName '_snakePlots_clustBy1.fig']);
            savefig(snakeTrace.snakePlot_clustPrePostFlight, [imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/snakePlots/' fileName '_snakePlot_clustPrePostFlight.fig']);
            savefig(snakeTrace.snakePlot_prefEachPrePostFlight, [imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/snakePlots/' fileName '_snakePlot_prefEachPrePostFlight.fig']);
            savefig(snakeTrace.snakePlot_prefAll, [imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/snakePlots/' fileName '_snakePlots_prefAll.fig']);
            savefig(snakeTrace.snakePlot_clustBy1PrePostFlight, [imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/snakePlots/' fileName '_snakePlot_clustBy1PrePostFlight.fig']);
            
        end
        
        %snake plots: raw fluorescent traces for each cell sorted by timing of
        %peak activity for the smoothed zscored mean across all trials in a cluster
        if plotROI3dayFlag == 1 && strcmp(extractBefore(sessionType,'-'),'fly')
            cd([imageFolders(kk).folder ,'/',imageFolders(kk).name,'/' analysis_Folder])
            mkdir('roi3day')
            cd([imageFolders(kk).folder,'/',imageFolders(kk).name,'/' analysis_Folder '/roi3day'])
            %load([subFolders(i).folder,'/',subFolders(i).name,'/',imageFolders(kk).name,'/' analysis_Folder '/',batName,'_',dateSesh,'_',sessionType,'_flightPaths.mat']);
            [roi3day] = ImBat_roi3day(videoData,flightPaths,alignment)
            
        end
        
        close all;
    end
    
    
    close all;
end