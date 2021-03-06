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

%topROI is top% of cells you want to look at
topROI = 60;
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
%place cells plot flags
plotPlaceCellsFlag = 1;
%snake/schnitz plot flags
plotSnakesFlag = 1;

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
            mkdir('analysis/ROI');
            % max projection from each session without ROI overlay
            [Ymax, Y, maxFig] = ImBat_Dff(videoData.Y,'batname',batName,'datesesh',dateSesh,'sessiontype',sessionType);
            hold on
            %save fig and tif of max projection
            set(findall(maxFig,'-property','FontSize'),'FontSize',20);
            savefig(maxFig,[imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/ROI/' fileName '_maxProject.fig']);
            saveas(maxFig, [imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/ROI/' fileName '_maxProject.tif']);
            % max projection with ROI overlay
            [ROIoverlay] = ImBat_ROIoverlay(cellData.results,videoData.Ysiz,centroidFlag,binaryMaskFlag,roiHeatFlag,'batname',batName,'datesesh',dateSesh,'sessiontype',sessionType);
            if centroidFlag == 1 || binaryMaskFlag == 1
                %save fig and tif of max projection
                set(findall(maxFig,'-property','FontSize'),'FontSize',20);
                savefig(maxFig,[imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/ROI/' fileName '_maxProjectROI.fig']);
                saveas(maxFig, [imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/ROI/' fileName '_maxProjectROI.tif']);
            end
            if roiHeatFlag == 1
                set(findall(ROIoverlay,'-property','FontSize'),'FontSize',20);
                savefig(ROIoverlay,[imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/ROI/' fileName '_maxProjectROIheatMap.fig']);
                saveas(ROIoverlay, [imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/ROI/' fileName '_maxProjectROIheatMap.tif']);
            end
            hold off
        end
        
        %plot flights in 3d and their clusters
        if plotFlightsFlag == 1 && strcmp(extractBefore(sessionType,'-'),'fly')
            mkdir('analysis/flights')
            %plot all flights in 3D
            if flightPathsAllFlag == 1 && strcmp(extractBefore(sessionType,'-'),'fly')
                [flightPaths] = ImBat_plotFlights(trackData,'batname',batName,'datesesh',dateSesh,'sessiontype',sessionType);
                save([imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/' fileName '_flightPaths.mat'],'flightPaths');
                set(findall(flightPaths.flightPathsAll,'-property','FontSize'),'FontSize',20);
                savefig(flightPaths.flightPathsAll,[imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/flights/' fileName '_flightPathsAll.fig']);
                saveas(flightPaths.flightPathsAll, [imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/flights/' fileName '_flightPathsAll.tif']);
                set(findall(flightPaths.flightPathsStartStop,'-property','FontSize'),'FontSize',20);
                savefig(flightPaths.flightPathsStartStop,[imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/flights/' fileName '_flightPathsAllStartStop.fig']);
                saveas(flightPaths.flightPathsStartStop, [imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/flights/' fileName '_flightPathsAllStartStop.tif']);
                set(findall(flightPaths.flightPathsClusterEach,'-property','FontSize'),'FontSize',20);
                savefig(flightPaths.flightPathsClusterEach,[imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/flights/' fileName '_flightPathsClusterEach.fig']);
                saveas(flightPaths.flightPathsClusterEach, [imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/flights/' fileName '_flightPathsClusterEach.tif']);
                set(findall(flightPaths.flightPathsClusterAll,'-property','FontSize'),'FontSize',20);
                savefig(flightPaths.flightPathsClusterAll,[imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/flights/' fileName '_flightPathsClusterAll.fig']);
                saveas(flightPaths.flightPathsClusterAll, [imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/flights/' fileName '_flightPathsClusterAll.tif']);
            end
            %plot flights to/from feeder in 3D
            if flightPathsFeederFlag == 1 && strcmp(extractBefore(sessionType,'-'),'fly')
                [flightPathsToFeeder, flightPathsFromFeeder,flightPathsClusterToFeederEach, flightPathsClusterToFeederAll, flightFeedersStartStop] = ImBat_plotFlightsToFeeder(trackData,'batname',batName,'datesesh',dateSesh,'sessiontype',sessionType);
                set(findall(flightPathsToFeeder,'-property','FontSize'),'FontSize',20);                
                savefig(flightPathsToFeeder,[imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/flights/' fileName '_flightPathsToFeeder.fig']);
                saveas(flightPathsToFeeder, [imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/flights/' fileName '_flightPathsToFeeder.tif']);
                set(findall(flightPathsFromFeeder,'-property','FontSize'),'FontSize',20); 
                savefig(flightPathsFromFeeder,[imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/flights/' fileName '_flightPathsFromFeeder.fig']);
                saveas(flightPathsFromFeeder, [imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/flights/' fileName '_flightPathsFromFeeder.tif']);
                set(findall(flightPathsClusterToFeederEach,'-property','FontSize'),'FontSize',20); 
                savefig(flightPathsClusterToFeederEach,[imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/flights/' fileName '_flightPathsClusterToFeederEach.fig']);
                saveas(flightPathsClusterToFeederEach, [imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/flights/' fileName '_flightPathsClusterToFeederEach.tif']);
                set(findall(flightPathsClusterToFeederAll,'-property','FontSize'),'FontSize',20); 
                savefig(flightPathsClusterToFeederAll,[imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/flights/' fileName '_flightPathsClusterToFeederAll.fig']);
                saveas(flightPathsClusterToFeederAll, [imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/flights/' fileName '_flightPathsClusterToFeederAll.tif']);
                save([imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/' fileName '_flightsToFromFeeders.mat'],...
                    'flightFeedersStartStop');
            end
            %plot flights over traces
            load([imageFolders(kk).folder,'/',imageFolders(kk).name,'/','analysis','/',fileName '_flightPaths.mat']);
            [flightVsVelocity,smoothAvgSpiking,smoothVelocity] = ImBat_plotFlightsVsCells(cellData,alignment,flightPaths,'batname',batName,'datesesh',dateSesh,'sessiontype',sessionType);
            save([imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/' fileName '_flightPaths.mat'],'smoothVelocity','smoothAvgSpiking','-append');
            set(findall(flightVsVelocity,'-property','FontSize'),'FontSize',20); 
            savefig(flightVsVelocity,[imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/flights/' fileName '_flightVsVelocity.fig']);
            saveas(flightVsVelocity, [imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/flights/' fileName '_flightVsVelocity.tif']);
        end
    
    %number of total flights and rewarded flights
    %numTotalFlights = length(flightPaths.flight_starts_idx);
    %numRewardedFlights = length(flightFeedersStartStop.flightToFeederStart);
    
    %variability of flights
    %covariance matrix

    %place cells
    %spike activity over flight trajectories
    if plotPlaceCellsFlag == 1 && strcmp(extractBefore(sessionType,'-'),'fly')
       mkdir('analysis/placeCells')
       cd([subFolders(i).folder,'/',subFolders(i).name,'/',imageFolders(kk).name,'/analysis/placeCells'])
       ImBat_PlaceCells_Tobias(flightPaths, cellData, alignment,'batname',batName,'datesesh',dateSesh,'sessiontype',sessionType)
       savefig(flightPathsToFeeder,[imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/flights/' fileName '_flightPathsToFeeder.fig']);
       saveas(flightPathsToFeeder, [imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/flights/' fileName '_flightPathsToFeeder.tif']);
 
    end

    %snake plots: raw fluorescent traces for each cell sorted by timing of
    %peak activity for the smoothed zscored mean across all trials in a cluster
    if plotSnakesFlag == 1 && strcmp(extractBefore(sessionType,'-'),'fly')
        cd([subFolders(i).folder,'/',subFolders(i).name,'/',imageFolders(kk).name,'/analysis'])
        mkdir('snakePlots')
        cd([subFolders(i).folder,'/',subFolders(i).name,'/',imageFolders(kk).name,'/analysis/snakePlots'])
        %load([subFolders(i).folder,'/',subFolders(i).name,'/',imageFolders(kk).name,'/analysis/',batName,'_',dateSesh,'_',sessionType,'_flightPaths.mat']);
        [snakeTrace] = ImBat_snakeData(cellData,flightPaths,alignment)        
        [snakeTrace] = ImBat_plotSnake(snakeTrace)
        saveas(snakeTrace.snakePlot_clust, [imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/snakePlots/' fileName '_snakePlots_clust.svg']);
        saveas(snakeTrace.snakePlot_clustOddEven, [imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/snakePlots/' fileName '_snakePlots_clustOddEven.svg']);
        saveas(snakeTrace.snakePlot_clustBy1, [imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/snakePlots/' fileName '_snakePlots_clustBy1.svg']);
        saveas(snakeTrace.snakePlot_prefEachPrePostFlight, [imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/snakePlots/' fileName '_snakePlot_prefEachPrePostFlight.svg']);
        saveas(snakeTrace.snakePlot_clustPrePostFlight, [imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/snakePlots/' fileName '_snakePlot_clustPrePostFlight.svg']);
        saveas(snakeTrace.snakePlot_clust, [imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/snakePlots/' fileName '_snakePlots_clustAll.tif']);
        saveas(snakeTrace.snakePlot_clustOddEven, [imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/snakePlots/' fileName '_snakePlots_clustOddEven.tif']);
        saveas(snakeTrace.snakePlot_clustBy1, [imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/snakePlots/' fileName '_snakePlots_clustBy1.tif']);
        saveas(snakeTrace.snakePlot_prefEachPrePostFlight, [imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/snakePlots/' fileName '_snakePlot_prefEachPrePostFlight.tif']);
        saveas(snakeTrace.snakePlot_clustPrePostFlight, [imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/snakePlots/' fileName '_snakePlot_clustPrePostFlight.tif']);
        saveas(snakeTrace.snakePlot_prefAll, [imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/snakePlots/' fileName '_snakePlot_prefAll.tif']);
        saveas(snakeTrace.snakePlot_prefAll, [imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/snakePlots/' fileName '_snakePlot_prefAll.svg']);
        saveas(snakeTrace.snakePlot_clustBy1PrePostFlight, [imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/snakePlots/' fileName '_snakePlot_clustBy1PrePostFlight.tif']);
        saveas(snakeTrace.snakePlot_clustBy1PrePostFlight, [imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/snakePlots/' fileName '_snakePlot_clustBy1PrePostFlight.svg']);

        save([imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/' fileName '_snakePlotData.mat'],'snakeTrace');
        savefig(snakeTrace.snakePlot_clust, [imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/snakePlots/' fileName '_snakePlots_clustAll.fig']);
        savefig(snakeTrace.snakePlot_clustOddEven, [imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/snakePlots/' fileName '_snakePlots_clustOddEven.fig']);
        savefig(snakeTrace.snakePlot_clustBy1, [imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/snakePlots/' fileName '_snakePlots_clustBy1.fig']);
        savefig(snakeTrace.snakePlot_clustPrePostFlight, [imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/snakePlots/' fileName '_snakePlot_clustPrePostFlight.fig']);
        savefig(snakeTrace.snakePlot_prefEachPrePostFlight, [imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/snakePlots/' fileName '_snakePlot_prefEachPrePostFlight.fig']);
        savefig(snakeTrace.snakePlot_prefAll, [imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/snakePlots/' fileName '_snakePlots_prefAll.fig']);
        savefig(snakeTrace.snakePlot_clustBy1PrePostFlight, [imageFolders(kk).folder '/' imageFolders(kk).name '/analysis/snakePlots/' fileName '_snakePlot_clustBy1PrePostFlight.fig']);

    end
    
    close all;    
    end
    
    
close all;
end