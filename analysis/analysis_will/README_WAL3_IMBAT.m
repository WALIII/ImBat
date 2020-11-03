% Imbat Multi-day Analysis Workflow %

%. colsolidate data across days for a single bat:
[ROI_Data] = ImBat_MultiDayAnalysis(BAT);



% Analysis

% LOAD:
% 1. ROI_Data
% 2. master_track_file
% 3. cell_registered_struct
% 4. aligned_data_struct



% Concatonate/Cluster data across days
flightPaths = ImBat_GroupFlights(ROI_Data,'mtf',master_track_file);         % just the flights
[CombinedROI] = ImBat_GroupCalcium(ROI_Data,cell_registered_struct,aligned_data_struct);


% Bassics of loaded and aligned data;
    % 1. Quality of ROI maps for each day, and across days
    % 2. Quality of flight data for each day, and across days
   
ImBat_PlotBasics(ROI_Data,flightPaths,CombinedROI); % Under construction...

%To DO: ImBat_PlotAdvanced % if ROIs are skipped, are they just not ID'd?



                % SUBFUNCTION: ImBat_Align_FC(CombinedROI,flightPaths,clust2use);
close all


% Figure Making: 
clust2use = 2; % seperate all functions by cluster: TO DO- add to batch...

% Stability over days:
ImBat_analysis_20201029(CombinedROI,flightPaths,clust2use)
 % find first day 

% Figure 1: plot Flights across days
ImBat_analysis_10212020(flightPaths,ROI_Data,CombinedROI,clust2use);

% STATs and tracking quality:
ScoreMatrix = ImBat_analysis_10212026(flightPaths,ROI_Data,CombinedROI,4);
ImBat_Tuning_Stability(ScoreMatrix); % migrate into this funciton 
