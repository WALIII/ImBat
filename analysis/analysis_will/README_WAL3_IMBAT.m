% Imbat Multi-day Analysis Workflow %

%. colsolidate data across days for a single bat:
[ROI_Data] = ImBat_MultiDayAnalysis(BAT);



% Analysis

% LOAD:
% 1. ROI_Data
% 2. master_track_file
% 3. cell_registered_struct
% 4. aligned_data_struct

% NOTE: this will load data subsets from ROI_Data depending on the CellRegistered mat file. 

% Concatonate/Cluster data across days
[CombinedROI,ROI_Data] = ImBat_GroupCalcium(ROI_Data,cell_registered_struct,aligned_data_struct);

% repair z bat flight data:
ROI_Data = ImBat_RepairFlightData(ROI_Data);

% now, we have a restricted set for ROI_Data, now cluster the flights:
flightPaths = ImBat_GroupFlights(ROI_Data,'mtf',master_track_file,'dist',1.2);         % just the flights
close all

[FlightAlignedROI] = ImBat_Align_FC(CombinedROI,flightPaths,2);
FlightAlignedROI_combined{1} = FlightAlignedROI;

% Bassics of loaded and aligned data;
    % 1. Quality of ROI maps for each day, and across days
    % 2. Quality of flight data for each day, and across days
   
ImBat_PlotBasics(ROI_Data,flightPaths,CombinedROI); % Under construction...

%To DO: ImBat_PlotAdvanced % if ROIs are skipped, are they just not ID'd?



                % SUBFUNCTION: ImBat_Align_FC(CombinedROI,flightPaths,clust2use);
close all


% Figure Making: 
clust2use = 3; % seperate all functions by cluster: TO DO- add to batch...

% Stability over days:
ImBat_analysis_20201029(CombinedROI,flightPaths,2)
 % find first day 

% Figure 1: plot Flights across days
ImBat_analysis_10212020(flightPaths,ROI_Data,CombinedROI,2);

% STATs and tracking quality:
ScoreMatrix = ImBat_analysis_10212026(flightPaths,ROI_Data,CombinedROI,2);
ImBat_Tuning_Stability(ScoreMatrix); % migrate into this funciton 

% plot heat map on flights:
ImBat_analysis_11062020(flightPaths,ROI_Data,CombinedROI,2);


% look at stability for all significant ROIs, on  first 3 flight clusters:
ImBat_analysis_11122020(CombinedROI,flightPaths);

% stats between variability in ROIs and in flight
out2 = ImBat_ROI_Behav_Correlation(FlightAlignedROI_combined)
