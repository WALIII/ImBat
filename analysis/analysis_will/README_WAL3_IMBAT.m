% Imbat Multi-day Analysis Workflow %

%. consolidate data across days for a single bat ( run in directory with your data):
[ROI_Data] = ImBat_MultiDayAnalysis(BAT); % BAT is a string, like 'za'

% Check extractions
[Extraction_Metadata] = ImBat_PreCheck(ROI_Data);

% Analysis

% LOAD:
% 1. ROI_Data
% 2. master_track_file
% 3. cell_registered_struct
% 4. aligned_data_struct

% NOTE: this will load data subsets from ROI_Data depending on the CellRegistered mat file. 


%% PRE-FLIGHT:
% Concatonate/Cluster data across days
[CombinedROI,ROI_Data] = ImBat_GroupCalcium(ROI_Data,cell_registered_struct,aligned_data_struct);

% repair z bat flight data:
ROI_Data = ImBat_RepairFlightData(ROI_Data);

% now, we have a restricted set for ROI_Data, now cluster the flights:
flightPaths = ImBat_GroupFlights(ROI_Data,'mtf',master_track_file,'dist',1.2);         % just the flights

close all


%% 
% find first day % Align Flight data to the top 3 flight clusters:
  for flight_cluster = 1:3;
    [FlightAlignedROI{flight_cluster}] = ImBat_Align_FC(CombinedROI,flightPaths,flight_cluster+1);
  end 

%% Save data
 mkdir('Saved_Data')
 save('Saved_Data/Aligned_Data.mat','flightPaths','CombinedROI','FlightAlignedROI','-v7.3');
 save('ROI_Data.mat', 'ROI_Data'); % Save fixed flights


% Basics of loaded and aligned data;
    % 1. Quality of ROI maps for each day, and across days
    % 2. Quality of flight data for each day, and across days
   
%ImBat_PlotBasics(ROI_Data,flightPaths,CombinedROI); % Under construction...

%To DO: ImBat_PlotAdvanced % if ROIs are skipped, are they just not ID'd?

% SUBFUNCTION: ImBat_Align_FC(CombinedROI,flightPaths,clust2use);


% Figure Making: 

% Stability over days:
% ImBat_analysis_20201029(CombinedROI,flightPaths,2)
 


% Get useful flights from each day
output = ImBat_Quantify_Flights(flightPaths,ROI_Data);

% Now, use the integrated function to find good ( or bad) cells. 
% FlightAlignedROI{1} is the first clustered flight:
ImBat_PlotAlignedROIs(FlightAlignedROI{1});



% Figure 1: plot Flights across days, and PLOT ROIs
ImBat_analysis_10212020(flightPaths,ROI_Data,CombinedROI,2);

% STATs and tracking quality:
ScoreMatrix = ImBat_analysis_10212026(flightPaths,ROI_Data,CombinedROI,2);
figure();
ImBat_Tuning_Stability(ScoreMatrix,'r'); % migrate into this funciton

% plot heat map on flights:
ImBat_analysis_11062020(flightPaths,ROI_Data,CombinedROI,2);


% look at stability for all significant ROIs, on  first 3 flight clusters:
ImBat_analysis_11122020(CombinedROI,flightPaths);

% stats between variability in ROIs and in flight
out2 = ImBat_ROI_Behav_Correlation(FlightAlignedROI_combined)




% % ROIs

% Manually selected ROIs to investigate 
cells2use = [ 24 26 39 48 49 51 61 70 82 87]; % Unstable ROIs to use

% Create Labeled dff projections:
ImBat_PlotTrackedMasks(ROI_Data,CombinedROI,cell_registered_struct,output.days2use{3}',cells2use); % # 3 is cluster 2 is 


%% Markov analysis
[out_markov] = ImBat_New_Markov(flightPaths);
ImBat_ProbSuffixTree(out_markov,5);

roi_2_use = 1;

ImBat_ClusterCalciumVar(FlightAlignedROI{1},roi_2_use); % check the flights
ImBat_PlotMarkov(out_markov,FlightAlignedROI{1},roi_2_use); % check the roi
