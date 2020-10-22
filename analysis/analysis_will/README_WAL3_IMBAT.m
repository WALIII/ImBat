% Imbat Multi-day Analysis Workflow %

%. colsolidate data across days for a single bat:
[ROI_Data] = ImBat_MultiDayAnalysis(BAT);

% upload master_track_file
%. CellReg

% Concatonate/Cluster data across days
flightPaths = ImBat_GroupFlights(ROI_Data,'mtf',master_track_file);         % just the flights
CombinedROI = ImBat_GroupCalcium(ROI_Data,CellReg2,flightPaths);            % Calcium Imaging Data
