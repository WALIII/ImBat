%% Markov analysis
% format/sort flights
[out_markov] = ImBat_New_Markov(flightPaths34);
% Create prob suff tree
p_min_input = 0.02;
ImBat_ProbSuffixTree(out_markov,5,p_min_input);
% Plotting:
roi_2_use = 1; % which ROI to examine

ImBat_ClusterCalciumVar(FlightAlignedROI{1},roi_2_use); % check if flights cluster based on ROI activity


ImBat_PlotMarkov(out_markov,FlightAlignedROI{1},roi_2_use); 


% now, we have a restricted set for ROI_Data, now cluster the flights:

% check the roi activity is different based on flight order