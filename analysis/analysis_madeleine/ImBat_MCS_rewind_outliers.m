function [c_s_r,flightPaths_r,num_clust_r] = ImBat_MCS_rewind_outliers(ROI_Data,master_track_file,cell_registered_struct,aligned_data_struct,outliers,rewind_value)

% Rewind all outlier flights to see where they came from
[CombinedROI,ROI_Data] = ImBat_GroupCalcium(ROI_Data,cell_registered_struct,aligned_data_struct);
ROI_Data = ImBat_RepairFlightData(ROI_Data);

flightPaths_r = ImBat_MCS_Groupflights_Rewind(ROI_Data,rewind_value,outliers,'mtf',master_track_file,'dist',1.8);
    
num_clust_r = size(flightPaths_r.clusterIndex,2);
    % Create vector of all flight cluster IDs 
true_times = flightPaths_r.AllFlightsMasterTime(flightPaths_r.flight_starts_idx(:));
[tts,tts_idx]  = sort(true_times(:));
c_s_r = flightPaths_r.id(tts_idx);
end
