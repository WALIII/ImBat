function [c_s_34,flightPaths34,num_clust,fd] = ImBat_MCS_Load_Flight_Sequence_Data(ROI_Data,master_track_file,cell_registered_struct,aligned_data_struct)


    % Does first few steps of Basic Processing Walkthrough (https://github.com/WALIII/ImBat/wiki/Basic-Processing-Walkthrough)
    % after the 4 files have been manually loaded

    [CombinedROI,ROI_Data] = ImBat_GroupCalcium(ROI_Data,cell_registered_struct,aligned_data_struct);
    ROI_Data = ImBat_RepairFlightData(ROI_Data);
    ROI_Data_34 = ROI_Data;

    % Construct the flightPaths struct (flightPaths34). This script calls
    % ImBat_flightAngelo_MCS which does the flight segregation.
    n_surv_clusters = 35;
    flightPaths34 = ImBat_GroupFlights(ROI_Data_34,'dist',1.8);
    
    num_clust = size(flightPaths34.clusterIndex,2);
    % Create vector of all flight cluster IDs 
    true_times = flightPaths34.AllFlightsMasterTime(flightPaths34.flight_starts_idx(:));
    [tts,tts_idx]  = sort(true_times(:));
    c_s_34 = flightPaths34.id(tts_idx);
    % Create vector of all flight days
    fd = flightPaths34.day(tts_idx);

end