function [CombinedROI,ROI_Data,flightPaths,c_s_34,FlightAlignedROI] = ImBat_MCS_Align_FlightPaths(aligned_data_struct,cell_registered_struct,master_track_file,ROI_Data,co)

% Load in the neurons!
% Concatonate/Cluster data across days
[CombinedROI,ROI_Data] = ImBat_GroupCalcium(ROI_Data,cell_registered_struct,aligned_data_struct);

% repair bat flight data ( important for z bats with bad tracking) :
ROI_Data = ImBat_RepairFlightData(ROI_Data);

% now, we have a restricted set for ROI_Data, now cluster the flights:
flightPaths = ImBat_GroupFlights(ROI_Data,'mtf',master_track_file,'dist',1.1);         % just the flights

% Sub in the c_s_34 overall clusters for the subset clusters
[ss_1,rr_1] = sort(flightPaths.flight_starts_idx);
c_s_34 = flightPaths.id(rr_1);

% Align to top clusters (co)
for flight_cluster = 1:5%size(flightPaths.clusterIndex,2)
    [FlightAlignedROI{flight_cluster}] = ImBat_Align_FC(CombinedROI,flightPaths,flight_cluster+1);
end 

end
