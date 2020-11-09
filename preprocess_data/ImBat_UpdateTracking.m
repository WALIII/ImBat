function ROI_Data = ImBat_UpdateTracking(ROI_Data,master_track_file)

% Wrapper for all tracking data

OBSOLETE USE 

for i = 1:size(ROI_Data,2);
 ROI_Data{i}.Alignment.out.flights_adjusted  = ImBat_Align_Tracking(ROI_Data{i}.Alignment.out.flights,master_track_file,ROI_Data{i}.date(3:end));
end
