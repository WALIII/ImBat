function [CombinedROI,ROI_Data]= ImBat_GroupCalcium_Single(ROI_Data)

% Save
CombinedROI.C_raw = ROI_Data{1}.ROIs.results.C_raw;
CombinedROI.C = ROI_Data{1}.ROIs.results.C;
CombinedROI.S = full(ROI_Data{1}.ROIs.results.S);
CombinedROI.timestamps =  ROI_Data{1}.Alignment.out.video_times; % align timestamps

% CombinedROI.ROI_all = ROI_all;
% CombinedROI.ROI = ROI;
CombinedROI.day_vector = ones(size(CombinedROI.C_raw,2),1);
