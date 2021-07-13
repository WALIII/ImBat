function [cell_registered_struct2] = ImBat_ManualROICorrection(

% mannually intervene and correct ROIs that appear misaligned

ROI2use = 1;
day2use = 1:3;
ImBat_PlotTrackedMasks(ROI_Data,CombinedROI,cell_registered_struct,day2use,ROI2use); % # 3 is cluster 2 is 

