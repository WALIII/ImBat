function [placeCells] = ImBat_extract_noResponseCells(placeCells)
nClusters = size(placeCells.pp_cells,1);
roiDir = dir('*ROI_*');
roiNum = zeros(1,length(roiDir));
for roi_i = 1:length(roiDir)
   %find all ROI numbers
   roiLong = extractAfter(roiDir(roi_i).name,'ROI_');
   roiNum(roi_i) = str2num(extractBefore(roiLong,'_cluster'));
    
end

roiNone = [];
for idx_i = 1:length(roiNum)
    try
consecIdx = roiNum(idx_i:idx_i+nClusters-2);
    catch
    end
    
    if sum(consecIdx)/3 == consecIdx(1)
        roiNone = [roiNone consecIdx(1)]; 
    end
end

placeCells.cells_none = roiNone;

