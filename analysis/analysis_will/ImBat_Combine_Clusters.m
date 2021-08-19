function [FlightAlignedROI_new2] = ImBat_Combine_Clusters(CombinedROI,flightPaths);


FL2use = [2:6]; % use the two largest flight clusters:




% extract the top 6 paths
for i = 1:6
[FlightAlignedROI2{i}] = ImBat_Align_FC(CombinedROI,flightPaths,i);
end

close all;
% concatonate flight aligned data




clear FlightAlignedROI_new2
for i = 1:size(FL2use,2);
    if i ==1; 
FlightAlignedROI_new2{1}.S = FlightAlignedROI2{FL2use(i)}.S;
FlightAlignedROI_new2{1}.ClustFlight_withPads  =  FlightAlignedROI2{FL2use(i)}.ClustFlight_withPads;  
FlightAlignedROI_new2{1}.cluster_idX  =  FlightAlignedROI2{FL2use(i)}.cluster_idX;  
FlightAlignedROI_new2{1}.cluster_ID  =  ones(1,size(FlightAlignedROI2{FL2use(i)}.ClustFlight_withPads,3))*FlightAlignedROI2{FL2use(i)}.clust_number;
FlightAlignedROI_new2{1}.CutCells_date =  FlightAlignedROI2{FL2use(i)}.CutCells_date;
    else 
FlightAlignedROI_new2{1}.S = cat(3,FlightAlignedROI_new2{1}.S, FlightAlignedROI2{FL2use(i)}.S);
FlightAlignedROI_new2{1}.ClustFlight_withPads  = cat(3,FlightAlignedROI_new2{1}.ClustFlight_withPads,FlightAlignedROI2{FL2use(i)}.ClustFlight_withPads); 
FlightAlignedROI_new2{1}.cluster_idX  = cat(2,FlightAlignedROI_new2{1}.cluster_idX,FlightAlignedROI2{FL2use(i)}.cluster_idX); 
FlightAlignedROI_new2{1}.cluster_ID  =  cat(2,FlightAlignedROI_new2{1}.cluster_ID,ones(1,size(FlightAlignedROI2{FL2use(i)}.ClustFlight_withPads,3))*FlightAlignedROI2{FL2use(i)}.clust_number);
FlightAlignedROI_new2{1}.CutCells_date  =  cat(2,FlightAlignedROI_new2{1}.CutCells_date,FlightAlignedROI2{FL2use(i)}.CutCells_date);

    end
end