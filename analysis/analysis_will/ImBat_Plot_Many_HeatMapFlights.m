function ImBat_Plot_Many_HeatMapFlights(CombinedROI,flightPaths);

% 11/13/2020
% WAL3


[FlightAlignedROI_1] = ImBat_Align_FC(CombinedROI,flightPaths,2);
[FlightAlignedROI_2] = ImBat_Align_FC(CombinedROI,flightPaths,3);
[FlightAlignedROI_3] = ImBat_Align_FC(CombinedROI,flightPaths,4);
[FlightAlignedROI_4] = ImBat_Align_FC(CombinedROI,flightPaths,5);


clear FlightAlignedROI;
close all

FlightAlignedROI{1} = FlightAlignedROI_1;
FlightAlignedROI{2} = FlightAlignedROI_2;
FlightAlignedROI{3} = FlightAlignedROI_3;
FlightAlignedROI{4} = FlightAlignedROI_4;


% TO DO seperate the xyz plotting to subpanels..
% To Do draw a line though the subpanel as identified burst times...

figure();
for i = 1:21;
ImBat_analysis_11112020(FlightAlignedROI,i);
pause();
clf('reset');
end
