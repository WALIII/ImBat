function ImBat_Plot_HeatCells(flightPaths,CombinedROI);

% FlightData = flightPaths.trajectoriesContinous;

ROIDat = CombinedROI.C_raw;
timevec = CombinedROI.timestamps;
exampFlight = flightPaths.trajectoriesContinous';
% % get ROI data conditioned...
% cell2use = 1;
% exampCell = ROIDat(cell2use,:);
% ave_fs = 30;
% ave_time=14:1/ave_fs:max(flightPaths.AllFlightsMasterTime);
% x = timevec';
% v = exampCell;
% xq = ave_time;
% 
% [~, ind] = unique(x); % ind = index of first occurrence of a repeated value 
% vq = interp1(x(ind)',v(ind)',xq','spline')';
% 
% figure();
% hold on;
% plot(vq,'r');
% plot(exampCell,'b');
% clear exampCell
% 
% exampCell = interp(vq,4);





  
  figure();
  
for ii = 1:60;


cell2use = ii;
  
    colormap(hot);


 % get ROI data conditioned...
cell2use = ii;
exampCell = ROIDat(cell2use,:);
ave_fs = 30;
ave_time=min(flightPaths.AllFlightsMasterTime):1/ave_fs:max(flightPaths.AllFlightsMasterTime);
x = timevec';
v = exampCell;
xq = ave_time;

[~, ind] = unique(x); % ind = index of first occurrence of a repeated value 
vq = interp1(x(ind)',v(ind)',xq','spline')';
exampCell = interp(vq,4);   
    
    
exampCell = exampCell(1:size(exampFlight,1));
exampCell(exampCell<-1) =0;
x = exampFlight(:,1)';
y = exampFlight(:,2)';
z = exampFlight(:,3)';
col = exampCell;  % This is the color, vary with x in this case.
% col = col-mode(col);
surface([x;x],[y;y],[z;z],[col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',2);
caxis([0 5])
pause();
clear exampCell
clf('reset');
end
