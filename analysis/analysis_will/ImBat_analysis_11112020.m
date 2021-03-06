function out = ImBat_analysis_11112020(FlightAlignedROI,roi2plot);

% Make heat maps of ROI data on Cells
plotIt = 0;
fl_counter = 0;
hold on;
PlotBoth = 0;
for ii = roi2plot;%1:size(CutCells,1);
hold on;
bound2use = 1:1400;

cell2use = ii;

for iii = 1:size(FlightAlignedROI,2);
    clear ClustFlight CutCells
hold on;
% get flight data, 
ClustFlight = FlightAlignedROI{iii}.ClustFlight_withPads;

% upsample the calcium
%CutCells = FlightAlignedROI{iii}.C_raw;
CutCells = FlightAlignedROI{iii}.S;



for i = 1: size(CutCells,3);
    
    if PlotBoth ==1;
    subplot(1,2,1)
    end
    bound2use = 1:1400;
    colormap(hot);
trial2use = i;
exampFlight = ClustFlight(bound2use,:,trial2use);
exampCell = CutCells(cell2use,:,trial2use);
exampCell = interp(exampCell,4);
exampCell = exampCell(1:size(exampFlight,1));

x = exampFlight(:,1)';
y = exampFlight(:,2)';
z = exampFlight(:,3)';
col = zscore(exampCell)-min(zscore(exampCell));  % This is the color, vary with x in this case.
col_raw = exampCell-min((exampCell));
if sum(abs(diff(col)))>1;
   if plotIt ==1;
p = surface([x;x],[y;y],[z;z],[col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',2);%'AlphaData',mat2gray(col), 'FaceAlpha','flat');
    set(p,'facealpha',0.2)
set(p,'edgealpha',0.2)
   end
fl_counter = fl_counter+1;
end
% TO DO set p == 0 if there is no detected activity... 
if PlotBoth ==1;
    subplot(1,2,2) 
bound2use = 1:1400;
        colormap(hot);
trial2use = i;
exampFlight = ClustFlight(bound2use,:,trial2use);
exampCell = CutCells(cell2use,:,trial2use);
exampCell = interp(exampCell,4);
exampCell = exampCell(1:size(exampFlight,1));

x = 1:length(exampFlight);
y = exampFlight(:,1)';
y2 = exampFlight(:,2)';
y3 = exampFlight(:,3)';
z = zeros(size(x));
col = zscore(exampCell)-min(zscore(exampCell));  % This is the color, vary with x in this case.
%col = exampCell.^(1);  % This is the color, vary with x in this case.
p2 = surface([x;x],[y;y],[z;z],[col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',2);
        set(p2,'facealpha',0.2)        
set(p2,'edgealpha',0.2)
 
 p3 = surface([x;x],[y2;y2],[z;z],[col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',2);
        set(p3,'facealpha',0.2)        
set(p3,'edgealpha',0.2)
 
 p4 = surface([x;x],[y3;y3],[z;z],[col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',2);
        set(p4,'facealpha',0.2)        
set(p4,'edgealpha',0.2)
 title(['ROI ', num2str(cell2use)]);
 end
 % Save data for export:
 out.data{i}.exampFlight = exampFlight;
 out.data{i}.col = col;
 out.data2.AllFlight(:,:,i)  = exampFlight;
 out.data2.AllCol(:,i) = col; % normalized traces
 out.data2.AllCol_raw(:,i)= col_raw; % raw spike counts
end
grid on;
end
% pause();
% clf('reset');
end

title(['ROI: ',num2str(roi2plot), ' , ', num2str(fl_counter), ' Flights total']);
fl_counter; % report how many flights there are
