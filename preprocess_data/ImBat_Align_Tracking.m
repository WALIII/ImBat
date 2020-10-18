function Location_adjusted = ImBat_Align_Tracking(Location,out,date2use);
% Load in Location data, adjust it to the day/day alignment from the
% Master_Tracking_File.mat

% WAL3
% 3/30/2020

% Get appropriate transform from the master tracking file:
Index = find(contains(out.date,date2use));

% Get apply transform from the master tracking file:
TEMP = pointCloud(Location);
ptCloudTformed = pctransform(TEMP,out.tform{Index});
Location_adjusted = ptCloudTformed.Location;

