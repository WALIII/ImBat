function master_track_file = ImBat_Create_Master_Track

% run in folder with calibration .c3d files 
% WAL3
% d03/31/2020



% Housekeeping:
out.Creation_date = datestr(now);

% STEP 1: get all .c3d files in folder:
DIR = pwd;
mov_listing=dir(fullfile(DIR,'*.c3d'));
mov_listing={mov_listing(:).name};

% select Reference Day
RefDay  = 23;


for i=1:length(mov_listing)

% STEP 2: extract Markers
    [path,file,ext]=fileparts(mov_listing{i});
    FILE = fullfile(DIR,mov_listing{i})

[Markers,~,~,~,~,~,~,~]=readC3D_analog(FILE);
out.date{i} = file(7:12); % read date from filename
out.Markers_Raw{i} = Markers; % save markers
end


% STEP 3: Align Markers to first day
disp('aligning all days...');
out.ReferenceDate = out.date{RefDay};
for i = 1:size(out.date,2)
[out.tform{i},out.Markers_Adjusted{i}] = ImBat_AlignPoints(out.Markers_Raw{RefDay},out.Markers_Raw{i});

% % STEP 4: Use tform to adjust markers ( for reference)
% 
% TEMP = pointCloud(out.Markers_Raw{i});
% ptCloudTformed = pctransform(TEMP,out.tform{i});
% out.Markers_Adjusted = ptCloudTformed.Location;

end

master_track_file = out;

save(['Master_Tracking_File_',datestr(now),'.mat'],'master_track_file');

