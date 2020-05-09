function [roi3day] = ImBat_roi3day(days1to3,varargin)

batName = [];
dateSesh = [];
sessionType = [];
loadFlag = 0;

% User inputs overrides
nparams=length(varargin);
for i=1:2:nparams
    switch lower(varargin{i})
        case 'batname'
            batName=varargin{i+1};
        case 'datesesh'
            dateSesh = varargin{i+1};
        case 'sessiontype'
            sessionType = varargin{i+1};
        case 'loadflag'
            loadFlag = varargin{i+1};
        case 'saveflag'
            saveFlag = varargin{i+1};
    end
end

datedir = dir('20*');
flyfolders = [];
for i = 1:length(days1to3)
   datadir = dir([days1to3{i} '\extracted\'  '*fly*']);
   dirFlags = [datadir.isdir];
    % Extract only those that are directories.
    flyFolders{i} = datadir(dirFlags);
   %enter newest processed folder
    processedFolders = dir([flyFolders{i}.folder,'/',flyFolders{i}.name,'/','processed*']);
    processedNewest = sort({processedFolders(end).name});
    processedNewest = char(processedNewest);
    analysisFolders = dir([flyFolders{i}.folder,'/',flyFolders{i}.name,'/','analysis*']);
    analysisNewest = sort({analysisFolders(end).name});
    analysisNewest = char(analysisNewest);
    trackDir = dir([flyFolders{i}.folder filesep flyFolders{i}.name filesep analysisNewest filesep '*_flightPaths.mat']);
    batName{i} = trackDir.name(1:3);
    seshTemp = extractAfter(trackDir.name,[batName{i} '_' days1to3{i} '_']);
    sessionType{i} = extractBefore(seshTemp,'_flightPaths.mat');
    results(i) = load([flyFolders{i}.folder filesep flyFolders{i}.name filesep processedNewest filesep 'results.mat']);
    track(i) = load([flyFolders{i}.folder filesep flyFolders{i}.name filesep analysisNewest filesep trackDir.name]);
    video(i) = load([pwd filesep days1to3{i} filesep 'Motion_corrected_Data.mat']);
    alignment(i) = load([pwd filesep days1to3{i} filesep 'Alignment.mat']);
    close all;
end






% ImBat_Dff
scaling = 4;

for i = 1:length(track)
    for p = 1:length(track(i).flightPaths.clusterIndex{2})
        %get imaging times of start and stop index converting from tracking to video times
        [minValueStart(p),closestIndexStart(p)] = min(abs(alignment(i).out.video_times-alignment.out.Location_time(flightPaths.flight_starts_idx(dur_i))));
        [minValueEnd(dur_i),closestIndexEnd(dur_i)] = min(abs(alignment.out.video_times-alignment.out.Location_time(flightPaths.flight_ends_idx(dur_i))));
        
        frames{i}(p) = video(i).Y(:,:,track(i).flightPaths.flight
    end
end


