function [data_3day] = ImBat_roi3day_loadData(days1to3,varargin)
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
%% load data for each of the days
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
    %load data
    trackDir = dir([flyFolders{i}.folder filesep flyFolders{i}.name filesep analysisNewest filesep '*_flightPaths.mat']);
    batName{i} = trackDir.name(1:3);
    seshTemp = extractAfter(trackDir.name,[batName{i} '_' days1to3{i} '_']);
    sessionType{i} = extractBefore(seshTemp,'_flightPaths.mat');
    results(i) = load([flyFolders{i}.folder filesep flyFolders{i}.name filesep processedNewest filesep 'results.mat']);
    track(i) = load([flyFolders{i}.folder filesep flyFolders{i}.name filesep analysisNewest filesep trackDir.name]);
    video(i) = load([pwd filesep days1to3{i} filesep 'Motion_corrected_Data.mat']);
    alignment(i) = load([flyFolders{i}.folder filesep flyFolders{i}.name filesep processedNewest filesep 'Alignment.mat']);
    close all;
end

data_3day.batName = batName;
data_3day.sessionType = sessionType;
data_3day.results = results;
data_3day.track = track;
data_3day.video = video;
data_3day.alignment = alignment;

