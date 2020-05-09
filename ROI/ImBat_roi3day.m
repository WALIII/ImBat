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

scaling = 4; %scaling for resizing the pixel/frame
%initialize frame vectors for video chunks
frames = cell(length(track),1);
framesCat = cell(length(track),1);
for i = 1:length(track)
    frames{i} = cell(length(track(i).flightPaths.clusterIndex{2}),1); 
    for p = 1:length(track(i).flightPaths.clusterIndex{2})
        %get imaging times of start and stop index converting from tracking to video times
        [minValueStart(p),closestIndexStart(p)] = min(abs(alignment(i).video_timesDS-alignment(i).out.Location_time(track(i).flightPaths.flight_starts_idx(p))));
        [minValueEnd(p),closestIndexEnd(p)] = min(abs(alignment(i).video_timesDS-alignment(i).out.Location_time(track(i).flightPaths.flight_ends_idx(p))));
        %add to frames vector and concatenate
        frames{i}{p} = video(i).Y(:,:,closestIndexStart(p):closestIndexEnd(p));
        framesCat{i} = cat(3,framesCat{i},frames{i}{p});
    end
% Filter movie
framesCat{i} = (convn(framesCat{i}, single(reshape([1 1 1] /10, 1, 1, [])), 'same'));
% Take median of movie
framesCat_med{i} = median(framesCat{i},3);
% Subtract median
framesCat_medNorm{i} = framesCat{i}-framesCat_med{i};
% take max
Ymax{i} = max(framesCat_medNorm{i},[],3);
Ymax{i} = imresize(Ymax{i},scaling);
maxFig = figure();
colormap(gray);
imagesc(Ymax{i});
set(gca,'YDir','normal');
hold on;    
    
end
%change into RGB and overlay colors
[RGB1,RGB2] = CaBMI_XMASS(Ymax{1},Ymax{2},Ymax{3});
figure(); image((RGB1(:,:,:)))


