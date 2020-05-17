function [framesData] = ImBat_roi3day(days1to3,varargin)
%%
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
%%
figure();
sgtitle([batName{1} ': cluster 2']);
scaling = 4; %scaling for resizing the pixel/frame
%initialize frame vectors for video chunks
frames = cell(length(track),1);
framesCat = cell(length(track),1);
framesCat4 = cell(length(track),1);
framesLength =[];
minValueStart = [];
minValueEnd = [];
closestIndexStart = [];
closestIndexEnd = [];
randInd = cell(length(track),1);
flightNums = cell(length(track),1);
%find min number of flights in each day to subsample later
for i = 1:length(track)
    flightNums{i} = track(i).flightPaths.clusterIndex{2}'; %get the indices for the flights in the cluster
    flightNumsLength(i) = length(flightNums{i}); %number of flights on each day
end
minFlightNums = min(flightNumsLength); %least number of flights across the 3 days
for i = 1:length(track)
    randInd{i} = randperm(flightNumsLength(i),minFlightNums); %generate random list with minlength of flights for random selection %flightNumsLength(i));
    for p = 1:minFlightNums %for each flight within cluster #2 %length(randInd)
        %get imaging times of start and stop index converting from tracking to video times
        [minValueStart{i}(p),closestIndexStart{i}(p)] = min(abs(alignment(i).video_timesDS-alignment(i).out.Location_time(track(i).flightPaths.flight_starts_idx(track(i).flightPaths.clusterIndex{2}(randInd{i}(p))))));
        [minValueEnd{i}(p),closestIndexEnd{i}(p)] = min(abs(alignment(i).video_timesDS-alignment(i).out.Location_time(track(i).flightPaths.flight_ends_idx(track(i).flightPaths.clusterIndex{2}(randInd{i}(p))))));
        framesLength{i}(p) = closestIndexEnd{i}(p) - closestIndexStart{i}(p);
    end
    framesMax(i) = max(framesLength{i}); %find max # frames across all flights
end
    framesMaxAll = max(framesMax);
for i = 1:length(track)
    frames{i} = cell(length(track(i).flightPaths.clusterIndex{2}),1); 
    for p = 1:minFlightNums %length(randInd)
        framesDiff{i}(p) = framesMaxAll - framesLength{i}(p);
        %add to frames vector and concatenate
        frames{i}{p} = video(i).Y(:,:,closestIndexStart{i}(p):(closestIndexEnd{i}(p)+framesDiff{i}(p)));
        framesCat{i} = cat(3,framesCat{i},frames{i}{p});
        framesCat4{i} = cat(4,framesCat4{i},frames{i}{p});
        subplot(4,3,i);
        plot3(track(i).flightPaths.pos(1,:,track(i).flightPaths.clusterIndex{2}(randInd{i}(p))),track(i).flightPaths.pos(2,:,track(i).flightPaths.clusterIndex{2}(randInd{i}(p))),track(i).flightPaths.pos(3,:,track(i).flightPaths.clusterIndex{2}(randInd{i}(p))));
        hold on;
    end
    framesCat4_mean{i} = mean(framesCat4{i},4);
    %make gaussian filter based on normcorre function
    gSig = 2; 
    gSiz = 5*gSig; 
    psf = fspecial('gaussian', round(2*gSiz), gSig);
    ind_nonzero = (psf(:)>=max(psf(:,1)));
    psf = psf-mean(psf(ind_nonzero));
    psf(~ind_nonzero) = 0;   % only use pixels within the center disk
    % Filter movie for max projection
    framesCat_conv{i} = imfilter(framesCat{i},psf,'symmetric');
    framesCat4_conv{i} = imfilter(framesCat4_mean{i},psf,'symmetric');
    
%     framesCat_conv{i} = (convn(framesCat{i}, single(reshape([1 1 1] /10, 1, 1, [])), 'same'));
     % Take median of movie
     framesCat_min{i} = median(framesCat_conv{i},3);
     framesCat4_min{i} = median(framesCat4_conv{i},3);
     % Subtract median
     framesCat_medNorm{i} = framesCat_conv{i}-framesCat_min{i};
     framesCat4_medNorm{i} = framesCat4_conv{i}-framesCat4_min{i};
    % take max
    Ymax{i} = max(framesCat_medNorm{i},[],3);
    Ymax{i} = imresize(Ymax{i},scaling);
    framesCat4_Ymax{i} = max(framesCat4_medNorm{i},[],3);
    framesCat4_Ymax{i} = imresize(framesCat4_medNorm{i},scaling);
    
    subplot(4,3,3+i)
    imagesc(Ymax{i});
    colormap('gray');
    set(gca,'YDir','normal');
    hold on;
    title([batName{i} ' ' days1to3{i} ' ' sessionType{i}]);
    
    %filter movie for PNR
    [metadata] = ImBat_defaults;
    [Cn{i}, PNR{i}, PNR_mov{i}] = ImBat_correlation_image(framesCat_medNorm{i},metadata);
    [Cn4{i}, PNR4{i}, PNR_mov4{i}] = ImBat_correlation_image(framesCat4_mean{i},metadata);
    subplot(4,3,i+6)
    imagesc(framesCat4_Ymax{i});
    colormap('gray');
    set(gca,'YDir','normal');
    hold on;
    title(['PNR ' batName{i} ' ' days1to3{i} ' ' sessionType{i}]);
    
end
%% align and overlay
%align all 3 days
[Ymax_aligned{1},Ymax_aligned{2},Ymax_aligned{3}] = ImBat_imageAlign(Ymax{1},Ymax{2},Ymax{3});
[PNR_aligned{1},PNR_aligned{2},PNR_aligned{3}] = ImBat_imageAlign(PNR{1},PNR{2},PNR{3});


%change into RGB and overlay colors
[RGB1,RGB2] = CaBMI_XMASS(Ymax{1},Ymax{2},Ymax{3});
[RGB3,RGB4] = CaBMI_XMASS(Ymax_aligned{1},Ymax_aligned{2},Ymax_aligned{3});
[RGB_PNR1,RGB_PNR2] = CaBMI_XMASS(PNR_aligned{1},PNR_aligned{2},PNR_aligned{3});

subplot(4,3,10)
image((RGB1(:,:,:)));
set(gca,'YDir','normal');
title([batName{1} ': ' days1to3{1} ' ' days1to3{2} ' ' days1to3{3}]);

subplot(4,3,11)
image((RGB3(:,:,:)));
set(gca,'YDir','normal');
title(['Aligned ' batName{1} ': ' days1to3{1} ' ' days1to3{2} ' ' days1to3{3}]);

subplot(4,3,12)
image((RGB_PNR1(:,:,:)));
set(gca,'YDir','normal');
title(['Aligned PNR ' batName{1} ': ' days1to3{1} ' ' days1to3{2} ' ' days1to3{3}]);

%% add to final data structure

framesData.framesLength = framesLength;
framesData.framesMax = framesMax;
framesData.frames = frames;
framesData.framesCat = framesCat;
framesData.framesCat4 = framesCat4;
framesData.Ymax = Ymax; 
framesData.framesCat4_Ymax = framesCat4_Ymax;
framesData.randInd = randInd;
framesData.flightNums = flightNums;
framesData.Ymax_aligned = Ymax_aligned;
framesData.PNR_aligned = PNR_aligned;
framesData.RGB1 = RGB1;
framesData.RGB3 = RGB3;
framesData.RGB_PNR1 = RGB_PNR1;
framesData.Cn = Cn;
framesData.PNR = PNR;
framesData.PNR_mov = PNR_mov;
framesData.Cn4 = Cn4;
framesData.PNR4 = PNR4;
framesData.PNR_mov4 = PNR_mov4;






