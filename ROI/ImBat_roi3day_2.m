function [framesData] = ImBat_roi3day_2(days1to3,varargin)
batName = [];
dateSesh = [];
sessionType = [];
loadFlag = 0;
minLimMult = 3; %min limit multiplier for max projections
maxLimMult = 0.6; %max limit multiplier for max projections
scaling = 4; %scaling for resizing the pixel/frame

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
%% concatenate, filter, max project, and mean 
figure();
sgtitle([batName{1} ': cluster 2']);
scaling = 4; %scaling for resizing the pixel/frame
%initialize frame vectors for video chunks
frames = cell(length(track),1);
framesCat = cell(length(track),1);
framesCat4 = cell(length(track),1);
framesCat4_conv = cell(length(track),1);
framesCat4_minNorm = cell(length(track),1);
framesCat4_min = cell(length(track),1);
Ymax = cell(length(track),1);
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
    %max number of frames across all flights and days
    framesMaxAll = max(framesMax);
    %make gaussian filter based on normcorre function
    gSig = 2; 
    gSiz = 5*gSig; 
    psf = fspecial('gaussian', round(2*gSiz), gSig);
    ind_nonzero = (psf(:)>=max(psf(:,1)));
    psf = psf-mean(psf(ind_nonzero));
    psf(~ind_nonzero) = 0;   % only use pixels within the center disk
    %filt_rad = 20;
    %filt_alpha = 2;
    %h=fspecial('gaussian',filt_rad,filt_alpha);
for i = 1:length(track)
    %filter whole day for full day stability max projection
    movFull_conv{i} = imfilter(video(i).Y,psf,'symmetric'); %filter
    movFull_min{i} = min(movFull_conv{i},[],3); %take min
    movFull_minNorm{i} = movFull_conv{i} - movFull_min{i}; %normalize by subtracting min
    movFull_Ymax{i} = max(movFull_minNorm{i},[],3); %take max
    movFull_Ymax{i} = imresize(movFull_Ymax{i},scaling); %resize the max   
    
    for p = 1:minFlightNums %length(randInd)
        framesDiff{i}(p) = framesMaxAll - framesLength{i}(p);
        %add to frames vector and concatenate
        frames{i}(:,:,:,p) = video(i).Y(:,:,closestIndexStart{i}(p):(closestIndexEnd{i}(p)+framesDiff{i}(p)));
        %max(max) concatenate all frames together in 1 long vector
        framesCat{i} = cat(3,framesCat{i},frames{i}(:,:,:,p));        
        %max(mean)
        framesCat4{i} = frames{i}; %cat(4,framesCat4{i},frames{i}{p}); )
        % Filter movie for max projection
        framesCat4_conv{i}(:,:,:,p) = imfilter(framesCat4{i}(:,:,:,p),psf,'symmetric');
        %framesCat4_conv{i}(:,:,:,p)=imfilter(framesCat4{i}(:,:,:,p),h,'circular');
        % Take median of movie
        framesCat4_min{i}(:,:,:,p) = min(framesCat4_conv{i}(:,:,:,p),[],3);
        % Subtract median
        framesCat4_minNorm{i}(:,:,:,p) = framesCat4_conv{i}(:,:,:,p)-framesCat4_min{i}(:,:,:,p);
        % take max of 1 trial
        Ymax_temp{i}(:,:,p) = max(framesCat4_minNorm{i}(:,:,:,p),[],3);
        Ymax{i}(:,:,p) = imresize(Ymax_temp{i}(:,:,p),scaling);
        
        %make pnr
        [metadata] = ImBat_defaults;
        [Cn4{i}(:,:,p), PNR4{i}(:,:,p), PNR_mov4{i}(:,:,:,p)] = ImBat_correlation_image(framesCat4{i}(:,:,:,p),metadata);      
        
        % plot flights
        subplot(7,3,i);
        plot3(track(i).flightPaths.pos(1,:,track(i).flightPaths.clusterIndex{2}(randInd{i}(p))),track(i).flightPaths.pos(2,:,track(i).flightPaths.clusterIndex{2}(randInd{i}(p))),track(i).flightPaths.pos(3,:,track(i).flightPaths.clusterIndex{2}(randInd{i}(p))));
        hold on;
    end
    %take average across all trials
    framesCat4_mean{i} = mean(framesCat4_minNorm{i},4); %mean filtered video of all trials for each day 
    Ymax_mean{i} = mean(Ymax{i},3); %mean of all the max projections for each trial
    PNR_mean{i} = mean(PNR4{i},3);
    
    % Filter movie for max projection with all trials concatenated
    % together, then filtered and averaged (max(max))
    framesCat3_conv{i} = imfilter(framesCat{i},psf,'symmetric');
    % Take median of movie
    framesCat3_med{i} = min(framesCat3_conv{i},[],3);
    % Subtract median
    framesCat3_medNorm{i} = framesCat3_conv{i}-framesCat3_med{i};
    % take max across the whole trial
    Ymax3{i} = max(framesCat3_medNorm{i},[],3); %max(max)
    Ymax3_resize{i} = imresize(Ymax3{i},scaling);
%% time plots
    filt2use = 5;
    exp2use = 2;
    [im1_rgb{i} norm_max_proj{i},I{i},idx_img{i}] = CABMI_allpxs(framesCat4_mean{i},'filt_rad',filt2use,'exp',exp2use);
    subplot(7,3,i+15)
    imagesc(I{i});
    set(gca,'YDir','normal');
    colorbar;
    title(['Across time ' batName{i} ' ' days1to3{i} ' ' sessionType{1}]);

%% max projection plots
    %full day max projection
    subplot(7,3,3+i)
    imagesc(movFull_Ymax{i},[min(min(movFull_Ymax{1}))*minLimMult max(max(movFull_Ymax{i}))*maxLimMult]);
    colormap('gray');
    set(gca,'YDir','normal');
    hold on;
    title(['Max Full Day ' batName{i} ' ' days1to3{i} ' ' sessionType{i}]);
    colorbar;
    
    %max projection with data filter>max project>mean
    subplot(7,3,6+i)
    imagesc(Ymax_mean{i},[min(min(Ymax_mean{1}))*minLimMult max(max(Ymax_mean{1}))*maxLimMult]);%[0.75 1.9]);
    colormap('gray');
    set(gca,'YDir','normal');
    hold on;
    title(['Max(mean) ' batName{i} ' ' days1to3{i} ' ' sessionType{i}]);
    colorbar;
    
    %PNR data PNR>mean
    subplot(7,3,9+i)
    %imagesc(framesCat4_Ymax{i});
    imagesc(PNR_mean{i},[min(min(PNR_mean{1}))*minLimMult max(max(PNR_mean{1}))*maxLimMult]);
    colormap('gray');
    set(gca,'YDir','normal');
    hold on;
    title(['PNR ' batName{i} ' ' days1to3{i} ' ' sessionType{1}]);
    colorbar;
    
    %max projection with data concat>filter>maxProject
    subplot(7,3,12+i)
    %imagesc(framesCat4_Ymax{i});
    imagesc(Ymax3_resize{i},[min(min(Ymax3_resize{1}))*minLimMult max(max(Ymax3_resize{1}))*maxLimMult]);
    colormap('gray');
    set(gca,'YDir','normal');
    hold on;
    title(['Max(max) ' batName{i} ' ' days1to3{i} ' ' sessionType{i}]);
    colorbar;
    hold off;
end
%% align and overlay 3 days
%align all 3 days
[Ymax_aligned{1},Ymax_aligned{2},Ymax_aligned{3}] = ImBat_imageAlign(Ymax3_resize{1},Ymax3_resize{2},Ymax3_resize{3});
[PNR_aligned{1},PNR_aligned{2},PNR_aligned{3}] = ImBat_imageAlign(PNR_mean{1},PNR_mean{2},PNR_mean{3});
%change into RGB and overlay colors
[RGB1,RGB2] = CaBMI_XMASS(Ymax3_resize{1},Ymax3_resize{2},Ymax3_resize{3},'normalize',2,'hl',[0.05,0.6]);
[RGB3,RGB4] = CaBMI_XMASS(Ymax_aligned{1},Ymax_aligned{2},Ymax_aligned{3},'normalize',2,'hl',[0.05,0.6]);
[RGB_PNR1,RGB_PNR2] = CaBMI_XMASS(PNR_aligned{1},PNR_aligned{2},PNR_aligned{3},'normalize',2,'hl',[0.05,0.6]);

%RGB of nonaligned max projections 
subplot(7,3,19)
image((RGB1(:,:,:)));
set(gca,'YDir','normal');
title(['Max(max) ' batName{1} ': ' days1to3{1} ' ' days1to3{2} ' ' days1to3{3}]);
colorbar;
%RGB of aligned max projections
subplot(7,3,20)
image((RGB3(:,:,:)));
set(gca,'YDir','normal');
title(['Aligned Max(max) ' batName{1} ': ' days1to3{1} ' ' days1to3{2} ' ' days1to3{3}]);
colorbar;
%RGB of aligned PNR
subplot(7,3,21)
image((RGB_PNR1(:,:,:)));
set(gca,'YDir','normal');
title(['Aligned PNR ' batName{1} ': ' days1to3{1} ' ' days1to3{2} ' ' days1to3{3}]);
colorbar;


%% add to final data structure

framesData.framesLength = framesLength;
framesData.framesMax = framesMax;
framesData.frames = frames;
framesData.framesCat = framesCat;
framesData.framesCat4 = framesCat4;
framesData.Ymax_temp = Ymax_temp; 
framesData.Ymax = Ymax;
framesData.randInd = randInd;
framesData.flightNums = flightNums;
framesData.Ymax_aligned = Ymax_aligned;
framesData.PNR_aligned = PNR_aligned;
framesData.RGB1 = RGB1;
framesData.RGB3 = RGB3;
framesData.RGB_PNR1 = RGB_PNR1;
%framesData.Cn = Cn;
%framesData.PNR = PNR;
%framesData.PNR_mov = PNR_mov;
framesData.Cn4 = Cn4;
framesData.PNR4 = PNR4;
framesData.PNR_mov4 = PNR_mov4;






