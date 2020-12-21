function ImBat_RGB_flightAlign(batId,fullSeshTag,day1,day2,day3)
%batId = 'Gal';
%fullSeshTag = 1, chooses to look at whole session for 3 day comparison
clustNum = 2;
saveFlag = 1;
rigidFlag = 1;
hFilt1 = 50; %this is for smoothing the image to determine background for dividing it out to normalize the pixel intensities on all layers
hFilt2 = 25; %this is for smoothing the actual image before the exponential to show the layers, use 25 for rgb overlay and 50 for the difference heat plots
clipRangeShow = [-0.15 0.15]; % [-0.4 0.4] clipping HL to use for the heat map plots
clipRangeOverlap = [0.05 .20]; %clipping HL to use for the RGB overlay
saveTag = ['rigid ' num2str(hFilt2) 'filt'];
dirAllTrials = pwd;
h1 = fspecial('disk',hFilt1); %normalize background smoothing
h2 = fspecial('disk',hFilt2); %smooth presentation image

%load first day placeCellStableROI data
if strcmp(batId,'Gal')
    %cd([dirTop(day1).folder filesep 'plots\200911-preDurPost cells across days']);
    load('Gal_200311to200324_activity_allTrials_allClusts_allTrials_sMat_newDff_newOrder.mat'); %load activity for pre,dur,post
elseif strcmp(batId,'Gen')
    load('Gen_200319to200324_activity_allTrials_allClusts_sMat_newDff_newOrder.mat'); %load activity for pre,dur,post
elseif strcmp(batId,'z2')
    load('Z2_200701to200712_activity_allTrials_allClusts_sMat_newDff_newOrder.mat');
end
%make saving directory
if saveFlag == 1
    saveDir1 = '/Volumes/Tobias_flig/topQualityData/analysis_done/plots/';
    %saveDir1 = '\\169.229.54.11\server_home\users\tobias\flight\data_processed\topQualityData\analysis_done\plots\';
    %saveDir1 = '/Users/periscope/Desktop/analysis/flight/plots/';
    %saveDir1 = 'C:\Users\tobias\Desktop\analysis\plots\';
    if ~exist([saveDir1 datestr(now,'yymmdd') filesep 'maxProjFlightAlign'])
        mkdir([saveDir1 datestr(now,'yymmdd') filesep 'maxProjFlightAlign']);
    else
        disp('You have been working today...');
    end
    saveDir = [saveDir1 datestr(now,'yymmdd') filesep 'maxProjFlightAlign' filesep];
end
%% pull out the correct data into the raw image variables
if fullSeshTag == 0 %comparing flight aligned max proj
    if isempty(day3)
        if day1 == day2 %if looking at the cluster flight aligned vs full session
            %set raw to the flight aligned frames and full session
            IM1_raw = activity_allTrials.maxMeanFrames_dur{clustNum}{day1};
            IM2_raw = activity_allTrials.YmaxFull{day1};
            IM3_raw = activity_allTrials.YmaxFull{day1}; %duplicate for now
            plotTitle = [saveTag ' ' batId ' clust ' num2str(clustNum) ': Day ' num2str(day1) ' (r) v full session (c)'];
        else %if looking at 1 day vs another day only
            %set raw to the flight aligned frames,
            IM1_raw = activity_allTrials.maxMeanFrames_dur{clustNum}{day1};
            IM2_raw = activity_allTrials.maxMeanFrames_dur{clustNum}{day2};
            IM3_raw = activity_allTrials.maxMeanFrames_dur{clustNum}{day2}; %duplicate for now
            plotTitle = [saveTag ' ' batId ' clust ' num2str(clustNum) ': Day ' num2str(day1) ' (r) v Day ' num2str(day2) ' (c)']
        end
    else %if looking at comparing 3 days
        IM1_raw = activity_allTrials.maxMeanFrames_dur{clustNum}{day1};
        IM2_raw = activity_allTrials.maxMeanFrames_dur{clustNum}{day2};
        IM3_raw = activity_allTrials.maxMeanFrames_dur{clustNum}{day3};
        plotTitle = [saveTag ' ' batId ' clust ' num2str(clustNum) ': Day ' num2str(day1) ' (r) v Day ' num2str(day2) ' (g) v Day ' num2str(day3) ' (b)'];
    end
else %comparing the full session
    if isempty(day3) %comparing 1 day vs another
        IM1_raw = activity_allTrials.YmaxFull{day1};
        IM2_raw = activity_allTrials.YmaxFull{day2};
        IM3_raw = activity_allTrials.YmaxFull{day2}; %duplicate for now
        plotTitle = [saveTag ' ' batId ' clust ' num2str(clustNum) ' Full Sesh: Day ' num2str(day1) ' (r) v Day ' num2str(day2) ' (c)'];
    else %comparing 3 days full
        IM1_raw = activity_allTrials.YmaxFull{day1};
        IM2_raw = activity_allTrials.YmaxFull{day2};
        IM3_raw = activity_allTrials.YmaxFull{day3};
        plotTitle = [saveTag ' ' batId ' Full Sesh: Day ' num2str(day1) ' (r) v Day ' num2str(day2) ' (g) v Day ' num2str(day3) ' (b)'];
    end
end

% IM_rawall = cell(3,1);
% for day_i = 1:3
%     IM_rawall{day_i} = zeros(size(YmaxFull{day_i},1),size(YmaxFull{day_i},2),3);
% for part_i = 1:3
% IM_rawall{day_i}(:,:,part_i) = YmaxFull{part_i+(day_i-1)*3};
% end
% end
% IM1_raw = mean(IM_rawall{1},3);
% IM2_raw = mean(IM_rawall{2},3);
% IM3_raw = mean(IM_rawall{3},3);
% IM1_raw = IM_rawall{3}(:,:,1);
% IM2_raw = IM_rawall{3}(:,:,2);
% IM3_raw = IM_rawall{3}(:,:,3);

%filter, and divide background to make the pixels all the same brightness
%resize and subtract background
IM1doub = imresize(double(IM1_raw),2);
IM2doub = imresize(double(IM2_raw),2);
IM3doub = imresize(double(IM3_raw),2);
% IM1_filt = IM1doub;
% IM2_filt = IM2doub;
% IM3_filt = IM3doub;

bground1=imfilter(IM1doub,h1,'replicate');%smoothdata(IM1doub,'gaussian',2);%
bground2=imfilter(IM2doub,h1,'replicate');%smoothdata(IM2doub,'gaussian',2);%
bground3=imfilter(IM3doub,h1,'replicate');%smoothdata(IM3doub,'gaussian',2);%
IM1_filt=IM1doub./(bground1+5);
IM2_filt=IM2doub./(bground2+5);
IM3_filt=IM3doub./(bground3+5);

%concatenate to make into grayscale and then break apart into correct sizes
IMcat = [IM1_filt IM2_filt IM3_filt];    IMcat_gray = mat2gray(IMcat);
IM1 = IMcat_gray(:,1:size(IMcat_gray,2)/3);
IM2 = IMcat_gray(:,(size(IMcat_gray,2)/3)+1:2*size(IMcat_gray,2)/3);
IM3 = IMcat_gray(:,(2*size(IMcat_gray,2)/3)+1:end);
%make sure all images are same size if overlapping later
minRow = min([size(IM1,1),size(IM2,1),size(IM3,1)]);
minCol = min([size(IM1,2),size(IM2,2),size(IM3,2)]);

%make copy of the image to heavily filter for image registration
image1 = IM1;
image2 = IM2;
image3 = IM3;
% imfilt1 = imfilter(image1,h,'replicate');
% imfilt2 = imfilter(image2,h,'replicate');
% imfilt3 = imfilter(image3,h,'replicate');
%  image1 = imfilt1.^2;
%  image2 = imfilt2.^2;
%  image3 = imfilt3.^2;
% Clip image ( binarize) above 25 percentile
image1(image1<0.25) =0;
image2(image2<0.25) =0;
image3(image3<0.25) =0;
% Clip image ( binarize) above 25 percentile
image1(image1>=0.75) =1;
image2(image2>=0.75) =1;
image3(image3>=0.75) =1;
% remove edges ( in case these are producing artifacts)
edgeCutoff = 70;
image1(1:edgeCutoff,:) = 0;
image1(:,1:edgeCutoff) = 0;
image1(end-edgeCutoff:end,:) = 0;
image1(:,end-edgeCutoff:end) = 0;
image2(1:edgeCutoff,:) = 0;
image2(:,1:edgeCutoff) = 0;
image2(end-edgeCutoff:end,:) = 0;
image2(:,end-edgeCutoff:end) = 0;
image3(1:edgeCutoff,:) = 0;
image3(:,1:edgeCutoff) = 0;
image3(end-edgeCutoff:end,:) = 0;
image3(:,end-edgeCutoff:end) = 0;

if isempty(day3)
    IM3 = [];
    image3 = [];
end
if rigidFlag == 1 %use the rigid alignment with affine
    try
        [imagesAligned] = ImBat_imageAlign(image1,image2,image3,IM1(1:minRow,1:minCol),IM2(1:minRow,1:minCol),IM3(1:minRow,1:minCol));
    catch
        [imagesAligned] = ImBat_imageAlign(image1,image2,image3,IM1(1:minRow,1:minCol),IM2(1:minRow,1:minCol),IM3);
    end
else %use nonrigid with demon flow
    [figNonrigid,imagesAligned] = ImBat_imageAlign_nonrigid(image1,image2,image3,IM1,IM2,IM3);
end
IM1_aligned = imagesAligned.IM1_aligned;
IM2_aligned = imagesAligned.IM2_aligned;
IM3_aligned = imagesAligned.IM3_aligned;

% IM1doub = imresize(double(IM1_aligned),2);
% IM2doub = imresize(double(IM2_aligned),2);
% IM3doub = imresize(double(IM3_aligned),2);
% bground1=imfilter(IM1doub,h,'replicate');%smoothdata(IM1doub,'gaussian',2);%
% bground2=imfilter(IM2doub,h,'replicate');%smoothdata(IM2doub,'gaussian',2);%
% bground3=imfilter(IM3doub,h,'replicate');%smoothdata(IM3doub,'gaussian',2);%
% IM1_bfilt=imresize(double(IM1_aligned),2)./(bground1+5);
% IM2_bfilt=imresize(double(IM2_aligned),2)./(bground2+5);
% IM3_bfilt=imresize(double(IM3_aligned),2)./(bground3+5);

%filter, then average the images and subtract each from that average to show difference
IM1_aligned_filt = imfilter(IM1_aligned,h2,'replicate');%(IM1_aligned,hFilt,'replicate');
IM2_aligned_filt = imfilter(IM2_aligned,h2,'replicate');%(IM2_aligned,hFilt,'replicate');
IM3_aligned_filt = imfilter(IM3_aligned,h2,'replicate');%(IM3_aligned,hFilt,'replicate');
IM1_aligned_filt = IM1_aligned_filt.^3;
IM2_aligned_filt = IM2_aligned_filt.^3;
IM3_aligned_filt = IM3_aligned_filt.^3;

%RGB the images
[aOverlap,bOverlap] = CaBMI_XMASS(IM1_aligned_filt,IM2_aligned_filt,IM3_aligned_filt,'hl',clipRangeOverlap);
%rgb the original unaligned for comparison
[aUnaligned,bUnaligned] = CaBMI_XMASS(IM1(1:minRow,1:minCol),IM2(1:minRow,1:minCol),IM3);

%quantify quartiles of overlap to see how good the days overlap
sizeIm = size(IM1_aligned_filt)/2;
quartRange(1,:) = [1,sizeIm(1),1,sizeIm(2)];
quartRange(2,:) = [1,sizeIm(1),sizeIm(2)+1,sizeIm(2)*2];
quartRange(3,:) = [sizeIm(1)+1,sizeIm(1)*2,1,sizeIm(2)];
quartRange(4,:) = [sizeIm(1)+1,sizeIm(1)*2,sizeIm(2)+1,sizeIm(2)*2];

for quart_i = 1:4
    %quartRange = [1+sizeIm(1)*(quart_i-1),quart_i*sizeIm(1),1+sizeIm(2)*(quart_i-1),quart_i*sizeIm(2)];
overlapQuart(quart_i,1) = sum(sum(abs(IM1_aligned_filt(quartRange(quart_i,1):quartRange(quart_i,2),quartRange(quart_i,3):quartRange(quart_i,4))-IM2_aligned_filt(quartRange(quart_i,1):quartRange(quart_i,2),quartRange(quart_i,3):quartRange(quart_i,4)))));
overlapQuart(quart_i,2) = sum(sum(abs(IM2_aligned_filt(quartRange(quart_i,1):quartRange(quart_i,2),quartRange(quart_i,3):quartRange(quart_i,4))-IM3_aligned_filt(quartRange(quart_i,1):quartRange(quart_i,2),quartRange(quart_i,3):quartRange(quart_i,4)))));
%overlapQuart(quart_i,3) = sum(sum(abs(IM1_aligned_filt(quartRange(quart_i,1):quartRange(quart_i,2),quartRange(quart_i,3):quartRange(quart_i,4))-IM3_aligned_filt(quartRange(quart_i,1):quartRange(quart_i,2),quartRange(quart_i,3):quartRange(quart_i,4)))));
medOverlapQuart(quart_i) = median(overlapQuart(quart_i,:));
end
medOverlapAll = round(median(overlapQuart,'all'));

%use the larger filter of 50 to smooth the images for subtraction/comparisons
IM1_aligned_filt_comp = imfilter(IM1_aligned,h1,'replicate');%(IM1_aligned,hFilt,'replicate');
IM2_aligned_filt_comp = imfilter(IM2_aligned,h1,'replicate');%(IM2_aligned,hFilt,'replicate');
IM3_aligned_filt_comp = imfilter(IM3_aligned,h1,'replicate');%(IM3_aligned,hFilt,'replicate');
IM1_aligned_filt_comp = IM1_aligned_filt_comp.^3;
IM2_aligned_filt_comp = IM2_aligned_filt_comp.^3;
IM3_aligned_filt_comp = IM3_aligned_filt_comp.^3;
IM_sum = IM1_aligned_filt_comp + IM2_aligned_filt_comp + IM3_aligned_filt_comp;
IM_mean = IM_sum/3;
IM1_diff = IM1_aligned_filt_comp - IM_mean;
IM2_diff = IM2_aligned_filt_comp - IM_mean;
IM3_diff = IM3_aligned_filt_comp - IM_mean;
IM_1diff2 = IM2_aligned_filt_comp - IM1_aligned_filt_comp;
%IM_1diff2 = IM_1diff2.^5;
IM_2diff3 = IM3_aligned_filt_comp - IM2_aligned_filt_comp;
%IM_2diff3 = IM_2diff3.*5;
IM_1diff3 = IM3_aligned_filt_comp - IM1_aligned_filt_comp;
%IM_1diff3 = IM_1diff3.*5;
%concatenate to take min/max and subtract from that
IM_combined(:,:,1) = IM1_aligned_filt_comp;
IM_combined(:,:,2) = IM2_aligned_filt_comp;
IM_combined(:,:,3) = IM3_aligned_filt_comp;
IM_combined_max = max(IM_combined,[],3);
IM_combined_min = min(IM_combined,[],3);
IM1_maxDiff = IM1_aligned_filt_comp - IM_combined_max;
IM2_maxDiff = IM2_aligned_filt_comp - IM_combined_max;
IM3_maxDiff = IM3_aligned_filt_comp - IM_combined_max;
IM1_minDiff = IM1_aligned_filt_comp - IM_combined_min;
IM2_minDiff = IM2_aligned_filt_comp - IM_combined_min;
IM3_minDiff = IM3_aligned_filt_comp - IM_combined_min;

%% load masks and warp to match registered max projection for each day
if rigidFlag == 1 %only if have the rigid twarp output
    maskDir1 = '\\169.229.54.11\server_home\users\tobias\flight\data_processed\topQualityData\analysis_done\new3\';
    maskDir1 = '/Volumes/Tobias_flig/topQualityData/analysis_done';
    days = [day1 day2 day3];
    ROI2plot = cell(length(days),1);
    ROI_coords = cell(length(days),1);
    centroid = cell(length(days),1);
    col = cell(length(days));
    for day_i = 1:length(days)
        %cd([maskDir1 batId]);
        cd([maskDir1]);
        dayDir = dir([batId(1:2) '*']);
        cd([dayDir(days(day_i)).name filesep 'extracted']);
        flyDir = dir('*fly-*extraction');
        cd(flyDir(end).name);
        processedDir = dir('processed*');
        load([processedDir(end).name filesep 'results.mat']);
        Atemp = full(results.A);
        resultsA{day_i} = Atemp;
        resultsCn{day_i} = results.Cn;
        resultsCraw{day_i} = results.C_raw;
        
        % Plot binary mask of all neurons in the A matrix
        %convert A matrix into full matrix
        Ysiz = size(resultsCn{day_i});
        nRois = size(resultsA{day_i},2);
        scaling = 10; %depends on size of frame and downsampling from extraction step
        %ROI2plot = (:,:,zeros(length(Atemp(1,:)));
        % get ROI centroids for top 30%;
        %ROI2plot{day_i} = zeros(Ysiz(1),Ysiz(2),nRois);
        %centroid{day_i} = zeros(nRois,2);
        col{day_i} = zeros(nRois,3);
        for roi_i = 1:nRois
            %create 3d matrix with all ROI heat maps
            ROI2plot{day_i}(:,:,roi_i) = mat2gray(reshape(Atemp(:,roi_i),Ysiz(1),Ysiz(2)));
            ROI2plot_scaled{day_i}(:,:,roi_i) = imresize(ROI2plot{day_i}(:,:,roi_i),scaling);
        end
        cd(dirAllTrials); %go back to home activitity all_trials directory
    end
    %register day ROIs 1 to day 2 alignment
    for day_i = 1:length(days)
        for roi_i = 1:size(resultsA{day_i},2)
            if day_i == 1
                ROI2plot_aligned{1}(:,:,roi_i) = imwarp(ROI2plot_scaled{1}(:,:,roi_i),imagesAligned.IM1_tform,'OutputView',imref2d(size(IM2)));
            elseif day_i == 2
                ROI2plot_aligned{2}(:,:,roi_i) = ROI2plot_scaled{2}(:,:,roi_i);
            elseif day_i == 3
                ROI2plot_aligned{3}(:,:,roi_i) = imwarp(ROI2plot_scaled{3}(:,:,roi_i),imagesAligned.IM3_tform,'OutputView',imref2d(size(IM2)));
            end
            binaryImage = imbinarize(ROI2plot_aligned{day_i}(:,:,roi_i));
            [y,x]=find(binaryImage);
            %get ROI coordinates
            ROI_coords{day_i}(roi_i,1) = {x};
            ROI_coords{day_i}(roi_i,2) = {y};
            %calculate centroids
            centroid{day_i}(roi_i,1) = mean(cell2mat(ROI_coords{day_i}(roi_i,1)));%*scaling;
            centroid{day_i}(roi_i,2) = mean(cell2mat(ROI_coords{day_i}(roi_i,2)));%*scaling;%get the centroid of mask
            %subplot(2,1,1);
            %imshow(ROI2plot_aligned{1}(:,:,roi_i));
            %subplot(2,1,2);
            %p = plot(ROI_coords{1}{roi_i,1},ROI_coords{1}{roi_i,2},'LineWidth',4);
            %clf;
        end
        col{day_i} = jet(size(resultsA{day_i},2));
    end
    %ROI_coords = smoothdata(ROI_coords,'gaussian',3); %filter the ROI coordinate mask so it is not so jagged
end






%% plots
%plot the unaligned vs aligned for comparison
plotDay123_unaligned = figure();
sgtitle(plotTitle);
subplot(1,2,1);
image((aUnaligned(:,:,:)));
title('un-aligned');
set(gca,'xticklabel',[],'yticklabel',[]);
subplot(1,2,2);
image((aOverlap(:,:,:)));
title('aligned');
set(gca,'xticklabel',[],'yticklabel',[]);

%plot just the aligned figure with all 3rgb overlapping
plotDay123_aligned = figure();
image(aOverlap(:,:,:));
title([plotTitle ' medDiff ' num2str(medOverlapAll)]);
set(gca,'xticklabel',[],'yticklabel',[]);
%set(gca,'YDir','normal');

%plot aligned RGB figure and each day with centroids on top
plotDay123_centroids = figure('units','normalized','outerposition',[0 0 1 0.6]);
sgtitle([plotTitle ' and daily centroids']);
ha = tight_subplot(1,4,[.02 .01],[.01 .08],[.01 .01]);
set(0,'CurrentFigure',plotDay123_centroids);
axes(ha(1)); imshow((aOverlap(:,:,:)),[]); title('RGB Overlap');
axes(ha(2)); imshow(IM1_aligned_filt,[]);
colormap(ha(2),gray);
hold on;
for roi_i = 1:size(ROI2plot_aligned{1},3)
        plot(ROI_coords{1}{roi_i,1},ROI_coords{1}{roi_i,2},'LineWidth',4,'Color',[col{1}(roi_i,:) 0.01]);
        p = text(centroid{1}(roi_i,1),centroid{1}(roi_i,2),num2str(roi_i));
        p.Color(1:3) = col{1}(size(ROI2plot_aligned{1},3) + 1 - roi_i,:);
end
title('Day 1 ROIs');
hold off
axes(ha(3)); imshow(IM2_aligned_filt,[]);
colormap(ha(3),gray);
hold on;
for roi_i = 1:size(ROI2plot_aligned{2},3)
    plot(ROI_coords{2}{roi_i,1},ROI_coords{2}{roi_i,2},'LineWidth',4,'Color',[col{2}(roi_i,:) 0.01]);
    p = text(centroid{2}(roi_i,1),centroid{2}(roi_i,2),num2str(roi_i));
    p.Color(1:3) = col{2}(size(ROI2plot_aligned{2},3) + 1 - roi_i,:);
end
title('Day 2 ROIs');
hold off
axes(ha(4)); imshow(IM3_aligned_filt,[]);
colormap(ha(4),gray);
hold on;
try
for roi_i = 1:size(ROI2plot_aligned{3},3)
        plot(ROI_coords{3}{roi_i,1},ROI_coords{3}{roi_i,2},'LineWidth',4,'Color',[col{3}(roi_i,:) 0.01]);
        p = text(centroid{3}(roi_i,1),centroid{3}(roi_i,2),num2str(roi_i));
        p.Color(1:3) = col{3}(size(ROI2plot_aligned{3},3) + 1 - roi_i,:);
end
catch
end
title('Day 3 ROIs');
hold off

%plot the rgb, mean, and each difference
plotDay123_diffs = figure('units','normalized','outerposition',[0 0 1 0.6]);
sgtitle([plotTitle ' and daily differences from mean']);
ha = tight_subplot(1,5,[.02 .01],[.01 .08],[.01 .01]);
set(0,'CurrentFigure',plotDay123_diffs);
axes(ha(1)); imshow(aOverlap(:,:,:),[]); title('RGB Overlap');
%colorbar('southoutside');
axes(ha(2)); imshow(IM_mean,[]); title('Mean Day 1, 2, 3');
colormap(ha(2),gray);
colorbar('southoutside');
axes(ha(3)); imshow(IM1_diff,clipRangeShow); title('Diff Day 1');
colormap(ha(3),fireice);
colorbar('southoutside');
axes(ha(4)); imshow(IM2_diff,clipRangeShow); title('Diff Day 2');
colormap(ha(4),fireice);
colorbar('southoutside');
axes(ha(5)); imshow(IM3_diff,clipRangeShow); title('Diff Day 3');
colormap(ha(5),fireice);
colorbar('southoutside');

%plot rgb and each progression of differences from 1 day to next
plotDay123_progression = figure('units','normalized','outerposition',[0 0 1 0.6]);
sgtitle([plotTitle ' and day to day progression']);
ha = tight_subplot(1,5,[.02 .01],[.01 .08],[.01 .01]);
set(0,'CurrentFigure',plotDay123_progression);
axes(ha(1)); imshow(aOverlap(:,:,:),[]); title('RGB Overlap');
%colorbar('southoutside');
axes(ha(2)); imshow(IM1_aligned_filt,[]); title('Day 1');
colormap(ha(2),gray);
colorbar('southoutside');
axes(ha(3)); imshow(IM_1diff2,clipRangeShow); title('Day 2 - Day 1');
colormap(ha(3),fireice);
colorbar('southoutside');
axes(ha(4)); imshow(IM_2diff3,clipRangeShow); title('Day 3 - Day 2');
colormap(ha(4),fireice);
colorbar('southoutside');
axes(ha(5)); imshow(IM_1diff3,clipRangeShow); title('Day 3 - Day 1');
colorbar('southoutside');
colormap(ha(5),fireice);

%% plot rgb and each progression of differences from 1 day to next
plotDay123_all = figure('units','normalized','outerposition',[0 0 1 1]);
sgtitle([plotTitle ' medDiff ' num2str(medOverlapAll)]);
ha = tight_subplot(4,5,[.02 .01],[.01 .08],[.01 .01]);
set(0,'CurrentFigure',plotDay123_all);
axes(ha(1)); imshow(aOverlap(:,:,:),[]); title('RGB Overlap');
%colorbar('southoutside');
axes(ha(2)); imshow(IM1_aligned_filt,[]); title('Day 1');
%colorbar('southoutside');
axes(ha(3)); imshow(IM_1diff2,clipRangeShow); title('Day 2 - Day 1');
%colorbar('southoutside');
axes(ha(4)); imshow(IM_2diff3,clipRangeShow); title('Day 3 - Day 2');
%colorbar('southoutside');
axes(ha(5)); imshow(IM_1diff3,clipRangeShow); title('Day 3 - Day 1');
%colorbar('southoutside');
colormap(ha(3),fireice);
colormap(ha(4),fireice);
colormap(ha(5),fireice);
axes(ha(6)); imshow(IM1_aligned_filt,[]);
colormap(ha(6),gray);
hold on;
for roi_i = 1:size(ROI2plot_aligned{1},3)
        plot(ROI_coords{1}{roi_i,1},ROI_coords{1}{roi_i,2},'LineWidth',4,'Color',[col{1}(roi_i,:) 0.01]);
        p = text(centroid{1}(roi_i,1),centroid{1}(roi_i,2),num2str(roi_i));
        p.Color(1:3) = col{1}(size(ROI2plot_aligned{1},3) + 1 - roi_i,:);
end
title('Day 1 ROIs');
hold off
axes(ha(7)); imshow(IM_mean,[]); title('Mean Day 1, 2, 3');
%colorbar('southoutside');
axes(ha(8)); imshow(IM1_diff,clipRangeShow); title('Day 1 - mean');
%colorbar('southoutside');
axes(ha(9)); imshow(IM2_diff,clipRangeShow); title('Day 2 - mean');
%colorbar('southoutside');
axes(ha(10)); imshow(IM3_diff,clipRangeShow); title('Day 3 - mean');
colormap(ha(8),fireice);
colormap(ha(9),fireice);
colormap(ha(10),fireice);%colorbar('southoutside');
axes(ha(11)); imshow(IM2_aligned_filt,[]);
colormap(ha(11),gray);
hold on;
for roi_i = 1:size(ROI2plot_aligned{2},3)
    plot(ROI_coords{2}{roi_i,1},ROI_coords{2}{roi_i,2},'LineWidth',4,'Color',[col{2}(roi_i,:) 0.01]);
    p = text(centroid{2}(roi_i,1),centroid{2}(roi_i,2),num2str(roi_i));
    p.Color(1:3) = col{2}(size(ROI2plot_aligned{2},3) + 1 - roi_i,:);
end
title('Day 2 ROIs');
hold off
axes(ha(12)); imshow(IM_combined_max,[]); title('Max Day 1, 2, 3');
%colorbar('southoutside');
axes(ha(13)); imshow(IM1_maxDiff,clipRangeShow); title('Day 1  - max');
%colorbar('southoutside');
axes(ha(14)); imshow(IM2_maxDiff,clipRangeShow); title('Day 2  - max');
%colorbar('southoutside');
axes(ha(15)); imshow(IM3_maxDiff,clipRangeShow); title('Day 3 - max');
colormap(ha(13),fireice);
colormap(ha(14),fireice);
colormap(ha(15),fireice);%colorbar('southoutside');
axes(ha(16)); imshow(IM3_aligned_filt,[]);
colormap(ha(16),gray);
hold on;
try
for roi_i = 1:size(ROI2plot_aligned{3},3)
        plot(ROI_coords{3}{roi_i,1},ROI_coords{3}{roi_i,2},'LineWidth',4,'Color',[col{3}(roi_i,:) 0.01]);
        p = text(centroid{3}(roi_i,1),centroid{3}(roi_i,2),num2str(roi_i));
        p.Color(1:3) = col{3}(size(ROI2plot_aligned{3},3) + 1 - roi_i,:);
end
catch
end
title('Day 3 ROIs');
hold off
axes(ha(17)); imshow(IM_combined_min,[]); title('Min Day 1, 2, 3');
%colorbar('southoutside');
axes(ha(18)); imshow(IM1_minDiff,clipRangeShow); title('Day 1 - min');
%colorbar('southoutside');
axes(ha(19)); imshow(IM2_minDiff,clipRangeShow); title('Day 2 - min');
%colorbar('southoutside');
axes(ha(20)); imshow(IM3_minDiff,clipRangeShow); title('Day 3 - min');
colormap(ha(18),fireice);
colormap(ha(19),fireice);
colormap(ha(20),fireice);%colorbar('southoutside');

%% plot the cnmfe time series only for the 3 days
plotCnmfeSeries = figure('units','normalized','outerposition',[0 0 1 1]);
sgtitle(plotTitle);
% subplot(4,4,1);
% image(aOverlap(:,:,:));
% title([plotTitle ' medDiff ' num2str(medOverlapAll)]);
% set(gca,'xticklabel',[],'yticklabel',[]);

for day_i = 1:length(resultsCraw)
    subplot(3,1,day_i);
for roi_i = 1:size(resultsCraw{day_i},1); 
    plot(1:size(resultsCraw{day_i},2),zscore(smoothdata(resultsCraw{day_i}(roi_i,:),'movmedian',30))+roi_i*6); 
    hold on;
end
dataTicks = [6:6:size(resultsCraw{day_i},1)*6];
set(gca,'YTick',dataTicks,'YTickLabel',[1:length(dataTicks)],'xlim',[0 size(resultsCraw{day_i},2)]);
xt = get(gca,'xtick');
set(gca,'XTick',xt, 'xticklabel',xt/(results.Fs*60));
title(['Day ' num2str(day_i)]);
ylabel([num2str(size(resultsCraw{day_i},1)) ' ROIs']);
if day_i == 3
    xlabel('Time (m)');
end
end

%% plot the time series from cnmfe next to RGB and daily overlap to check that cells are selected properly
%plot rgb overlap
plotOverlapTimeSeries = figure('units','normalized','outerposition',[0 0 1 1]);
sgtitle([plotTitle ' medDiff ' num2str(medOverlapAll)]);
%ha = tight_subplot(4,4,[.02 .01],[.01 .08],[.01 .01]);
set(0,'CurrentFigure',plotOverlapTimeSeries);
ax1 = subplot(4,4,1); imagesc(aOverlap(:,:,:)); title('RGB Overlap');
set(gca,'XTick',[],'YTick',[],'XTickLabel',[],'YTickLabel',[]);
%colorbar('southoutside');
ax2 = subplot(4,4,2); imagesc(IM_1diff2,clipRangeShow); title('Day 2 - Day 1');
set(gca,'XTick',[],'YTick',[],'XTickLabel',[],'YTickLabel',[]);
%colorbar('southoutside');
ax3 = subplot(4,4,3); imagesc(IM_2diff3,clipRangeShow); title('Day 3 - Day 2');
set(gca,'XTick',[],'YTick',[],'XTickLabel',[],'YTickLabel',[]);
%colorbar('southoutside');
ax4 = subplot(4,4,4); imagesc(IM_1diff3,clipRangeShow); title('Day 3 - Day 1');
set(gca,'XTick',[],'YTick',[],'XTickLabel',[],'YTickLabel',[]);
colormap(ax2,fireice);
colormap(ax3,fireice);
colormap(ax4,fireice);
ax5 = subplot(4,4,5); imagesc(IM1_aligned_filt);
colormap(ax5,gray);
set(gca,'XTick',[],'YTick',[],'XTickLabel',[],'YTickLabel',[]);
hold on;
for roi_i = 1:size(ROI2plot_aligned{1},3)
        plot(ROI_coords{1}{roi_i,1},ROI_coords{1}{roi_i,2},'LineWidth',4,'Color',[col{1}(roi_i,:) 0.01]);
        p = text(centroid{1}(roi_i,1),centroid{1}(roi_i,2),num2str(roi_i));
        p.Color(1:3) = col{1}(size(ROI2plot_aligned{1},3) + 1 - roi_i,:);
end
title('Day 1 ROIs');
hold off
subplot(4,4,6:8);
for roi_i = 1:size(resultsCraw{1},1); 
    plot(1:size(resultsCraw{1},2),zscore(smoothdata(resultsCraw{1}(roi_i,:),'movmedian',30))+roi_i*6); 
    hold on;
end
dataTicks = [6:6:size(resultsCraw{1},1)*6];
set(gca,'YTick',dataTicks,'YTickLabel',[1:length(dataTicks)],'xlim',[0 size(resultsCraw{1},2)]);
xt = get(gca,'xtick');
set(gca,'XTick',xt, 'xticklabel',xt/(results.Fs*60));
title(['Day 1']);
ylabel([num2str(size(resultsCraw{1},1)) ' ROIs']);
hold off
ax9 = subplot(4,4,9); imagesc(IM2_aligned_filt);
colormap(ax9,gray);
set(gca,'XTick',[],'YTick',[],'XTickLabel',[],'YTickLabel',[]);
hold on;
for roi_i = 1:size(ROI2plot_aligned{2},3)
        plot(ROI_coords{2}{roi_i,1},ROI_coords{2}{roi_i,2},'LineWidth',4,'Color',[col{2}(roi_i,:) 0.01]);
        p = text(centroid{2}(roi_i,1),centroid{2}(roi_i,2),num2str(roi_i));
        p.Color(1:3) = col{2}(size(ROI2plot_aligned{2},3) + 1 - roi_i,:);
end
title('Day 2 ROIs');
hold off
subplot(4,4,10:12); 
for roi_i = 1:size(resultsCraw{2},1); 
    plot(1:size(resultsCraw{2},2),zscore(smoothdata(resultsCraw{2}(roi_i,:),'movmedian',30))+roi_i*6); 
    hold on;
end
dataTicks = [6:6:size(resultsCraw{2},1)*6];
set(gca,'YTick',dataTicks,'YTickLabel',[1:length(dataTicks)],'xlim',[0 size(resultsCraw{2},2)]);
xt = get(gca,'xtick');
set(gca,'XTick',xt, 'xticklabel',xt/(results.Fs*60));
title(['Day 2']);
ylabel([num2str(size(resultsCraw{2},1)) ' ROIs']);
hold off
ax13 = subplot(4,4,13); imagesc(IM3_aligned_filt);
colormap(ax13,gray);
set(gca,'XTick',[],'YTick',[],'XTickLabel',[],'YTickLabel',[]);
hold on;
try
for roi_i = 1:size(ROI2plot_aligned{3},3)
        plot(ROI_coords{3}{roi_i,1},ROI_coords{3}{roi_i,2},'LineWidth',4,'Color',[col{3}(roi_i,:) 0.01]);
        p = text(centroid{3}(roi_i,1),centroid{3}(roi_i,2),num2str(roi_i));
        p.Color(1:3) = col{3}(size(ROI2plot_aligned{3},3) + 1 - roi_i,:);
end
catch
end
title('Day 3 ROIs');
hold off
subplot(4,4,14:16);
for roi_i = 1:size(resultsCraw{3},1); 
    plot(1:size(resultsCraw{3},2),zscore(smoothdata(resultsCraw{3}(roi_i,:),'movmedian',30))+roi_i*6); 
    hold on;
end
dataTicks = [6:6:size(resultsCraw{3},1)*6];
set(gca,'YTick',dataTicks,'YTickLabel',[1:length(dataTicks)],'xlim',[0 size(resultsCraw{3},2)]);
xt = get(gca,'xtick');
set(gca,'XTick',xt, 'xticklabel',xt/(results.Fs*60));
title(['Day 3']);
ylabel([num2str(size(resultsCraw{3},1)) ' ROIs']);
hold off
%% ask which cells and days times series to compare
inputGo = [];
inputGo = input('Go?');
while isempty(inputGo)
    inputDay = input('What days would you like to compare?');
inputROI = input('What ROIs would you like to compare?');
Ysiz = size(resultsCn{1}); %size of the image frame
Rtime = [];
Rmask = [];
timeseries = [];
binaryMask = [];
binaryMaskAlign = [];
%find duration of each day and concatenate the day/roi numbers for titles 
for comp_i = 1:length(inputDay)
    dur(comp_i) = length(resultsCraw{inputDay(comp_i)});
    dayROI{comp_i} = [num2str(inputDay(comp_i)) '.' num2str(inputROI(comp_i))];
end
minDurDay = min(dur); %find min duration of the days to limit the time series
%extract the time series and binary mask of each Day/ROI combo 
for comp_i = 1:length(inputDay)
    timeseries(:,comp_i) = zscore(smoothdata(resultsCraw{inputDay(comp_i)}(inputROI(comp_i),1:minDurDay),2,'movmedian',30))';
    
    roiMask = imresize(mat2gray(reshape(resultsA{inputDay(comp_i)}(:,inputROI(comp_i)),Ysiz(1),Ysiz(2))),10);
    binaryMask(:,:,comp_i) = imbinarize(roiMask);%mat2gray(reshape(resultsA{inputDay(comp_i)}(:,inputROI(comp_i)),Ysiz(1),Ysiz(2))));
    
    %warp the binary mask according to the twarps from the image alignment
    if inputDay(comp_i) == 1
    binaryMaskAlign(:,:,comp_i) = imwarp(binaryMask(:,:,comp_i),imagesAligned.IM1_tform,'OutputView',imref2d(size(IM2)));
    elseif inputDay(comp_i) == 2
        binaryMaskAlign(:,:,comp_i) = binaryMask(:,:,comp_i);
    elseif inputDay(comp_i) == 3
    binaryMaskAlign(:,:,comp_i) = imwarp(binaryMask(:,:,comp_i),imagesAligned.IM3_tform,'OutputView',imref2d(size(IM2)));    
    end
end
%find correlation coefficients and pvalues for timeseries
[Rtime,Ptime] = corrcoef(timeseries);
%find correlation coefficients for 2d spatial masks
for corr_i = 1:length(inputDay)
    for corr_ii = corr_i:length(inputDay)
    Rmask(corr_i,corr_ii) = corr2(binaryMaskAlign(:,:,corr_i),binaryMaskAlign(:,:,corr_ii));
    Rmask(corr_ii,corr_i) = corr2(binaryMaskAlign(:,:,corr_ii),binaryMaskAlign(:,:,corr_i));
    end  
end
%overlay the RGB 
if length(inputDay) < 3
    [aMask,bMask] = CaBMI_XMASS(binaryMaskAlign(:,:,1),binaryMaskAlign(:,:,2),binaryMaskAlign(:,:,2));
    titleMask = ['RGB overlap ' dayROI{1} '(r) ' dayROI{2} '(c)'];
elseif length(inputDay) >= 3
    [aMask,bMask] = CaBMI_XMASS(binaryMaskAlign(:,:,1),binaryMaskAlign(:,:,2),binaryMaskAlign(:,:,3));
    titleMask = ['RGB overlap ' dayROI{1} '(r) ' dayROI{2} '(g) ' dayROI{3} '(b)'];
end
%plot the similarity matrix of time series, mask, and the time series
plotCompareTimeSeries = figure('units','normalized','outerposition',[0 0 1 0.5]);
sgtitle([batId ' Comparing ROIs']);
subplot(1,7,1);
imagesc(Rtime,[0 1]); %plot similarity matrix of time series
set(gca,'XTick',[1:length(inputDay)],'XTickLabel',dayROI,'YTick',[1:length(inputDay)],'YTickLabel',dayROI);
xlabel('Day.ROI');
ylabel('Day.ROI');
title('Correlation of Time series');
colorbar('southoutside');

subplot(1,7,2:4); %plot time series of selected ROIS
for plot_i = 1:length(inputDay)
    plot(1:minDurDay,timeseries(:,plot_i)+plot_i*6);
    hold on;
end
dataTicks = [6:6:length(inputDay)*6];
set(gca,'YTick',dataTicks,'YTickLabel',dayROI,'xlim',[0 minDurDay]);
xt = get(gca,'xtick');
set(gca,'XTick',xt, 'xticklabel',xt/(results.Fs*60));
title(['ROI Time Series']);
%ylabel('ROIs (Day.ROI)');
xlabel('Time (m)');
hold off
subplot(1,7,5);
imagesc(Rmask,[0 1]); %plot similarity matrix of masks
set(gca,'XTick',[1:length(inputDay)],'XTickLabel',dayROI,'YTick',[1:length(inputDay)],'YTickLabel',dayROI);
xlabel('Day.ROI');
ylabel('Day.ROI');
title('Correlation of ROI masks');
colorbar('southoutside');
subplot(1,7,6:7);
imshow(aMask);
set(gca,'xtick',[],'ytick',[]);
title(titleMask);
%concatenate the day/rois into 1 title so it can be used for saving title
dayROITitle = []; 
for i = 1:length(dayROI); 
    dayROITitle = [dayROITitle '_' dayROI{i}]; 
end
if saveFlag == 1
   saveas(plotCompareTimeSeries,[saveDir filesep batId '_corrROIS_days' dayROITitle '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.tif']);
   savefig(plotCompareTimeSeries,[saveDir filesep batId '_corrROIS_days' dayROITitle '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);  
end
inputGo = input('Go?');
end
%compData1 = ;

%%
if saveFlag == 1
    if strcmp(batId,'Gal')
        saveas(plotDay123_unaligned,[saveDir filesep 'Gal_200311and20_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.tif']);
        savefig(plotDay123_unaligned,[saveDir filesep 'Gal_200311and20_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
        saveas(plotDay123_aligned,[saveDir filesep 'Gal_200311and20_overLap_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.tif']);
        savefig(plotDay123_aligned,[saveDir filesep 'Gal_200311and20_overLap_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
        saveas(plotDay123_diffs,[saveDir filesep 'Gal_200311and20_diffs_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.tif']);
        savefig(plotDay123_diffs,[saveDir filesep 'Gal_200311and20_diffs_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
        saveas(plotDay123_progression,[saveDir filesep 'Gal_200311and20_progression_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.tif']);
        savefig(plotDay123_progression,[saveDir filesep 'Gal_200311and20_progression_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
        saveas(plotDay123_centroids,[saveDir filesep 'Gal_200311and20_centroids_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.tif']);
        savefig(plotDay123_centroids,[saveDir filesep 'Gal_200311and20_centroids_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
        saveas(plotDay123_all,[saveDir filesep 'Gal_200311and20_allPlots_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.tif']);
        savefig(plotDay123_all,[saveDir filesep 'Gal_200311and20_allPlots_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
        saveas(plotOverlapTimeSeries,[saveDir filesep 'Gal_200311and20_timeSeriesOverlap_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.tif']);
        savefig(plotOverlapTimeSeries,[saveDir filesep 'Gal_200311and20_timeSeriesOverlap_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
        saveas(plotCnmfeSeries,[saveDir filesep 'Gal_200311and20_timeSeries_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.tif']);
        savefig(plotCnmfeSeries,[saveDir filesep 'Gal_200311and20_timeSeries_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
                if rigidFlag == 0
            saveas(figNonrigid,[saveDir filesep 'Gal_200311and20_nonrigidAlign_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.tif']);
            savefig(figNonrigid,[saveDir filesep 'Gal_200311and20_nonrigidAlign_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
        end
    elseif strcmp(batId,'Gen')
        saveas(plotDay123_unaligned,[saveDir filesep 'Gen_200319and24_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.tif']);
        savefig(plotDay123_unaligned,[saveDir filesep 'Gen_200319and24_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
        saveas(plotDay123_aligned,[saveDir filesep 'Gen_200319and24_overLap_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.tif']);
        savefig(plotDay123_aligned,[saveDir filesep 'Gen_200319and24_overLap_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
        saveas(plotDay123_diffs,[saveDir filesep 'Gen_200319and24_diffs_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.tif']);
        savefig(plotDay123_diffs,[saveDir filesep 'Gen_200319and24_diffs_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
        saveas(plotDay123_progression,[saveDir filesep 'Gen_200319and24_progression_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.tif']);
        savefig(plotDay123_progression,[saveDir filesep 'Gen_200319and24_progression_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
        saveas(plotDay123_all,[saveDir filesep 'Gen_200319and24_allPlots_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.tif']);
        savefig(plotDay123_all,[saveDir filesep 'Gen_200319and24_allPlots_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
        saveas(plotDay123_centroids,[saveDir filesep 'Gen_200319and24_centroids_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.tif']);
        savefig(plotDay123_centroids,[saveDir filesep 'Gen_200319and24_centroids_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
        saveas(plotOverlapTimeSeries,[saveDir filesep 'Gen_200319and24_overlapTimeSeries_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.tif']);
        savefig(plotOverlapTimeSeries,[saveDir filesep 'Gen_200319and24_overlapTimeSeries_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
        saveas(plotCnmfeSeries,[saveDir filesep 'Gen_200319and24_timeSeries_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.tif']);
        savefig(plotCnmfeSeries,[saveDir filesep 'Gen_200319and24_timeSeries_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
        elseif strcmp(batId,'z2')
        saveas(plotDay123_unaligned,[saveDir filesep 'Zo2_200701and12_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.tif']);
        savefig(plotDay123_unaligned,[saveDir filesep 'Zo2_200701and12_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
        saveas(plotDay123_aligned,[saveDir filesep 'Zo2_200701and12_overLap_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.tif']);
        savefig(plotDay123_aligned,[saveDir filesep 'Zo2_200701and12_overLap_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
        saveas(plotDay123_diffs,[saveDir filesep 'Zo2_200701and12_diffs_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.tif']);
        savefig(plotDay123_diffs,[saveDir filesep 'Zo2_200701and12_diffs_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
        saveas(plotDay123_progression,[saveDir filesep 'Zo2_200701and12_progression_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.tif']);
        savefig(plotDay123_progression,[saveDir filesep 'Zo2_200701and12_progression_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
        saveas(plotDay123_all,[saveDir filesep 'Zo2_200701and12_allPlots_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.tif']);
        savefig(plotDay123_all,[saveDir filesep 'Zo2_200701and12_allPlots_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
        saveas(plotDay123_centroids,[saveDir filesep 'Zo2_200701and12_centroids_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.tif']);
        savefig(plotDay123_centroids,[saveDir filesep 'Zo2_200701and12_centroids_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
        saveas(plotOverlapTimeSeries,[saveDir filesep 'Zo2_200701and12_overlapTimeSeries_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.tif']);
        savefig(plotOverlapTimeSeries,[saveDir filesep 'Zo2_200701and12_overlapTimeSeries_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
        saveas(plotCnmfeSeries,[saveDir filesep 'Zo2_200701and12_timeSeries_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.tif']);
        savefig(plotCnmfeSeries,[saveDir filesep 'Zo2_200701and12_timeSeries_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
        if rigidFlag == 0
            saveas(figNonrigid,[saveDir filesep 'Gen_200319and24_nonrigidAlign_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.tif']);
            savefig(figNonrigid,[saveDir filesep 'Gen_200319and24_nonrigidAlign_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
        end
    end
end


