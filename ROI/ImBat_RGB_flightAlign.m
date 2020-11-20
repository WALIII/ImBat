function ImBat_RGB_flightAlign(batId,fullSeshTag,day1,day2,day3)
%batId = 'Gal';
%fullSeshTag = 1, chooses to look at whole session for 3 day comparison
clustNum = 2;
saveFlag = 1;
rigidFlag = 1;
saveTag = 'rigid';
dirAllTrials = pwd;
h = fspecial('disk',50);

%load first day placeCellStableROI data
if strcmp(batId,'Gal')
    %cd([dirTop(day1).folder filesep 'plots\200911-preDurPost cells across days']);
    load('Gal_200311to200324_activity_allTrials_allClusts_allTrials_sMat_newDff_newOrder.mat'); %load activity for pre,dur,post
elseif strcmp(batId,'Gen')
    load('Gen_200319to200324_activity_allTrials_allClusts_sMat_newDff_newOrder.mat'); %load activity for pre,dur,post
end
%make saving directory
if saveFlag == 1
    saveDir1 = '\\169.229.54.11\server_home\users\tobias\flight\data_processed\topQualityData\analysis_done\plots\';
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
        plotTitle = [saveTag ' ' batId ' clust ' num2str(clustNum) ' Full Sesh: Day ' num2str(day1) ' (r) v Day ' num2str(day2) ' (g) v Day ' num2str(day3) ' (b)'];
    end
end


%filter, and divide background to make the pixels all the same brightness
%resize and subtract background
IM1doub = imresize(double(IM1_raw),2);
IM2doub = imresize(double(IM2_raw),2);
IM3doub = imresize(double(IM3_raw),2);
bground1=imfilter(IM1doub,h,'replicate');%smoothdata(IM1doub,'gaussian',2);%
bground2=imfilter(IM2doub,h,'replicate');%smoothdata(IM2doub,'gaussian',2);%
bground3=imfilter(IM3doub,h,'replicate');%smoothdata(IM3doub,'gaussian',2);%
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
% Clip image ( binarize) above 25 percentile
image1(image1<0.25) =0;
image2(image2<0.25) =0;
image3(image3<0.25) =0;
% Clip image ( binarize) above 25 percentile
image1(image1>=0.75) =1;
image2(image2>=0.75) =1;
image3(image3>=0.75) =1;
% remove edges ( in case these are producing artifacts)
edgeCutoff = 60;
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
%RGB the images
[aOverlap,bOverlap] = CaBMI_XMASS(IM1_aligned,IM2_aligned,IM3_aligned);
%rgb the original unaligned for comparison
[aUnaligned,bUnaligned] = CaBMI_XMASS(IM1(1:minRow,1:minCol),IM2(1:minRow,1:minCol),IM3);
%filter, then average the images and subtract each from that average to show difference
hFilt = fspecial('disk',15);
IM1_aligned_filt = imfilter(IM1_aligned,hFilt,'replicate');
IM2_aligned_filt = imfilter(IM2_aligned,hFilt,'replicate');
IM3_aligned_filt = imfilter(IM3_aligned,hFilt,'replicate');
IM_sum = IM1_aligned_filt + IM2_aligned_filt + IM3_aligned_filt;
IM_mean = IM_sum/3;
IM1_diff = IM1_aligned_filt - IM_mean;
IM2_diff = IM2_aligned_filt - IM_mean;
IM3_diff = IM3_aligned_filt - IM_mean;
IM_1diff2 = IM2_aligned_filt - IM1_aligned_filt;
IM_2diff3 = IM3_aligned_filt - IM2_aligned_filt;
IM_1diff3 = IM3_aligned_filt - IM1_aligned_filt;

IM_combined(:,:,1) = IM1_aligned_filt;
IM_combined(:,:,2) = IM2_aligned_filt;
IM_combined(:,:,3) = IM3_aligned_filt;
IM_combined_max = max(IM_combined,[],3);
IM_combined_min = min(IM_combined,[],3);
IM1_maxDiff = IM1_aligned_filt - IM_combined_max;
IM2_maxDiff = IM2_aligned_filt - IM_combined_max;
IM3_maxDiff = IM3_aligned_filt - IM_combined_max;
IM1_minDiff = IM1_aligned_filt - IM_combined_min;
IM2_minDiff = IM2_aligned_filt - IM_combined_min;
IM3_minDiff = IM3_aligned_filt - IM_combined_min;

%% load masks and warp to match registered max projection for each day
if rigidFlag == 1 %only if have the rigid twarp output
    maskDir1 = '\\169.229.54.11\server_home\users\tobias\flight\data_processed\topQualityData\analysis_done\new3\';
    days = [day1 day2 day3];
    ROI2plot = cell(length(days),1);
    ROI_coords = cell(length(days),1);
    centroid = cell(length(days),1);
    col = cell(length(days));
    for day_i = 1:length(days)
        cd([maskDir1 batId]);
        dayDir = dir([batId(1:2) '*']);
        cd([dayDir(days(day_i)).name filesep 'extracted']);
        flyDir = dir('*fly-*extraction');
        cd(flyDir(end).name);
        processedDir = dir('processed*');
        load([processedDir(end).name filesep 'results.mat']);
        Atemp = full(results.A);
        resultsA{day_i} = Atemp;
        resultsCn{day_i} = results.Cn;
        
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

%%
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
image((aOverlap(:,:,:)))
title(plotTitle);
set(gca,'xticklabel',[],'yticklabel',[]);
%set(gca,'YDir','normal');

%plot aligned RGB figure and each day with centroids on top
plotDay123_centroids = figure('units','normalized','outerposition',[0 0 1 0.6]);
sgtitle([plotTitle ' and daily centroids']);
ha = tight_subplot(1,4,[.02 .01],[.01 .08],[.01 .01]);
set(0,'CurrentFigure',plotDay123_centroids);
axes(ha(1)); imshow((aOverlap(:,:,:)),[]); title('RGB Overlap');
axes(ha(2)); imshow(IM1_aligned,[]);
colormap(ha(2),gray);
hold on;
for roi_i = 1:size(ROI2plot_aligned{1},3)
        plot(ROI_coords{1}{roi_i,1},ROI_coords{1}{roi_i,2},'LineWidth',4,'Color',[col{1}(roi_i,:) 0.01]);
        p = text(centroid{1}(roi_i,1),centroid{1}(roi_i,2),num2str(roi_i));
        p.Color(1:3) = col{1}(size(ROI2plot_aligned{1},3) + 1 - roi_i,:);
end
title('Day 1 ROIs');
hold off
axes(ha(3)); imshow(IM2_aligned,[]);
colormap(ha(3),gray);
hold on;
for roi_i = 1:size(ROI2plot_aligned{2},3)
    plot(ROI_coords{2}{roi_i,1},ROI_coords{2}{roi_i,2},'LineWidth',4,'Color',[col{2}(roi_i,:) 0.01]);
    p = text(centroid{2}(roi_i,1),centroid{2}(roi_i,2),num2str(roi_i));
    p.Color(1:3) = col{2}(size(ROI2plot_aligned{2},3) + 1 - roi_i,:);
end
title('Day 2 ROIs');
hold off
axes(ha(4)); imshow(IM3_aligned,[]);
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
axes(ha(3)); imshow(IM1_diff,[-0.6 0.6]); title('Diff Day 1');
colormap(ha(3),fireice);
colorbar('southoutside');
axes(ha(4)); imshow(IM2_diff,[-0.6 0.6]); title('Diff Day 2');
colormap(ha(4),fireice);
colorbar('southoutside');
axes(ha(5)); imshow(IM3_diff,[-0.6 0.6]); title('Diff Day 3');
colormap(ha(5),fireice);
colorbar('southoutside');

%plot rgb and each progression of differences from 1 day to next
plotDay123_progression = figure('units','normalized','outerposition',[0 0 1 0.6]);
sgtitle([plotTitle ' and day to day progression']);
ha = tight_subplot(1,5,[.02 .01],[.01 .08],[.01 .01]);
set(0,'CurrentFigure',plotDay123_progression);
axes(ha(1)); imshow(aOverlap(:,:,:),[]); title('RGB Overlap');
%colorbar('southoutside');
axes(ha(2)); imshow(IM1_aligned,[]); title('Day 1');
colormap(ha(2),gray);
colorbar('southoutside');
axes(ha(3)); imshow(IM_1diff2,[-0.6 0.6]); title('Day 2 - Day 1');
colormap(ha(3),fireice);
colorbar('southoutside');
axes(ha(4)); imshow(IM_2diff3,[-0.6 0.6]); title('Day 3 - Day 2');
colormap(ha(4),fireice);
colorbar('southoutside');
axes(ha(5)); imshow(IM_1diff3,[-0.6 0.6]); title('Day 3 - Day 1');
colorbar('southoutside');
colormap(ha(5),fireice);

%% plot rgb and each progression of differences from 1 day to next
plotDay123_all = figure('units','normalized','outerposition',[0 0 1 1]);
sgtitle([plotTitle ' comparing days']);
ha = tight_subplot(4,5,[.02 .01],[.01 .08],[.01 .01]);
set(0,'CurrentFigure',plotDay123_all);
axes(ha(1)); imshow(aOverlap(:,:,:),[]); title('RGB Overlap');
%colorbar('southoutside');
axes(ha(2)); imshow(IM1_aligned,[]); title('Day 1');
%colorbar('southoutside');
axes(ha(3)); imshow(IM_1diff2,[-0.6 0.6]); title('Day 2 - Day 1');
%colorbar('southoutside');
axes(ha(4)); imshow(IM_2diff3,[-0.6 0.6]); title('Day 3 - Day 2');
%colorbar('southoutside');
axes(ha(5)); imshow(IM_1diff3,[-0.6 0.6]); title('Day 3 - Day 1');
%colorbar('southoutside');
colormap(ha(3),fireice);
colormap(ha(4),fireice);
colormap(ha(5),fireice);
axes(ha(6)); imshow(IM1_aligned,[]);
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
axes(ha(8)); imshow(IM1_diff,[-0.6 0.6]); title('Day 1 - mean');
%colorbar('southoutside');
axes(ha(9)); imshow(IM2_diff,[-0.6 0.6]); title('Day 2 - mean');
%colorbar('southoutside');
axes(ha(10)); imshow(IM3_diff,[-0.6 0.6]); title('Day 3 - mean');
colormap(ha(8),fireice);
colormap(ha(9),fireice);
colormap(ha(10),fireice);%colorbar('southoutside');
axes(ha(11)); imshow(IM2_aligned,[]);
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
axes(ha(13)); imshow(IM1_maxDiff,[-0.6 0.6]); title('Day 1  - max');
%colorbar('southoutside');
axes(ha(14)); imshow(IM2_maxDiff,[-0.6 0.6]); title('Day 2  - max');
%colorbar('southoutside');
axes(ha(15)); imshow(IM3_maxDiff,[-0.6 0.6]); title('Day 3 - max');
colormap(ha(13),fireice);
colormap(ha(14),fireice);
colormap(ha(15),fireice);%colorbar('southoutside');
axes(ha(16)); imshow(IM3_aligned,[]);
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
axes(ha(18)); imshow(IM1_minDiff,[-0.6 0.6]); title('Day 1 - min');
%colorbar('southoutside');
axes(ha(19)); imshow(IM2_minDiff,[-0.6 0.6]); title('Day 2 - min');
%colorbar('southoutside');
axes(ha(20)); imshow(IM3_minDiff,[-0.6 0.6]); title('Day 3 - min');
colormap(ha(18),fireice);
colormap(ha(19),fireice);
colormap(ha(20),fireice);%colorbar('southoutside');
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
        if rigidFlag == 0
            saveas(figNonrigid,[saveDir filesep 'Gen_200319and24_nonrigidAlign_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.tif']);
            savefig(figNonrigid,[saveDir filesep 'Gen_200319and24_nonrigidAlign_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
        end
    end
end


%%

%
%     if fullSeshTag == 0
%         if isempty(day3) %comparing flight aligned max proj
%             if day1==day2 %if looking at the cluster flight aligned vs full session
%
%                 plotTitle = [batId ' clust ' num2str(clustNum) ': Day ' num2str(day1) ' (r) v full session (c)'];
%             else %if looking at 1 day vs another day only
%                 IM1_raw = activity_allTrials.maxMeanFrames_dur{clustNum}{day1};
%                 IM2_raw = activity_allTrials.maxMeanFrames_dur{clustNum}{day2};
%                 %normalize the 2 images to each other
%                 IMcat = [IM1_raw IM2_raw];
%                 IMcat_gray = mat2gray(IMcat);
%                 IM1half = IMcat_gray(:,1:size(IMcat_gray,2)/2);
%                 IM1preFilt = IM1half.*2;
%                 IM1 = smoothdata(IM1preFilt,'gaussian',2);
%                 IM2preFilt = IMcat_gray(:,(size(IMcat_gray,2)/2)+1:end);
%                 IM2 = smoothdata(IM2preFilt,'gaussian',2);
%                 plotTitle = [batId ' clust ' num2str(clustNum) ': Day ' num2str(day1) ' (r) v Day ' num2str(day2) ' (c)']
%             end
%             IM3 = [];
%             % Alignment script
%             rowTemp = [size(IM1,1),size(IM2,1),size(IM3,1)];
%             colTemp = [size(IM1,2),size(IM2,2),size(IM3,2)];
%             minRow = min(rowTemp(rowTemp>0));
%             minCol = min(colTemp(colTemp>0));
%             [IM1_aligned,IM2_aligned] = ImBat_imageAlign(IM1(1:minRow,1:minCol),IM2(1:minRow,1:minCol),IM3);
%             IM3_aligned = [];
%             % Will RGB script
%             [a1,b1] = CaBMI_XMASS(IM1(1:minRow,1:minCol),IM2(1:minRow,1:minCol),IM3);
%             [a2,b2] = CaBMI_XMASS(IM1_aligned,IM2_aligned,IM3_aligned);
%         else %if looking at comparing 3 days
%             IM1_raw = activity_allTrials.maxMeanFrames_dur{clustNum}{day1};%mat2gray(activity_allTrials.maxMeanFrames_dur{clustNum}{day1});
%             IM2_raw = activity_allTrials.maxMeanFrames_dur{clustNum}{day2};%mat2gray(activity_allTrials.maxMeanFrames_dur{clustNum}{day2});
%             IM3_raw = activity_allTrials.maxMeanFrames_dur{clustNum}{day3};%mat2gray(activity_allTrials.maxMeanFrames_dur{clustNum}{day3});
%             %double size of data so everything isn't washed out by removing background
%             IM1doub = imresize(double(IM1_raw),2);
%             IM2doub = imresize(double(IM2_raw),2);
%             IM3doub = imresize(double(IM3_raw),2);
%             %find background of the image with a fspecial filter from above
%             bground1=imfilter(IM1doub,h,'replicate');%smoothdata(IM1doub,'gaussian',2);%
%             bground2=imfilter(IM2doub,h,'replicate');%smoothdata(IM2doub,'gaussian',2);%
%             bground3=imfilter(IM3doub,h,'replicate');%smoothdata(IM3doub,'gaussian',2);%
%             %remove the background and amplify low intensity pixels so they are same power as high intensity pixels
%             IM1_filt=IM1doub./(bground1+1);
%             IM2_filt=IM2doub./(bground2+1);
%             IM3_filt=IM3doub./(bground3+1);
%             %convert to grayscale and cut back apart into correct segments
%             IMcat = [IM1_filt IM2_filt IM3_filt];
%             IMcat_gray = mat2gray(IMcat);
%             IM1 = IMcat_gray(:,1:size(IMcat_gray,2)/3);
%             IM2 = IMcat_gray(:,(size(IMcat_gray,2)/3)+1:2*size(IMcat_gray,2)/3);
%             IM3 = IMcat_gray(:,(2*size(IMcat_gray,2)/3)+1:end);
%             %         test=imresize(double(IM1preFilt),1);
%             %         x = [1 5 10 15 20 25];
%             %         % Plot
%             %         figure();
%             %         for i = 1: 6
%             %             subplot(2,3,i)
%             %             h=fspecial('disk',x(i));
%             %             bground=imfilter(test,h);
%             %             test=test./(bground+1);
%             %             imagesc(test);
%             %             title(['disk = ',num2str(x(i))]);
%             %         end
%             minRow = min([size(IM1,1),size(IM2,1),size(IM3,1)]);
%             minCol = min([size(IM1,2),size(IM2,2),size(IM3,2)]);
%             % Alignment script
%             [IM1_aligned,IM2_aligned,IM3_aligned] = ImBat_imageAlign(IM1(1:minRow,1:minCol),IM2(1:minRow,1:minCol),IM3(1:minRow,1:minCol));
%             % Will RGB script
%             [a1,b1] = CaBMI_XMASS(IM1(1:minRow,1:minCol),IM2(1:minRow,1:minCol),IM3(1:minRow,1:minCol));
%             [a2,b2] = CaBMI_XMASS(IM1_aligned,IM2_aligned,IM3_aligned);
%             plotTitle = [batId ' clust ' num2str(clustNum) ': Day ' num2str(day1) ' (r) v Day ' num2str(day2) ' (g) v Day ' num2str(day3) ' (b)'];
%         end
%     elseif fullSeshTag == 1
%         %if looking at comparing 3 days full sessions
%         %IM1 = mat2gray(activity_allTrials.YmaxFull{day1});
%         %IM2 = mat2gray(activity_allTrials.YmaxFull{day2});
%         %IM3 = mat2gray(activity_allTrials.YmaxFull{day3});
%         IM1_raw = activity_allTrials.YmaxFull{day1};
%         IM2_raw = activity_allTrials.YmaxFull{day2};
%         IM3_raw = activity_allTrials.YmaxFull{day3};
%         %IM1_scale = imresize(IM1_raw,1);
%         %IM2_scale = imresize(IM2_raw,1);
%         %IM3_scale = imresize(IM3_raw,1);
%         IM1doub = imresize(double(IM1_raw),2);
%         IM2doub = imresize(double(IM2_raw),2);
%         IM3doub = imresize(double(IM3_raw),2);
%         bground1=imfilter(IM1doub,h,'replicate');%smoothdata(IM1doub,'gaussian',2);%
%         bground2=imfilter(IM2doub,h,'replicate');%smoothdata(IM2doub,'gaussian',2);%
%         bground3=imfilter(IM3doub,h,'replicate');%smoothdata(IM3doub,'gaussian',2);%
%         IM1_filt=IM1doub./(bground1+5);
%         IM2_filt=IM2doub./(bground2+5);
%         IM3_filt=IM3doub./(bground3+5);
%         IMcat = [IM1_filt IM2_filt IM3_filt];    IMcat_gray = mat2gray(IMcat);
%         IM1 = IMcat_gray(:,1:size(IMcat_gray,2)/3);
%         IM2 = IMcat_gray(:,(size(IMcat_gray,2)/3)+1:2*size(IMcat_gray,2)/3);
%         IM3 = IMcat_gray(:,(2*size(IMcat_gray,2)/3)+1:end);
%         %IM1 = smoothdata(IM1preFilt,'gaussian',2);
%         %IM2 = smoothdata(IM2preFilt,'gaussian',2);
%         %IM3 = smoothdata(IM3preFilt,'gaussian',2);
%         minRow = min([size(IM1,1),size(IM2,1),size(IM3,1)]);
%         minCol = min([size(IM1,2),size(IM2,2),size(IM3,2)]);
%         image1 = IM1;
%         image2 = IM2;
%         image3 = IM3;
%         % Clip image ( binarize) above 25 percentile
%         image1(image1<0.25) =0;
%         image2(image2<0.25) =0;
%         image3(image3<0.25) =0;
%         % Clip image ( binarize) above 25 percentile
%         image1(image1>=0.75) =1;
%         image2(image2>=0.75) =1;
%         image3(image3>=0.75) =1;
%         % remove edges ( in case these are producing artifacts)
%         image1(1:10,1:10) = 0;
%         image2(1:10,1:10) = 0;
%         image3(1:10,1:10) = 0;
%         image1(end-10:end,end-10:end) = 0;
%         image2(end-10:end,end-10:end) = 0;
%         image3(end-10:end,end-10:end) = 0;
%
%         [figNonrigid,imagesAligned] = ImBat_imageAlign_nonrigid(image1,image2,image3,IM1,IM2,IM3);
%         % Alignment script
%         %[IM1_aligned,IM2_aligned,IM3_aligned] = ImBat_imageAlign(IM1(1:minRow,1:minCol),IM2(1:minRow,1:minCol),IM3(1:minRow,1:minCol));
%         % Will RGB script
%         [a1,b1] = CaBMI_XMASS(IM1(1:minRow,1:minCol),IM2(1:minRow,1:minCol),IM3(1:minRow,1:minCol));
%         [a2,b2] = CaBMI_XMASS(imagesAligned.IM1_aligned,imagesAligned.IM2imagesAligned.,IM3_aligned);
%         plotTitle = ['nonrigid alpha' num2str(imagesAligned.alpha) ' ' batId ' clust ' num2str(clustNum) ' Full Sesh: Day ' num2str(day1) ' (r) v Day ' num2str(day2) ' (g) v Day ' num2str(day3) ' (b)'];
%     end
%
%
%
%

