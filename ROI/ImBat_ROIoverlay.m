function [ROIoverlay,goodCellIndex,badCellIndex] = ImBat_ROIoverlay(results,Ysiz,centroidFlag,binaryMaskFlag,roiHeatFlag,varargin)

global topROI

%manual inputs
%centroidFlag = 1;
%binaryMaskFlag = 1;
%roiHeatFlag = 1;
roiHeatFlagIndiv = 0;
scaling = 8; %depends on size of frame and downsizing from extraction step
topROILocal = topROI * 0.01; %look at first x% of ROIs

batName = [];
dateSesh = [];
sessionType = [];

% User inputs overrides
nparams=length(varargin);
for i=1:2:nparams
    switch lower(varargin{i})
        case 'batName'
            batName=varargin{i+1};
        case 'dateSesh'
            dateSesh = varargin{i+1};
        case 'sessionType'
            sessionType = varargin{i+1};
        case 'roiheatflagindiv'
            roiHeatFlagIndiv = varargin{i+1};
    end
end

% Plot binary mask of all neurons in the A matrix
%convert A matrix into full matrix
Atemp = full(results.A);
%ROI2plot = (:,:,zeros(length(Atemp(1,:)));
% get ROI centroids for top 30%;
for i = 1:round(length(results.A(1,:)))
    %create 3d matrix with all ROI heat maps
    ROI2plot(:,:,i) = mat2gray(reshape(Atemp(:,i),Ysiz(1),Ysiz(2)));
    %binarize the coordinates into mask
    binaryImage = imbinarize(mat2gray(reshape(Atemp(:,i),Ysiz(1),Ysiz(2))));
    [y,x]=find(binaryImage);
    %get ROI coordinates
    ROI_coords(i,1) = {x*scaling};
    ROI_coords(i,2) = {y*scaling};
    %calculate centroids
    centroid(i,1) = mean(ROI_coords{i,1});%*scaling;
    centroid(i,2) = mean(ROI_coords{i,2});%*scaling;%get the centroid of mask
end
%ROI_coords = smoothdata(ROI_coords,'gaussian',3); %filter the ROI coordinate mask so it is not so jagged
hold on
% modify labels for tick marks
xticks = get(gca,'xtick');
yticks = get(gca,'ytick');
scaling  = 1.1; %1.1um per pixel
newlabelsX = arrayfun(@(ax) sprintf('%g', scaling * ax), xticks, 'un', 0);
newlabelsY = arrayfun(@(ay) sprintf('%g', scaling * ay), yticks, 'un', 0);
set(gca,'xticklabel',newlabelsX,'yticklabel',newlabelsY);
title(['Max Projection: ' batName ' ' dateSesh ' ' sessionType]);
xlabel('um'); ylabel('um');
col = jet(length(ROI_coords));

%plot centroid over top of max projection
if centroidFlag == 1
    for i = 1:length(ROI_coords)
        try
            p = plot(centroid(i,1),centroid(i,2),'o');
            p.Color(1:3) = col(i,:);
            
        catch
        end
    end
    hold off
end
%plot binary mask over top of max projection
if binaryMaskFlag == 1
    hold on
    for i = 1:length(ROI_coords)
        try
            p = plot(ROI_coords{i,1},ROI_coords{i,2},'LineWidth',8);
            p.Color(4) = 0.4;
        catch
        end
    end
    hold off
end

%plot heat maps of all ROI
if roiHeatFlag == 1
    roiHeatMax = max(ROI2plot,[],3); %plot all of the ROI heat maps
    roiHeatMax = imresize(roiHeatMax, scaling);
    ROIoverlay = figure();
    imagesc(roiHeatMax);
    hold on
    %set(gca,'YDir','normal');
    xticks = get(gca,'xtick');
    yticks = get(gca,'ytick');
    scaling  = 1.1; %1.1um per pixel
    newlabelsX = arrayfun(@(ax) sprintf('%g', scaling * ax), xticks, 'un', 0);
    newlabelsY = arrayfun(@(ay) sprintf('%g', scaling * ay), yticks, 'un', 0);
    set(gca,'xticklabel',newlabelsX,'yticklabel',newlabelsY);
    title(['Max Projection ROI: ' batName ' ' dateSesh ' ' sessionType]);
    xlabel('um'); ylabel('um');
end

%plot heat maps of individual ROI
if roiHeatFlagIndiv == 1
    goodCellIndex = [];
    badCellIndex = [];
    for cell_i = 1:round(length(results.A(1,:)))
        %roiHeatMax = max(ROI2plot,[],3); %plot all of the ROI heat maps
        roiHeatMax = max(ROI2plot(:,:,cell_i),[],3); %plot individual heat maps for each ROI 1 at a time in for loop
        roiHeatMax = imresize(roiHeatMax, scaling);
        ROIoverlay = figure();
        imagesc(roiHeatMax);
        hold on
        %set(gca,'YDir','normal');
        xticks = get(gca,'xtick');
        yticks = get(gca,'ytick');
        scaling  = 1.1; %1.1um per pixel
        newlabelsX = arrayfun(@(ax) sprintf('%g', scaling * ax), xticks, 'un', 0);
        newlabelsY = arrayfun(@(ay) sprintf('%g', scaling * ay), yticks, 'un', 0);
        set(gca,'xticklabel',newlabelsX,'yticklabel',newlabelsY);
        title(['Max Projection ROI ' cell_i ': ' batName ' ' dateSesh ' ' sessionType]);
        xlabel('um'); ylabel('um');
        
        prompt = 'Is this a good cell? [1]=Y, [2]=N \n';
        str = input(prompt,'s');
        if strcmp(str,'1')
            goodCellIndex = [goodCellIndex cell_i];
        else
            badCellIndex = [badCellIndex cell_i];
        end
        %pause
        %close
    end
    save([pwd '/goodCellIdx.mat'],'goodCellIndex','badCellIndex');
end