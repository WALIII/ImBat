function        ImBat_inspect_noTuneCells(placeCells,results)
 batName = 'Ge';
 dateSesh = '200319';
Atemp = full(results.A); %spatial masks
% Plot binary mask of all neurons in the A matrix
        %convert A matrix into full matrix
        Ysiz = size(results.Cn);
        nRois = length(placeCells.cells_none);
        scaling = 10; %depends on size of frame and downsampling from extraction step
        col = zeros(nRois,3);

        
figCheckRoi = figure();

for roi_i = placeCells.cells_none
        subplot(4,1,1:3);
    imagesc(imresize(results.Cn,scaling));
    colormap(gray);
    hold on;
    
    %plot the map of the cell
    sgtitle([batName ' ' dateSesh ' ROI #' num2str(roi_i)]);
    roiMask = imresize(mat2gray(reshape(Atemp(:,roi_i),Ysiz(1),Ysiz(2))),10);
    binaryMask = imbinarize(roiMask);%mat2gray(reshape(resultsA{inputDay(roi_i)}(:,inputROI(roi_i)),Ysiz(1),Ysiz(2))));
    
    [y,x]=find(binaryMask);
    %get ROI coordinates
    ROI_coords(roi_i,1) = {x};%*scaling};
    ROI_coords(roi_i,2) = {y};%*scaling};
    %calculate centroids
    centroid(roi_i,1) = mean(ROI_coords{roi_i,1});%*scaling;
    centroid(roi_i,2) = mean(ROI_coords{roi_i,2});%*scaling;%get the centroid of mask
    %plot the mask
    p = plot(ROI_coords{roi_i,1},ROI_coords{roi_i,2},'LineWidth',4);
    title(['ROI Mask']);
    xticklabels([]);
    yticklabels([]);
   
    %plot the time series from each cell
    subplot(4,1,4);
        timeseries(:,roi_i) = zscore(smoothdata(results.C_raw(roi_i,1:end),2,'movmedian',30))';
    plot(timeseries(:,roi_i));%+plot_i*6);
    xt = get(gca,'xtick');
    set(gca,'XTick',xt, 'xticklabel',xt/(results.Fs*60));
    title(['ROI Time Series']);
    %ylabel('ROIs (Day.ROI)');
    xlabel('Time (m)');
    hold off
end