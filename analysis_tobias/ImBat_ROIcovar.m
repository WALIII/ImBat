dir

batName = 'Gio';
dateSesh = '200407';
sessionType = 'fly-1';
distThresh = 10; %number of pixels to check if the cells are close enough to be considered same cell
corrThresh = 0.7; %max correlation of time series if cells are very close

distThresh = 10; %number of pixels to check if the cells are close enough to be considered same cell
corrThresh = 0.65; %max correlation of time series if cells are very close
scaling = 5; % 5x but depends on size of frame and downsampling from extraction step
minLim = 0.7; %limits for correlation imagesc
maxLim = 1; %limits for correlation imagesc
distThresh = distThresh * scaling;

g = dir('G*');
z = dir('Z*');
dirTop = vertcat(g,z); %find all folders in top quality directory

ROI_duplicate = cell(length(dirTop),1); %make cell for indices of duplicated ROIS
ROI_unique = cell(length(dirTop),1); %make cell for unique indices
distance = cell(length(dirTop),1); %make cell for recording pairwise distances of ROIS
timeCorr = cell(length(dirTop),1);%make cell for recording correlation between time series of nearby ROIs
ROI_coords = cell(length(dirTop),1);
centroid = cell(length(dirTop),1);


plot_ROI_refined = figure('units','normalized','outerposition',[0 0 0.9 0.9]);
for d = 1:length(dirTop)-2

    try %extract metadata names and enter processed folder
        cd([dirTop(d).name filesep 'extracted']);
        flyFolders = dir('*fly*extraction');
        batName = flyFolders(end).name(1:3);
        dateSesh = flyFolders(end).name(5:10);
        sessionType = flyFolders(end).name(12:16);
        
        cd(flyFolders(end).name);
        dirProcessed = dir('processed_*');
        if strcmp(batName(1),'G')
            cd(dirProcessed(1).name);
        else
            cd(dirProcessed(end).name);
        end
    catch
        cd(dirTop(d).name);
        flyFolders = dir('*fly*extraction');
        batName = flyFolders(end).name(1:3);
        dateSesh = flyFolders(end).name(5:10);
        sessionType = flyFolders(end).name(12:16);
        
        cd(flyFolders(end).name);
        dirProcessed = dir('processed_*');
        cd(dirProcessed(end).name);
    end
    
    load('results.mat');
    %convert A matrix into full matrix
    Atemp = full(results.A);
    Ysiz = size(results.Cn);
    %make array for coordinates of each cell
    ROI_coords{d} = cell(length(results.A(1,:)),2);
    centroid{d} = zeros(length(results.A(1,:)),2);
    
    % get ROI centroids for each cell;
    for i = 1:round(length(results.A(1,:)))
        %create 3d matrix with all ROI heat maps
        ROI2plot(:,:,i) = mat2gray(reshape(Atemp(:,i),Ysiz(1),Ysiz(2)));
        %binarize the coordinates into mask
        binaryImage = imbinarize(mat2gray(reshape(Atemp(:,i),Ysiz(1),Ysiz(2))));
        [y,x]=find(binaryImage);
        %get ROI coordinates
        ROI_coords{d}(i,1) = {x*scaling};
        ROI_coords{d}(i,2) = {y*scaling};
        %calculate centroids
        centroid{d}(i,1) = mean(ROI_coords{d}{i,1});%*scaling;
        centroid{d}(i,2) = mean(ROI_coords{d}{i,2});%*scaling;%get the centroid of mask
    end
    
    %compare euclidian distance between centroids and the timeseries correlation between
    %close ROIS
    distance{d} = zeros(round(length(results.A(1,:))));
    timeCorr{d} = zeros(round(length(results.A(1,:)))); 
    for ii = 1:round(length(results.A(1,:)))
        for iii = ii+1:round(length(results.A(1,:)))
            distance{d}(ii,iii) = sqrt((centroid{d}(iii,1)-centroid{d}(ii,1))^2 + (centroid{d}(iii,2)-centroid{d}(ii,2))^2);   %find distance between centroids
            if distance{d}(ii,iii) < distThresh %if closer than 10 pixels
                R = corrcoef(results.C_raw(ii,:),results.C_raw(iii,:)); %take correlation between two time series
                timeCorr{d}(ii,iii) = R(1,2);
                if timeCorr{d}(ii,iii) > corrThresh %if correlation is greater than 0.8
                    if snr(zscore(smoothdata(results.C_raw(ii,:),'movmedian',3)),results.Fs) > snr(zscore(smoothdata(results.C_raw(iii,:),'movmedian',3)),results.Fs) %unique index is the time series with higher snr
                        ROI_duplicate{d} = [ROI_duplicate{d} iii]; %duplicated index is low snr signal
                    else
                        ROI_duplicate{d} = [ROI_duplicate{d} ii];
                    end
                end
            end
        end
    end
    %identify unique indices that were not duplicated
    for u = 1:round(length(results.A(1,:)))
        if ismember(u,ROI_duplicate{d}) == 0
            ROI_unique{d} = [ROI_unique{d} u];
        end
    end
%%    
    %plot correlation image raw
    sgtitle([num2str(distThresh/scaling) ' pixels :: ' num2str(corrThresh) ' corr']);
    ax1 = subplot(length(results.C(:,1))+2,4,1:4:(length(results.C(:,1))*2));
    %subplot(length(results.C(:,1)),3,1:3:(length(results.C(:,1))*1.5)-2)
    imagesc(imresize(results.Cn,scaling),[minLim maxLim]); 
    colormap(ax1,gray);
    title([batName ' ' dateSesh ' ' sessionType ': ' num2str(length(ROI_unique{d})) '/' num2str(length(results.A(1,:))) ' ROI'])
    set(gca,'YDir','normal');
    axis 'off'
    drawnow;
    
    %plot correlation image with unique cells in blue and replicates in red
    remAdjust = rem(round(length(results.C(:,1))*2),8);
    if remAdjust == 0;
        adjust2 = 5;
        adjust4 = 6;
    elseif remAdjust == 6;
        adjust2 = 7;
        adjust4 = 8;
    elseif remAdjust == 4;
        adjust2 = 9;
        adjust4 = 10;
    elseif remAdjust == 2;
        adjust2 = 11;
        adjust4 = 12;
    end
    ax2 = subplot(length(results.C(:,1))+2,4,round(length(results.C(:,1))*2)+adjust2:4:length(results.C(:,1))*4);
    %subplot(length(results.C(:,1)),3,round(length(results.C(:,1))*1.5)+2:3:length(results.C(:,1))*3)
    imagesc(imresize(results.Cn,scaling),[minLim maxLim]); 
    colormap(ax2,gray);
    %title([batName ' ' dateSesh ' ' sessionType ': ' num2str(length(ROI_unique{d})) '/' num2str(length(results.A(1,:))) ' ROI'])
    set(gca,'YDir','normal');
    axis 'off' %'tight' 'equal'
    hold on
    for i = 1:length(ROI_unique{d}) %plot unique ROIs in blue
        try
            p = plot(ROI_coords{d}{ROI_unique{d}(i),1},ROI_coords{d}{ROI_unique{d}(i),2},'b','LineWidth',4);
            p.Color(4) = 0.1;
            t = text(centroid{d}(ROI_unique{d}(i),1),centroid{d}(ROI_unique{d}(i),2),num2str(ROI_unique{d}(i)));
            %t.Color(1:3) = col(i,:);
        catch
        end
    end
    for ii = 1:length(ROI_duplicate{d}) %plot duplicated ROIs in red
        try
            p = plot(ROI_coords{d}{ROI_duplicate{d}(ii),1},ROI_coords{d}{ROI_duplicate{d}(ii),2},'r','LineWidth',4);
            p.Color(4) = 0.1;
            t = text(centroid{d}(ROI_duplicate{d}(ii),1),centroid{d}(ROI_duplicate{d}(ii),2),num2str(ROI_duplicate{d}(ii)));
            %t.Color(1:3) = col(ii,:);        
        catch
        end
    end
    hold off
    drawnow;
    
    %plot distance covariance matrix
    distFull{d} = distance{d} + distance{d}.';
    ax3 = subplot(length(results.C(:,1))+2,4,2:4:round(length(results.C(:,1))*2)-8);
    %subplot(length(results.C(:,1)),3,2:3:round(length(results.C(:,1))*1.5)-6)
    imagesc(distFull{d});
    colormap(ax3,parula);
    colorbar;
    title('Distance between ROI pairs (um)');
    set(gca,'xtick',1:4:length(results.C(:,1)),'ytick',1:4:length(results.C(:,1)));
    ylabel('ROI #');
    xlabel([]);%('ROI #');
    drawnow; 
    
    %plot time correlation covariance matrix
    timeCorrFull{d} = timeCorr{d} + timeCorr{d}.';
    ax4 = subplot(length(results.C(:,1))+2,4,round(length(results.C(:,1))*2)+adjust4:4:length(results.C(:,1))*4);
    %subplot(length(results.C(:,1)),3,round(length(results.C(:,1))*1.5)+9:3:length(results.C(:,1))*3)
    imagesc(timeCorrFull{d});
    colormap(ax4,parula);
    colorbar;
    title('Correlation coefficient of nearby ROI pairs');
    set(gca,'xtick',1:4:length(results.C(:,1)),'ytick',1:4:length(results.C(:,1)));
    ylabel('ROI #');
    xlabel('ROI #');
    drawnow;
    
    %plot traces for each ROI
    subArea1 = 3:4:4*length(results.C(:,1));
    subArea2 = 4:4:4*length(results.C(:,1));
    subAreaComb = [subArea1,subArea2];
    ax5 = subplot(length(results.C(:,1))+2,4,subAreaComb');%3:3:length(results.C(:,1))*3)
    %subplot(length(results.C(:,1)),3,3:3:length(results.C(:,1))*3)
    hold on
    normSmoothData = zeros(length(results.C(:,1)),length(results.C(1,:)));
    for i = 1:length(results.C(:,1)) %[9 31 16 27 26 29 41 46 54 61 68 98 86 87 151]%
        normSmoothData(i,:) = zscore(smoothdata(results.C_raw(i,:),'movmedian',3)); %normalize each ROI trace
        try
            if sum(ismember(ROI_unique{d},i)) > 0
            plot(1:length(results.C_raw(1,:)),(zscore(smoothdata(results.C_raw(i,:),'movmedian',3)))+i*6,'b') %may have to tweak the +i*6 at the end
            else
                         plot(1:length(results.C_raw(1,:)),(zscore(smoothdata(results.C_raw(i,:),'movmedian',3)))+i*6,'r') %may have to tweak the +i*6 at the end   
            end
        catch
        end
    end 
    title(['ROI Activity: ' batName ' ' dateSesh ' ' sessionType])
    ylabel(['z-score dff: ' num2str(length(results.C(:,1))) ' ROIs'])
    %yticklabel(gca,
    xlabel('Time (min)')
    set(gca,'ytick',6:12:length(results.C(:,1))*6,'yticklabel',1:2:length(results.C(:,1)));
    set(gca,'xtick',1:(5*60*results.Fs):length(results.C(1,:)),'xticklabel',1:5:round(length(results.C(1,:))/(60*results.Fs)))
    xlim([1 length(results.C(1,:))])
    ylim([1 (length(results.C(:,1))+2)*results.Fs])
    axis 'tight'
    hold off
    drawnow;
    
    %save fig and tif of max projection
    %set(findall(maxFig,'-property','FontSize'),'FontSize',20);
    savefig(plot_ROI_refined,['/Users/periscope/Desktop/analysis/ROI_refined/plot_ROIrefined_' batName dateSesh sessionType '_' datestr(now,'yymmdd-hhMMss') '.fig']);
    saveas(plot_ROI_refined, ['/Users/periscope/Desktop/analysis/ROI_refined/plot_ROIrefined_' batName dateSesh sessionType '_' datestr(now,'yymmdd-hhMMss') '.tif']);
    saveas(plot_ROI_refined, ['/Users/periscope/Desktop/analysis/ROI_refined/plot_ROIrefined_' batName dateSesh sessionType '_' datestr(now,'yymmdd-hhMMss') '.svg']);
    clf;
    cd(dirTop(d).folder);
end

for u = 1:round(length(results.A(1,:)))
    if ismember(u,ROI_duplicate) == 0
        ROI_unique = [ROI_unique u];
    end
end

figure();
   imagesc(imresize(results.Cn,scaling)); colormap(gray);
   set(gca,'YDir','normal');
    hold on
    for f = 1:length(ROI_unique)
        try
            p = plot(ROI_coords{ROI_unique(f),1},ROI_coords{ROI_unique(f),2},'b','LineWidth',4);
            p.Color(4) = 0.2;
        catch
        end
    end
    for p = 1:length(ROI_duplicate)
        try
            p = plot(ROI_coords{ROI_duplicate(p),1},ROI_coords{ROI_duplicate(p),2},'r','LineWidth',4);
            p.Color(4) = 0.2;
        catch
        end
    end
    hold off 
    title([batName ' ' dateSesh ' ' sessionType ': d(' num2str(distThresh) '), c(' num2str(corrThresh) ')']);
ROI_refined.ROI_duplicate = ROI_duplicate;
ROI_refined.ROI_unique = ROI_unique;
ROI_refined.distance = distance;
ROI_refined.distFull = distFull;
ROI_refined.timeCorr = timeCorr;
ROI_refined.timeCorrFull = timeCorrFull;
ROI_refined.distThresh = distThresh/scaling; %number of pixels to check if the cells are close enough to be considered same cell
ROI_refined.corrThresh = corrThresh; %max correlation of time series if cells are very close
ROI_refined.scaling = scaling; % 5x but depends on size of frame and downsampling from extraction step
ROI_refined.minLim = minLim; %limits for correlation imagesc
ROI_refined.maxLim = maxLim; 
ROI_refined.ROI_coords = ROI_coords;
ROI_refined.centroid = centroid;


save(['/Users/periscope/Desktop/analysis/ROI_refined/ROI_refined_' num2str(distThresh/scaling) '_' num2str(corrThresh) '_' datestr(now,'yyMMdd-hhmmss') '.mat'],'ROI_refined');