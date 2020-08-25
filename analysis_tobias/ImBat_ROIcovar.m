function [ROI_refined] = ImBat_ROIcovar
%this function checks distance between all pairs of ROIs and removes
%duplicate cells that are within the distance threshold and have a high
%correlation. It takes the cell that has the higher SNR between two
%duplicates. Plot the covariance matrices and the distribution of
%correlation coefficients.

saveFlag = 1; %do you want to save the figures and output structure?
saveDir1 = '\\169.229.54.11\server_home\users\tobias\flight\data_processed\topQualityData\analysis_done\plots\';
% Check if folder exists
if exist([saveDir1 datestr(now,'yymmdd')])>0;
    disp('Youve been working today..');
else
    mkdir([saveDir1 datestr(now,'yymmdd')])
end
saveDir = [saveDir1 datestr(now,'yymmdd') '\'];


distThresh = 10; %number of pixels to check if the cells are close enough to be considered same cell
corrThresh = 0.65; %max correlation of time series if cells are very close
scaling = 5; % 5x but depends on size of frame and downsampling from extraction step
minLim = 0.7; %limits for correlation imagesc
maxLim = 1; %limits for correlation imagesc
distThresh = distThresh * scaling;

g = dir('Ge*');
z = dir('Z1*');
dirTop = vertcat(g,z); %find all folders in top quality directory

ROI_duplicate = cell(length(dirTop),1); %make cell for indices of duplicated ROIS
ROI_unique = cell(length(dirTop),1); %make cell for unique indices
distance = cell(length(dirTop),1); %make cell for recording pairwise distances of ROIS
timeCorrClose = cell(length(dirTop),1);%make cell for recording correlation between time series of nearby ROIs
ROI_coords = cell(length(dirTop),1);
centroid = cell(length(dirTop),1);
timeCorrAll = cell(length(dirTop),1); %correlation coefficients between all ROIs pairwise
corrIndClose =  cell(length(dirTop),1); %make cell for indices of correlations between close cells
corrIndAll = cell(length(dirTop),1); %make cell for indices of correlations across all pairwise cells
timeCorrS = cell(length(dirTop),1); %make cell for indices of correlations across all pairwise cells with S matrix
timeCorrCloseS = cell(length(dirTop),1); %make cell for indices of correlations across all close cells with S matrix
corrIndCloseS =  cell(length(dirTop),1); %make cell for indices of correlations between close cells with S matrix
corrIndAllS = cell(length(dirTop),1); %make cell for indices of correlations across all pairwise cells with S matrix
spikeSmooth = cell(length(dirTop),1); %make cell for smoothed S matrix

for d = 1:length(dirTop)
    plot_ROI_refined = figure('units','normalized','outerposition',[0 0 0.9 0.9]);
    try %extract metadata names and enter processed folder
        cd([dirTop(d).name filesep 'extracted'])
        flyFolders = dir('*fly*extraction');
        batName{d} = flyFolders(end).name(1:3);
        dateSesh{d} = flyFolders(end).name(5:10);
        sessionType{d} = flyFolders(end).name(12:16);
        
        cd(flyFolders(end).name);
        dirProcessed = dir('processed_*');
        if strcmp(batName{d}(1),'G')
            cd(dirProcessed(end).name);
        else
            cd(dirProcessed(end).name);
        end
    catch
        cd(dirTop(d).name);
        flyFolders = dir('*fly*extraction');
        batName{d} = flyFolders(end).name(1:3);
        dateSesh{d} = flyFolders(end).name(5:10);
        sessionType{d} = flyFolders(end).name(12:16);
        
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
    timeCorrAll{d} = zeros(round(length(results.A(1,:))));
    timeCorrClose{d} = zeros(round(length(results.A(1,:)))); 
    for ii = 1:round(length(results.A(1,:)))
        for iii = ii+1:round(length(results.A(1,:)))
            distance{d}(ii,iii) = sqrt((centroid{d}(iii,1)-centroid{d}(ii,1))^2 + (centroid{d}(iii,2)-centroid{d}(ii,2))^2);   %find distance between centroids
            R1 = corrcoef(results.C_raw(ii,:),results.C_raw(iii,:)); %take correlation between two time series
            timeCorrAll{d}(ii,iii) = R1(1,2);
            
            w = gausswin(10); %smooth the S matrix by gaussian window of 10
            spikeSmooth{d}(ii,:) = filter(w,1,full(results.S(ii,:)));
            spikeSmooth{d}(iii,:) = filter(w,1,full(results.S(iii,:)));
            RS = corrcoef(spikeSmooth{d}(ii,:),spikeSmooth{d}(iii,:)); %take correlation between two smooth close series
            timeCorrS{d}(ii,iii) = RS(1,2);
            if distance{d}(ii,iii) < distThresh %if closer than 10 pixels
                %R = corrcoef(results.C_raw(ii,:),results.C_raw(iii,:)); %take correlation between two time series
                timeCorrClose{d}(ii,iii) = R1(1,2);
                %RSclose = corrcoef(spikeSmooth,spikeSmooth); %take correlation between two smooth close series
                timeCorrCloseS{d}(ii,iii) = RS(1,2);
                if timeCorrClose{d}(ii,iii) > corrThresh %if correlation is greater than 0.8
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
    title([batName{d} ' ' dateSesh{d} ' ' sessionType{d} ': ' num2str(length(ROI_unique{d})) '/' num2str(length(results.A(1,:))) ' ROI'])
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
    %title([batName{d} ' ' dateSesh{d} ' ' sessionType{d} ': ' num2str(length(ROI_unique{d})) '/' num2str(length(results.A(1,:))) ' ROI'])
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
    timeCorrFull{d} = timeCorrClose{d} + timeCorrClose{d}.';
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
            plot(1:length(results.C_raw(1,:)),normSmoothData(i,:)+i*6,'b') %may have to tweak the +i*6 at the end
            else
            plot(1:length(results.C_raw(1,:)),normSmoothData(i,:)+i*6,'r') %may have to tweak the +i*6 at the end   
            end
        catch
        end
    end 
    title(['ROI Activity: ' batName{d} ' ' dateSesh{d} ' ' sessionType{d}])
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
    
    %find all correlation indices across all pairs and close pairs of cells
    corrIndAll{d} = find(abs(timeCorrAll{d})>0);
    corrIndClose{d} = find(abs(timeCorrClose{d})>0); 
    corrIndAllS{d} = find(abs(timeCorrS{d})>0);
    corrIndCloseS{d} = find(abs(timeCorrCloseS{d})>0);
    %plot distributions of correlation coeffiecients for all and close ROIs
    plot_corrDist = figure();
    subplot(1,2,1);
    histogram(timeCorrAll{d}(corrIndAll{d}));
    hold on;
    title(['Dist R ALL ROIs: ' batName{d} ' ' dateSesh{d} ' ' sessionType{d}]);
    xlabel('Correlation coefficients');
    ylabel('# ROI pairs');
    subplot(1,2,2);
    histogram(timeCorrClose{d}(corrIndClose{d}),8,'FaceColor','r');
    hold on;
    title('Dist R CLOSE ROIs (cRaw matrix)');
    xlabel('Correlation coefficients');
    ylabel('# ROI pairs');
    hold off;
    
    %plot distributions of correlation coeffiecients for all and close ROIs
    plot_corrDistS = figure();
    subplot(1,2,1);
    histogram(timeCorrS{d}(corrIndAllS{d}));
    hold on;
    title(['Dist R ALL ROIs: ' batName{d} ' ' dateSesh{d} ' ' sessionType{d}]);
    xlabel('Correlation coefficients');
    ylabel('# ROI pairs');
    subplot(1,2,2);
    histogram(timeCorrCloseS{d}(corrIndCloseS{d}),8,'FaceColor','r');
    hold on;
    title('Dist R CLOSE ROIs (S matrix)');
    xlabel('Correlation coefficients');
    ylabel('# ROI pairs');
    hold off;
    
    if saveFlag == 1
    %save fig and tif of max projection
    %set(findall(maxFig,'-property','FontSize'),'FontSize',20);
    savefig(plot_ROI_refined,[saveDir 'plot_corrDist_' batName{d} dateSesh{d} sessionType{d} '_' datestr(now,'yymmdd-hhMMss') '.fig']);
    saveas(plot_ROI_refined, [saveDir 'plot_corrDist_' batName{d} dateSesh{d} sessionType{d} '_' datestr(now,'yymmdd-hhMMss') '.tif']);
    saveas(plot_ROI_refined, [saveDir 'plot_corrDist_' batName{d} dateSesh{d} sessionType{d} '_' datestr(now,'yymmdd-hhMMss') '.svg']);
    
    %save fig and tif of correlation distributions
    %set(findall(maxFig,'-property','FontSize'),'FontSize',20);
    savefig(plot_corrDist,[saveDir 'plot_corrDist_' batName{d} dateSesh{d} sessionType{d} '_' datestr(now,'yymmdd-hhMMss') '.fig']);
    saveas(plot_corrDist, [saveDir 'plot_corrDist_' batName{d} dateSesh{d} sessionType{d} '_' datestr(now,'yymmdd-hhMMss') '.tif']);
    saveas(plot_corrDist, [saveDir 'plot_corrDist_' batName{d} dateSesh{d} sessionType{d} '_' datestr(now,'yymmdd-hhMMss') '.svg']);
    %save fig and tif of correlation distributions
    %set(findall(maxFig,'-property','FontSize'),'FontSize',20);
    savefig(plot_corrDistS,[saveDir 'plot_corrDistS_' batName{d} dateSesh{d} sessionType{d} '_' datestr(now,'yymmdd-hhMMss') '.fig']);
    saveas(plot_corrDistS, [saveDir 'plot_corrDistS_' batName{d} dateSesh{d} sessionType{d} '_' datestr(now,'yymmdd-hhMMss') '.tif']);
    saveas(plot_corrDistS, [saveDir 'plot_corrDistS_' batName{d} dateSesh{d} sessionType{d} '_' datestr(now,'yymmdd-hhMMss') '.svg']);
    
    end
    close all;
    cd(dirTop(d).folder);
end

% for u = 1:round(length(results.A(1,:)))
%     if ismember(u,ROI_duplicate) == 0
%         ROI_unique = [ROI_unique u];
%     end
% end
% 
% figure();
%    imagesc(imresize(results.Cn,scaling)); colormap(gray);
%    set(gca,'YDir','normal');
%     hold on
%     for f = 1:length(ROI_unique)
%         try
%             p = plot(ROI_coords{ROI_unique(f),1},ROI_coords{ROI_unique(f),2},'b','LineWidth',4);
%             p.Color(4) = 0.2;
%         catch
%         end
%     end
%     for p = 1:length(ROI_duplicate)
%         try
%             p = plot(ROI_coords{ROI_duplicate(p),1},ROI_coords{ROI_duplicate(p),2},'r','LineWidth',4);
%             p.Color(4) = 0.2;
%         catch
%         end
%     end
%     hold off 
%     title([batName{d} ' ' dateSesh{d} ' ' sessionType{d} ': d(' num2str(distThresh) '), c(' num2str(corrThresh) ')']);
ROI_refined.batName = batName;
ROI_refined.dateSesh = dateSesh;
ROI_refined.sessionType = sessionType;
ROI_refined.ROI_duplicate = ROI_duplicate;
ROI_refined.ROI_unique = ROI_unique;
ROI_refined.distance = distance;
ROI_refined.distFull = distFull;
ROI_refined.timeCorrAll = timeCorrAll;
ROI_refined.timeCorrClose = timeCorrClose;
ROI_refined.timeCorrFull = timeCorrFull;
ROI_refined.timeCorrS = timeCorrS;
ROI_refined.timeCorrCloseS = timeCorrCloseS;
ROI_refined.corrIndAll = corrIndAll;
ROI_refined.corrIndClose = corrIndClose;
ROI_refined.corrIndAllS = corrIndAllS;
ROI_refined.corrIndCloseS = corrIndCloseS;
ROI_refined.distThresh = distThresh/scaling; %number of pixels to check if the cells are close enough to be considered same cell
ROI_refined.corrThresh = corrThresh; %max correlation of time series if cells are very close
ROI_refined.scaling = scaling; % 5x but depends on size of frame and downsampling from extraction step
ROI_refined.minLim = minLim; %limits for correlation imagesc
ROI_refined.maxLim = maxLim; 
ROI_refined.ROI_coords = ROI_coords;
ROI_refined.centroid = centroid;

%if saveFlag == 1
save([saveDir 'ROI_refined_c05' num2str(distThresh/scaling) '_' num2str(corrThresh) '_' datestr(now,'yyMMdd-hhmmss') '.mat'],'ROI_refined');
%save(['/Users/periscope/Desktop/analysis/ROI_refined/ROI_refined_' num2str(distThresh/scaling) '_' num2str(corrThresh) '_' datestr(now,'yyMMdd-hhmmss') '.mat'],'ROI_refined');
%end