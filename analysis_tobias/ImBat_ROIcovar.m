
distThresh = 10; %number of pixels to check if the cells are close enough to be considered same cell
corrThresh = 0.65; %max correlation of time series if cells are very close
scaling = 1; % 5x but depends on size of frame and downsampling from extraction step

g = dir('G*');
z = dir('Z*');
dirTop = vertcat(g,z); %find all folders in top quality directory

ROI_duplicate = cell(length(dirTop),1);
ROI_unique = cell(length(dirTop),1);
distance = cell(length(dirTop),1);%zeros(round(length(results.A(1,:))),1);
timeCorr = cell(length(dirTop),1);%zeros(round(length(results.A(1,:))),1);

for d = 1:length(dirTop)-2
    
    
    try
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
    % get ROI centroids for each cell;
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
    
    %compare euclidian distance between centroids and the timeseries correlation between
    %close ROIS
    distance{d} = zeros(round(length(results.A(1,:))));
    timeCorr{d} = zeros(round(length(results.A(1,:))));
    
    for ii = 1:round(length(results.A(1,:)))
        for iii = ii+1:round(length(results.A(1,:)))
            distance{d}(ii,iii) = sqrt((centroid(iii,1)-centroid(ii,1))^2 + (centroid(iii,2)-centroid(ii,2))^2);   %find distance between centroids
            if distance{d}(ii,iii) < distThresh %if closer than 10 pixels
                R = corrcoef(results.C_raw(ii,:),results.C_raw(iii,:)); %take correlation between two time series
                timeCorr{d}(ii,iii) = R(1,2);
                if timeCorr{d}(ii,iii) > corrThresh %if correlation is greater than 0.8
                    ROI_duplicate{d} = [ROI_duplicate{d} iii];
                end
            end
        end
    end
    
    for u = 1:round(length(results.A(1,:)))
        if ismember(u,ROI_duplicate{d}) == 0
            ROI_unique{d} = [ROI_unique{d} u];
        end
    end
    cd(dirTop(d).folder);
end

ROI_refined.ROI_duplicate = ROI_duplicate;
ROI_refined.ROI_unique = ROI_unique;
ROI_refined.distance = distance;
ROI_refined.timeCorr = timeCorr;

save(['/Users/periscope/Desktop/analysis/ROI_refined/ROI_refined_' datestr(now,'yyMMdd-hhmmss') '.mat'],'ROI_refined');