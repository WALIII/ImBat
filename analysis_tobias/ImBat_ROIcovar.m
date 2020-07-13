distThresh = 10; %number of pixels to check if the cells are close enough to be considered same cell
corrThresh = 0.8; %max correlation of time series if cells are very close
scaling = 1; % 5x but depends on size of frame and downsampling from extraction step
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
ROI_duplicate = [];
ROI_unique = [];
distance = zeros(round(length(results.A(1,:))),1);
timeCorr = zeros(round(length(results.A(1,:))),1);

for ii = 1:round(length(results.A(1,:)))
    for iii = ii+1:round(length(results.A(1,:)))
        distance(ii,iii) = sqrt((centroid(iii,1)-centroid(ii,1))^2 + (centroid(iii,2)-centroid(ii,2))^2);   %find distance between centroids
        if distance(ii,iii) < distThresh %if closer than 10 pixels
            R = corrcoef(results.C_raw(ii,:),results.C_raw(iii,:)); %take correlation between two time series
            timeCorr(ii,iii) = R(1,2);
            if timeCorr(ii,iii) > corrThresh %if correlation is greater than 0.8    
                ROI_duplicate = [ROI_duplicate iii];
            end
        end
    end
end

for u = 1:round(length(results.A(1,:)))
    if ismember(u,ROI_duplicate) == 0
        ROI_unique = [ROI_unique u];
    end
end