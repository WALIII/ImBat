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
for ii = 1:round(length(results.A(1,:)))
    for iii = ii+1:round(length(results.A(1,:)))
    if sqrt((centroid(iii,1)-centroid(ii,1))^2 + (centroid(iii,2)-centroid(ii,2)^2)) < 10 %if closer than 10 pixels
        if corrcoef(results.C_raw(ii),results.C_raw(iii)) > 0.8 %if correlation is greater than 0.8
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