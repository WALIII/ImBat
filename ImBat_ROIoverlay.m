function ImBat_ROIoverlay(maxFig)

% Plot binary mask of all neurons in the A matrix
%convert A matrix into full matrix
Atemp = full(results.A);
% get ROI centroids;
for i = 1:length(results.A(1,:))
ROI2plot = mat2gray(reshape(Atemp(:,i),Ysiz(1),Ysiz(2)));
%binarize the coordinates into mask
binaryImage = imbinarize(mat2gray(reshape(Atemp(:,i),Ysiz(1),Ysiz(2))));
[y,x]=find(binaryImage);
%get ROI coordinates
ROI_coords(i,1) = {x};
ROI_coords(i,2) = {y};
%calculate centroids
centroid(i,1) = mean(ROI_coords{i,1});
centroid(i,2) = mean(ROI_coords{i,2});%get the centroid of mask
end

% Just PLot Centroids:
% figure(); 
% hold on;
% imagesc(ROI2plot);
plot(centroid(:,1),centroid(:,2),'o');


% % PLot the actual binary mask:
% figure(); 
% hold on;
% imagesc(ROI2plot);
figure();
for i = 1:length(ROI_coords)
p = plot(ROI_coords{i,1},ROI_coords{i,2},'LineWidth',5)
p.Color(4) = 0.2;
hold on
end