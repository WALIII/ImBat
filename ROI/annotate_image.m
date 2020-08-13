function ROI = annotate_image(im)
an = Annotator(im);
an.wait();
masks = an.getMasks();


for i = 1: size(masks,3)

    [row,col] = find(masks(:,:,i)== 1);
    
    ROI.coordinates{1,i} = [col row];
	ROI.stats(i).Centroid=mean(ROI.coordinates{i});
	ROI.stats(i).Diameter=max(pdist(ROI.coordinates{i},'euclidean'));
	k=convhull(ROI.coordinates{i}(:,1),ROI.coordinates{i}(:,2));
	ROI.stats(i).ConvexHull=ROI.coordinates{i}(k,:);
    
end

ROI.type = 'Image';
im2 = imread(im);
ROI.reference_image = im2;
b = im(1:end-4);
b = strcat(b,'_ROI_DATA.mat');
save(b,'ROI')
end

