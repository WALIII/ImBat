function ImBat_Analysis_071019(ROI_Data,day);
% All Pix dispatch

% Get the flight data:
 if exist('Y') ==0; % load in Y from local directory
     disp( 'Y matrix is being loaded from local directory...');
     load([ROI_Data{day}.date,'/',ROI_Data{day}.folder,'/Motion_corrected_Data_DS.mat'])
 end
 

close all
numClust = 3; 
for i = 1:numClust
% all pix images
[~, Ydata{i}, ClustFlight{i}] = ImBat_Analysis_070119(ROI_Data,'Y',Y,'clust',i,'day',day);
Ydata2{i} = squeeze(mean(Ydata{i},4));
Ydata3{i} = Ydata2{i}- mean(Ydata2{i},3);
end

close all
% Plotting
disp( 'plotting flights...');

exp2use = 2;
filt2use = 3;

figure();
for i = 1:numClust
[im1_rgb norm_max_proj{i},I{i},idx_img{i}] = CABMI_allpxs(imresize(Ydata3{i}(:,:,150:250),3),'filt_rad',filt2use,'exp',exp2use);
end


[RGB1 RGB2] = CaBMI_XMASS(norm_max_proj{1},norm_max_proj{2},norm_max_proj{3},'HL',[0.00 0.95]);



figure();
for i = 1:numClust
subplot(1,numClust,i);
imshow(I{i});
end

figure();
col = {'r','g','b'};
subplot(1,2,1);
hold on;
for i = 1:numClust
    for ii = 1:size(ClustFlight{i},3)
    plot3(ClustFlight{i}(:,1,ii),ClustFlight{i}(:,2,ii),ClustFlight{i}(:,3,ii),'Color',col{i});
    end
end
grid on;
view(-49,24)
subplot(122);
imshow(squeeze(RGB1));

