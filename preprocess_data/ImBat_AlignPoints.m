function [tform,Markers2_aligned] = ImBat_AlignPoints(Markers1,Markers2);
% ImBat_AlignPoints.m

% Alignes markers2(moving), to markers1(fixed)



% 3D registration transfrom from daily room alignment 

M1_temp = squeeze(median(Markers1(:,:,:),1));
M2_temp = squeeze(median(Markers2(:,:,:),1));

outlierBound = 2000;
% remove outliers:
M1 = M1_temp( M1_temp(:,1)<-outlierBound | M1_temp(:,1)>=outlierBound | M1_temp(:,2)<-outlierBound | M1_temp(:,2)>=outlierBound,:);
M2 = M2_temp( M2_temp(:,1)<-outlierBound | M2_temp(:,1)>=outlierBound | M2_temp(:,2)<-outlierBound | M2_temp(:,2)>=outlierBound,:);


% make euclidean distance matrix:

% for i = 1:size(M1,1)
%     for ii = 1:size(M1,1)
%         GG = cat(1,M1(i,:),M2(ii,:));
%         out = squareform(pdist(GG));
%         dist2(i,ii) = out(2)
%     end
% end
% 
% figure(); imagesc(dist2);

% plot the offsets:
figure();
subplot(121);
grid on;
hold on;
for i = 1: size(M2,1);
    plot3(M2(i,1),M2(i,2),M2(i,3),'bo');
end
for i = 1: size(M1,1);
    plot3(M1(i,1),M1(i,2),M1(i,3),'ro');
    %     plot3(M4(i,1),M4(i,2),M4(i,3),'go');
end
title('UN Aligned Reference Points');
xlim([-4000 4000]);
ylim([-4000 4000]);


% turn XYZ matrix to pointcloud
moving = pointCloud(M2);
fixed = pointCloud(M1);
try
[tform,movingReg] = pcregistericp(moving,fixed)
catch
   disp('whoops'); 
end

ptCloudTformed = pctransform(moving,tform);


M2_t = ptCloudTformed.Location;

% plot the offsets:
subplot(122);
grid on;
hold on;
for i = 1:size(M2,1);
    plot3(M1(i,1),M1(i,2),M1(i,3),'ro');
    plot3(M2_t(i,1),M2_t(i,2),M2_t(i,3),'bo');
end
for i = 1:size(M1,1);
    plot3(M1(i,1),M1(i,2),M1(i,3),'ro');
end
title('Aligned Reference Points');
xlim([-4000 4000]);
ylim([-4000 4000]);

Markers2_aligned = M2_t; % these are the 'cleaned' and aligned.


