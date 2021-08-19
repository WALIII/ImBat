

% compare 2D rate maps across days


% Generate the unique Flight Cluster:
[FlightAlignedROI_rando] = ImBat_Align_FC(CombinedROI,flightPaths,1);
FlightAlignedROI_rando1{1} = FlightAlignedROI_rando;


% Grab 2D rate maps for the top flights 
for i = 1:50;
[R{ii}(i,:),M2S{i}] = ImBat_Compare2D_ratemeaps_acrossDays(FlightAlignedROI,i,1);
[R{ii}(i,:),M2Sb{i}] = ImBat_Compare2D_ratemeaps_acrossDays(FlightAlignedROI,i,2);
[R{ii}(i,:),M2Sc{i}] = ImBat_Compare2D_ratemeaps_acrossDays(FlightAlignedROI,i,3);

[R2{ii}(i,:),M2S2{i}] = ImBat_Compare2D_ratemeaps_acrossDays(FlightAlignedROI_rando1,i,1);
close all
end


cell2use = 1
figure(); 
for i = 1:5;
subplot(3,5,i);
imagesc(M2S{cell2use}(:,:,i));
title(['Day ',num2str(i)]);
if i ==1;
    ylabel(['ROI ', num2str(cell2use), ' Top FLight path']);
end
axis equal
axis tight

subplot(3,5,i+5);
imagesc(M2Sb{cell2use}(:,:,i));
title(['Day ',num2str(i)]);
if i ==1;
    ylabel(['ROI ', num2str(cell2use), ' Second FLight path ']);
end
axis equal
axis tight

subplot(3,5,i+10);
imagesc(M2S2{cell2use}(:,:,i));
title(['Day ',num2str(i)]);
if i ==1;
    ylabel(['ROI ', num2str(cell2use), ' Unique FLights']);
end
axis equal
axis tight

end


% check for reward or boundry representation
for i = 1:50;
   Gu(:,:,i) = mat2gray(mean(M2S2{i}(:,:,:),3));
   Gc(:,:,i) = mat2gray(mean(M2S{i}(:,:,:),3));
   Gc2(:,:,i) = mat2gray(mean(M2Sb{i}(:,:,:),3));
   Gc3(:,:,i) = mat2gray(mean(M2Sc{i}(:,:,:),3));
end

meanGu = mean( Gu(:,:,:),3);
meanGc = mean( Gc(:,:,:),3);
meanGc2 = mean( Gc2(:,:,:),3);
meanGc3 = mean( Gc3(:,:,:),3);

figure(); 
imagesc(meanGc);
title('Place field largest clustered flight');

figure(); 
imagesc(meanGu);
title('Place field density for Unique flights');


figure(); 
imagesc((meanGc+meanGc2+meanGc3)/3);
title('Place field density for top3 flights');


figure(); 
subplot(1,2,1)
imagesc(meanGc2);
subplot(1,2,2)
imagesc(meanGc3);
title('Place field largest clustered flight');



% check similarity of fields across flight types:

figure(); 

cell2use = 2
figure(); 
for i = 1:5;
subplot(2,5,i);
imagesc(M2S{cell2use}(:,:,i));
title(['Day ',num2str(i)]);
if i ==1;
    ylabel(['ROI ', num2str(cell2use), ' Structured FLights']);
end


subplot(2,5,i+5);
imagesc(M2Sc{cell2use}(:,:,i));
title(['Day ',num2str(i)]);
if i ==1;
    ylabel(['ROI ', num2str(cell2use), ' Unique FLights']);
end
end


% take the difference

% compare
counter = 1;
for cell2use = 1:50
%     for i = 1:5
YYY(counter) = corr2(squeeze(mean(M2S{cell2use}(:,:,:),3)),squeeze(mean(M2Sb{cell2use}(:,:,:),3)));
YYY2(counter) = corr2(squeeze(mean(M2S{cell2use}(:,:,:),3)),squeeze(mean(M2S2{cell2use}(:,:,:),3)));

counter = counter+1;
%     end
end
figure();
hold on;
histogram(1-YYY,'BinWidth',0.1,'Normalization','probability')
histogram(1-YYY2,'BinWidth',0.1,'Normalization','probability')
title('2D corr of structured-structured (blue) vs structured-random (orange)')
xlabel('1-r')
ylabel('Relative frequency');
