
function ImBat_2D_overdays(CombinedROI,flightPaths,FlightAlignedROI);
% compare 2D rate maps across days


% Generate the combined Flight Cluster:
[FlightAlignedROI_new2] = ImBat_Combine_Clusters(CombinedROI,flightPaths,1:6);

close all

ii = 1;
maxCell = size(FlightAlignedROI{1}.S,1);%50;
% Grab 2D rate maps for the top flights 
for i = 1:maxCell;
    disp(['processing ROI ', num2str(i)]);
[R{ii}(i,:),M2S{i},Occ1{i}] = ImBat_Compare2D_ratemeaps_acrossDays(FlightAlignedROI,i,1); % primary cluster
[R1{ii}(i,:),M2S2{i},Occ2{i}] = ImBat_Compare2D_ratemeaps_acrossDays(FlightAlignedROI_new2,i,1); % all flights
%[R2{ii}(i,:),M2Sc{i},Occ3{i}] = ImBat_Compare2D_ratemeaps_acrossDays(FlightAlignedROI,i,3);
%[R3{ii}(i,:),M2S2{i},Occ4{i}] = ImBat_Compare2D_ratemeaps_acrossDays(FlightAlignedROI_rando1,i,1);
close all
end



roi2plot = 1;

figure(); 
for i = 1:5
    clear G oc2use imAlpha
subplot(2,5,i);
G = squeeze(M2S{roi2plot}(:,:,i));
G(G==0) = NaN;
% imagesc(G) 
oc2use = squeeze(Occ1{roi2plot}(:,:,i));
imAlpha=ones(size(oc2use));
imAlpha(isnan(oc2use))=0;
imagesc(G,'AlphaData',imAlpha);
set(gca,'color',0.9*[1 1 1]);
title(['Light ROI',num2str(roi2plot)])
colormap(jet);

axis equal
ylim([0 50])
xlim([0 50]);
end


% check for reward or boundry representation
for i = 1:maxCell;
   Gu(:,:,i) = mat2gray(sum(M2S{i}(:,:,:),3)./(sum(Occ1{i}(:,:,:),3)+.05));
   Gc(:,:,i) = mat2gray(sum(M2S2{i}(:,:,:),3)./(sum(Occ2{i}(:,:,:),3)+.05));
%    Gc2(:,:,i) = mat2gray(mean(M2Sb{i}(:,:,:),3));
%    Gc3(:,:,i) = mat2gray(mean(M2Sc{i}(:,:,:),3));
end

meanGu = nanmean( Gu(:,:,:),3);
meanGc = nanmean( Gc(:,:,:),3);
meanGu(meanGu ==1) = NaN;
meanGc(meanGc ==1) = NaN;


figure(); 
imagesc(meanGc);
title('Place field largest clustered flight');

figure(); 
imagesc(meanGu);
title('Place field density for Unique flights');


%%% look at behavioral consistancy 

%ImBat_ConsistMap(x1,y1,centersF);



%%% To DO: Stats, shuffling statistics. 

% shuff:
id=randi(size(meanGc,2),1,size(meanGc,1));
meanGu_2=cell2mat(arrayfun(@(x) circshift(meanGc(x,:),[1 id(x)]),(1:numel(id))','un',0));

figure();
% plot the hist
for i = 1: 45;
TrueSI = nansum(meanGu(:,i)/sum(abs(1-isnan(meanGu(:,i)))));
figure();
hold on;
histogram(nansum(meanGu_2,1)./sum(abs(1-isnan(meanGu(:,:)))));
hold on;
plot([TrueSI TrueSI],[0 30],'LineWidth',10);
title('Shuffled SI vs true SI');
legend('shuffeled','true');
pause();
clf
end


% figure(); 
% imagesc((meanGc+meanGc2+meanGc3)/3);
% title('Place field density for top3 flights');


% figure(); 
% subplot(1,2,1)
% imagesc(meanGc2);
% subplot(1,2,2)
% imagesc(meanGc3);
% title('Place field largest clustered flight');



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
imagesc(M2S2{cell2use}(:,:,i));
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
