function [NormRateMat, occupancyMap]=  ImBat_2dHeatMap(flights,spikes,val)
% Where do indirect neurons 'live' in 2D space?
% Load in the extracted direct neurons, and the indirect cells

% get the indirect cell you want, and find when it is active

filt = 1.3;

% pre-allocate maps with Nans
NormRateMat = NaN(50,50);
occupancyMap = NaN(50,50);


stdx = 1;

if isempty(spikes)
    spikes =1;
end

% ROI = roi_ave.S_dec(cell,:);
% % figure(); plot(ROI)
% 
% [pks, locs] = findpeaks(double(ROI),'MinPeakDistance',10,'MinPeakHeight',std(ROI)*2);
% 
% if isempty(pks)
%     disp('empty matrix detected...');
%     pks = 10;
%     locs = 10;
% end
% figure();
% hold on;
% plot(ROI);
% plot(locs,pks,'*');


% Reconstruct the cursor
locs = spikes;
pks = ones(length(locs),1);

TData2.cursorA = flights(1,:); %(zscore(roi_ave2.interp_dff(1,:)) + zscore(roi_ave2.interp_dff(2,:)))  - (zscore(roi_ave2.interp_dff(3,:)) + zscore(roi_ave2.interp_dff(4,:)));
TData2.cursorB = flights(2,:); %(zscore(roi_ave2.interp_dff(5,:)) + zscore(roi_ave2.interp_dff(6,:)))  - (zscore(roi_ave2.interp_dff(7,:)) + zscore(roi_ave2.interp_dff(8,:)));

% Fake cursor ( pick random cells as 'cursor')
% TData2.cursorA = (zscore(roi_ave.F_dff(1,:)) + zscore(roi_ave.F_dff(3,:)))  - (zscore(roi_ave.F_dff(3,:)) + zscore(roi_ave.F_dff(4,:)));
% TData2.cursorB = (zscore(roi_ave.F_dff(5,:)) + zscore(roi_ave.F_dff(6,:)))  - (zscore(roi_ave.F_dff(7,:)) + zscore(roi_ave.F_dff(8,:)));


counter = 1;

% index into the 2D point in the direct cells, save the index
Rperm = randi(size(TData2.cursorA,2),1,size(locs,2)); % random timestamps
for i = 1:size(locs,2)

% Bin the location
for ii = 1:val;
BinR(1,counter) = TData2.cursorA(:,locs(i));
BinR(2,counter) = -TData2.cursorB(:,locs(i));

% Randomize...
rBinR(1,counter) = TData2.cursorA(:,Rperm(i));
rBinR(2,counter) = -TData2.cursorB(:,Rperm(i));

counter = counter+1;
end
% % Get the angle
% u = TData2.cursorA(:,(locs(i)-1):(locs(i)+1));
% v = TData2.cursorB(:,(locs(i)-1):(locs(i)+1));
% ThetaInDegrees(:,i) = atan2d(norm(cross(u,v)),dot(u,v));
end



% col = jet(round(max(pks)*100));
% %%% Make Scatter Plot
% figure(1);
% plot(TData2.cursorA,TData2.cursorB);
% hold on;
% for i = 1:size(locs,2)
% plot(BinR(1,i),BinR(2,i),'*','color',col(round(pks(i)*100),:));
% hold on;
% end
% 
% plot(zeros(11,1),-5:5,'k');
% plot(-5:5,zeros(11,1),'k');


% calculate occupancy:
x1 = TData2.cursorA(:,:);
y1 = -TData2.cursorB(:,:);
%h=fspecial('gaussian',filt,filt);
x1 = [ -3000 x1 3000];
y1 = [-3000 y1 3000];
[valuesF, centersF] = hist3([x1(:) y1(:)],[50 50]);
valuesF = valuesF+1;
occ2use = valuesF;
valuesF2 = imgaussfilt(valuesF,filt);
%imagesc((valuesF2));



% Make heatplot

%figure(2);
%h=fspecial('gaussian',filt,filt);
try
x = BinR(1,:);
y = BinR(2,:);
catch
    disp(' no activity in this bin, skipping..');
    return
end

x = [ -3000 x 3000];
y = [-3000 y 3000];
[values, centers] = hist3([x(:) y(:)],[50 50]);
values(1,1) = 0;
values(50,50) = 0;
%values2=imfilter(values,h,'circular');
values2 = imgaussfilt(values,filt);
%imagesc((values2));
NormRateMat = (values2./valuesF2);
NormRateMat = imgaussfilt(NormRateMat,filt);

occupancyMap = occ2use;
occupancyMap(1,1) = 0;
occupancyMap(50,50) = 0;
occupancyMap = imgaussfilt(occupancyMap,filt);
occupancyMap(occupancyMap>1) = 2;
occupancyMap(occupancyMap<1.1) = NaN;

% imagesc(NormRateMat);

% %filter artifacts
% filt2 = 2;
% h2=fspecial('gaussian',filt2,filt2);
% values2=imfilter(values2,h2,'circular');


% indexa = -5:.1:5;

% imagesc(indexa,indexa,zeros(100,100));
% hold on;
% imagesc(centers{2}([1 end]),centers{1}([1 end]),values2);
% hold on
% 
% plot(zeros(11,1),-5:5,'k');
% plot(-5:5,zeros(11,1),'k');
% 
% plot(ones(11,1)*stdx,-5:5,'--c','LineWidth',1.2);
% plot(-5:5,ones(11,1)*-stdx,'--c','LineWidth',1.2);
% 

%%% Get the angle




%%% Make Difference Images with null distribution:

% % Image one
% figure(3);
% imagesc(indexa,indexa,zeros(100,100));
% hold on;
% imagesc(centers{2}([1 end]),centers{1}([1 end]),values2);
% hold on
% 
% ylim([-5 5]);
% xlim([-5 5]);
% colormap(gray);
% 
% hold off
% 
% IMa = getframe(); % Get zoomed portion that is visible.
% IM = double(squeeze(IMa.cdata(:,:,1)));
% 
% %%% Image 2
% 
% 
% % Make random Image
% x = rBinR(2,:);
% y = rBinR(1,:);
% 
% [values, centers] = hist3([x(:) y(:)],[Sz Sz]);
% values2=imfilter(values,h,'circular');
% 
% %filter artifacts
% values2=imfilter(values2,h2,'circular');
% 
% 
% 
% imagesc(indexa,indexa,zeros(100,100));
% hold on;
% imagesc(centers{2}([1 end]),centers{1}([1 end]),values2);
% hold on
% 
% ylim([-5 5]);
% xlim([-5 5]);
% colormap(gray);
% 
% hold off
% 
% IMa = getframe(); % Get zoomed portion that is visible.
% IM2 = double(squeeze(IMa.cdata(:,:,1)));
% 
% 
% X = CaBMI_3D_compare(mat2gray(IM),mat2gray(IM2));
% [h,p1] = kstest2(BinR(2,:),rBinR(2,:),'Alpha',0.01);
% [h,p2] = kstest2(BinR(1,:),rBinR(1,:),'Alpha',0.01);
% title(['p =   ', num2str(p1), '  randomized control  p =   ', num2str(p2)]);
% ylim([-5 5]);
% xlim([-5 5]);



% figure();
% subplot(1,2,1);
% hold on;
% h2 = histogram(rBinR(1,:),'FaceColor','m');
% h1 = histogram(BinR(1,:),'FaceColor','g');
% h1.Normalization = 'probability';
% h1.BinWidth = 0.5;
% h2.Normalization = 'probability';
% h2.BinWidth = 0.5;
% subplot(1,2,2);
% hold on;
% h2 = histogram(rBinR(2,:),'FaceColor','m');
% h1 = histogram(BinR(2,:),'FaceColor','g');
% h1.Normalization = 'probability';
% h1.BinWidth = 0.5;
% h2.Normalization = 'probability';
% h2.BinWidth = 0.5;

% Generage figure




% figure(5);
% histogram(ThetaInDegrees);
