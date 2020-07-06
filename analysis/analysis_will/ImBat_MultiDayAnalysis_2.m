function ImBat_MultiDayAnalysis_2(CellReg2,ROI_Data);


% Plot 'place cells' across days.


% 
topCells = size(CellReg2.cell_to_index_map,1);

  mkdir('PlaceCellsComp/fig'); % save as a figure file in local dir
   mkdir('PlaceCellsComp/jpg'); % Save as .jpg in local dir

offset = 0.0;
col = {'or','og','ob'};

for ii = 1:topCells;; % for each cell
    clear IDX
    for ix = 1: 3;
%         [a, b] =  find(CellReg2.cell_to_index_map(:,ix) == ii); 
%         IDX{ix} = a;
IDX{ix} = CellReg2.cell_to_index_map(ii,ix);
if IDX{ix} ==0;
   IM(:,:,ix) = zeros(ROI_Data{ix}.Ysiz_DS(1)*2,ROI_Data{ix}.Ysiz_DS(2)*2);
else
    IM(:,:,ix) = CellReg2.spatial_footprints_corrected{ix,:}(IDX{ix},:,:);
    %IM(:,:,ix) = reshape(full(ROI_Data{ix}.ROIs.results.A(:,IDX{ix})),ROI_Data{ix}.Ysiz_DS(1),ROI_Data{ix}.Ysiz_DS(2));
end
    end
    [RGB1 RGB2] = CaBMI_XMASS(IM(:,:,1),IM(:,:,2),IM(:,:,3));
figure(1); 
clf
subplot(1,2,1);
imagesc(squeeze(RGB1))
title('ROI Mask on days 1(r) 2(g) and 3(b)');
%     IDX
%     
subplot(1,2,2);
      hold on;
      A = ROI_Data{1, 1}.Alignment;
      B = ROI_Data{1, 2}.Alignment;
      C = ROI_Data{1, 3}.Alignment;
plot3(A.out.Location2(:,1),A.out.Location2(:,2),A.out.Location2(:,3),'Color',[0.7 0.1 0.1]);% plot the flight trajectory in space
plot3(B.out.Location2(:,1),B.out.Location2(:,2),B.out.Location2(:,3),'Color',[0.1 0.7 0.1]);% plot the flight trajectory in space
plot3(C.out.Location2(:,1),C.out.Location2(:,2),C.out.Location2(:,3),'Color',[0.1 0.1 0.7]);% plot the flight trajectory in space


for iii = 1:3
    clear Spike_times xy LX LY LZ
    if (IDX{iii}) ==0; ; 
    else
[~,xy] = find(ROI_Data{1, iii}.ROIs.results.S(IDX{iii},:)>0.1);  % get time neuron is active
Spike_times = ROI_Data{1, iii}.Alignment.out.video_times2(xy)-offset; % convert this to 'spike time'


try % this 'try/catch' is to avoid crashing if cells are not active in plotting window...
for i = 1:size(Spike_times,1)
try
    % Find the closest 'Location time' to the 'Spike time'
[minValue(:,i),closestIndex(:,i)] = min(abs(ROI_Data{1, iii}.Alignment.out.Location_time-Spike_times(i)));

LX(i) = ROI_Data{1, iii}.Alignment.out.flights(closestIndex(:,i),1);
LY(i) = ROI_Data{1, iii}.Alignment.out.flights(closestIndex(:,i),2);
LZ(i) = ROI_Data{1, iii}.Alignment.out.flights(closestIndex(:,i),3);
catch % we need this if Spiketime occurs before/after the location tracking was on..
continue
end
end

% display, in the title, how many bursts there were:
disp([num2str(size(LZ,2)-sum(isnan(LZ))),' Bursts in flight'])
scatter3(LX,LY,LZ,50,col{iii},'filled');
title(['Cell no ',num2str(ii),'  ',num2str(size(LX)),' Bursts in flight']);
xlim([-3000 3000]);
ylim([-3000 3000]);
zlim([-3000 3000]);

[Maps{iii}] = ImBat_2DheatMap(LX,LY,LZ);

catch % if cell was not active...
    disp('cell not active');
    LX = 0;
    LY = 0;
    LZ = 0;
 [Maps{iii}] = ImBat_2DheatMap(LX,LY,LZ);

end

    end


% Clear the buffer for the next cell:
clear LX LY LZ closestIndex Spike_times
end


    set(gcf, 'InvertHardcopy', 'off')

% Save 'place cells' as jpg and fig files..
saveas(gcf,['PlaceCellsComp/fig/','Cell_',num2str(ii)]);
saveas(gcf,['PlaceCellsComp/jpg/','Cell_',num2str(ii),'.jpg']);

% Plot place cell stability...
figure(2);
TitIdx = {'X-Y', 'X-Z', ' Y-Z'};
counter = 1;
for i = 1:3 % all days
    
    for ix = 1:3; % XY,YZ,XZ
    subplot(3,3,counter)
    try
    imagesc(Maps{i}{ix}.centers{2}([1 end]),Maps{i}{ix}.centers{1}([1 end]),Maps{i}{ix}.values2);
     Apix(:,:,ix,i) = Maps{i}{ix}.values2;
    catch
     Apix(:,:,ix,i) = zeros(41,41);
    end
    
    title(['day ',num2str(i),' ', TitIdx{ix},' Projection']);
    counter = counter+1;
    end
end


title2 = {'X-Y projection','X-Z projection','Y-Z projection'};
% figure(3); 
% for ixi = 1:3;
%         
%   [RGB1 RGB2] = CaBMI_XMASS(Apix(:,:,ixi,1),Apix(:,:,ixi,2),Apix(:,:,ixi,3));
%   subplot(1,3,ixi);
% image(squeeze(RGB1));
% title(title2(ixi));
% 
% end





clear Maps
 clf(figure(1))
clf(figure(2))
clf(figure(3))
%clf(figure(4)); 



  

end


  