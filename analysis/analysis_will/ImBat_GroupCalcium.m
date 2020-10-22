function [CombinedROI] = ImBat_GroupCalcium(ROI_Data,CellReg2);            % Calcium Imaging Data

cell_registered_struct = CellReg2

INDEX = cell_registered_struct.cell_to_index_map;

clear roi_combined temp_C_raw new_C_raw

for i = 1: size(ROI_Data,2); % for all days
    
    for ID = 1:size(cell_registered_struct.cell_to_index_map,1) % for all ROIs
        
        if i ==1 % if the first day
            if INDEX(ID,i) ==0; % if cell was not detected, pad with zeros
                new_C_raw(ID,:) = zeros(size(ROI_Data{i}.ROIs.results.C_raw(1,:)));
%                 new_C(ID,:) = zeros(size(ROI_Data{i}.ROIs.results.C(1,:)));
%                 new_S(ID,:) = zeros(size(ROI_Data{i}.ROIs.results.C(1,:)));
%                 new_timestamps(ID,:) = zeros(size(ROI_Data{i}.ROIs.results.C_raw(1,:)));

            else
                new_C_raw(ID,:) = ROI_Data{i}.ROIs.results.C_raw(INDEX(ID,i),:);
            end
            
        else % concat if not the first day
            
            if INDEX(ID,i) ==0; % if cell was not detected, pad with zeros
                temp_C_raw(ID,:) = zeros(size(ROI_Data{i}.ROIs.results.C_raw(1,:)));
            else
                temp_C_raw(ID,:) = ROI_Data{i}.ROIs.results.C_raw(INDEX(ID,i),:);
            end
        end
      end
      
        
        if  i ==1
            new_C_raw = new_C_raw;
        else
            
            new_C_raw = cat(2,new_C_raw,temp_C_raw);
            clear temp_C_raw
        end
        
end




%% A MATRIX of ROI Masks
clear ROI ROI_temp ROI_all RGBim2 ROI2plot RGBim
resize_factor = 4; 
% now, Plot a projection for all 'unique' ROIs on this interval:
for ID = 1:size(cell_registered_struct.cell_to_index_map,1) % for all ROIs
    for i = 1: size(ROI_Data,2); % for all days
        if INDEX(ID,i) ==0;
            ROI_temp(:,:,i) = imresize(zeros(size(squeeze(cell_registered_struct.spatial_footprints_corrected{1}(1,:,:)))),resize_factor);
        else
            ROI_temp(:,:,i) = imresize(squeeze(cell_registered_struct.spatial_footprints_corrected{i}(INDEX(ID,i),:,:)),resize_factor);
        end
    end
  %  figure(); for i = 1: 5; subplot(1,5,i); imagesc(squeeze(ROI_temp(:,:,i))); end

    ROI_all{ID} = ROI_temp;
    ROI(:,:,ID) = mat2gray(squeeze(mean(ROI_temp,3)));
end

figure(); imagesc(squeeze(sum(ROI,3)));


ROI2plot = ROI;
scaling  = 1.1; %1.1um per pixel

%% Color overlay Image:
roiHeatMax = max(ROI2plot,[],3); %plot all of the ROI heat maps
roiHeatMax = imresize(roiHeatMax, scaling);
for i = 1: size(ROI2plot,3)
    rand_col = randi(100)./100; % initialize rand seed color
    
    for ii = 1:3
        rand_col_prev = rand_col;
        rand_col = randi(100)./100;
        while abs(rand_col-rand_col_prev)<0.05 % get distinct colors
            rand_col = randi(100)./100;
        end
        
        RGBim(:,:,ii,i) = imresize(ROI2plot(:,:,i),scaling)*rand_col;
    end
end
RGBim2 = squeeze(max(RGBim,[],4));


figure()
imagesc(RGBim2)
title('All Unique ROIs');




% Plot to look for artifacts:

 
figure(); 
for i = 1:80;
hold on;
for ii=1:5;
subplot(2,5,ii)
 imagesc(squeeze(ROI_all{i}(:,:,ii)));
end
subplot(2,5,6:10)
plot(smooth(new_C_raw(i,:),30)); ylim([-4 20]); 
pause(); 
clf('reset');

end










% 
% % Group Unique ROI data together across days, based on 
% 
% for i = 1: size(ROI_Data,2)
% 
% if i ==1;
% A = 
% C = 
% C_r = 
% S = 
% 
% else
% 
% A = 
% C = 
% C_r = 
% S = 
% 
% 


