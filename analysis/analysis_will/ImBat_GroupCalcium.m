function [CombinedROI,ROI_Data] = ImBat_GroupCalcium(ROI_Data,cell_registered_struct,aligned_data_struct);            % Calcium Imaging Data
% Combine all ROI data to a dsingle timeseries

% WAL3
% d10/21/2020



plotting = 0;

% make sure the data is sorted...
for i = 1:size(aligned_data_struct.file_names,2)
    G(i) = str2num(aligned_data_struct.file_names{i}(end-9:end-4));
end

[G1, ord] = sort(G);
% find data in ROI_Data:
counter = 1;
for i = 1:size(ROI_Data,2)
Gcomp(i) = str2num(ROI_Data{i}.date(3:end));
end

% match indexes:
for i = 1: size(G,2);
[idx,loc] = ismember(G,Gcomp);
end

for i = 1:size(loc,2)
ROI_Data2{i} = ROI_Data{loc(i)};
end
% ROI Data now needs to go back in extraction order...

% find data in ROI_Data:
counter = 1;
for i = 1:size(ROI_Data2,2)
G_backSort(i) = str2num(ROI_Data2{i}.date(3:end));
end

[G1, ord_R] = sort(G_backSort);
for i = 1:size(ord_R,2);
ROI_Data3{i} = ROI_Data2{ord_R(i)};
end

clear ROI_Data ROI_Data2;
ROI_Data = ROI_Data3;
clear ROI_Data2;
% [G1, ord] = sort(G);

for i = 1:size(cell_registered_struct.p_same_registered_pairs,1)
 CombinedROI.p_same_registered_pairs{i} = cell_registered_struct.p_same_registered_pairs{i}(ord,ord);
end

cell_registered_struct.cell_to_index_map = cell_registered_struct.cell_to_index_map(:,ord);
cell_registered_struct.spatial_footprints_corrected = cell_registered_struct.spatial_footprints_corrected(ord);
INDEX = cell_registered_struct.cell_to_index_map;

clear roi_combined temp_C_raw new_C_raw

for i = 1: size(ROI_Data,2); % for all days
    % Timestamp alignment first
    
    % ROI alignment next
    for ID = 1:size(cell_registered_struct.cell_to_index_map,1) % for all ROIs
        if i ==1 % if the first day
            if INDEX(ID,i) ==0; % if cell was not detected, pad with zeros
                new_C_raw(ID,:) = zeros(size(ROI_Data{i}.ROIs.results.C_raw(1,:)));
                new_C(ID,:) = zeros(size(ROI_Data{i}.ROIs.results.C(1,:)));
                new_S(ID,:) = zeros(size(ROI_Data{i}.ROIs.results.S(1,:)));
            else
                new_C_raw(ID,:) = ROI_Data{i}.ROIs.results.C_raw(INDEX(ID,i),:);
                new_C(ID,:) = ROI_Data{i}.ROIs.results.C(INDEX(ID,i),:);
                new_S(ID,:) = ROI_Data{i}.ROIs.results.S(INDEX(ID,i),:);
            end
        else % concat if not the first day
            if INDEX(ID,i) ==0; % if cell was not detected, pad with zeros
                temp_C_raw(ID,:) = zeros(size(ROI_Data{i}.ROIs.results.C_raw(1,:)));
                temp_C(ID,:) = zeros(size(ROI_Data{i}.ROIs.results.C(1,:)));
                temp_S(ID,:) = zeros(size(ROI_Data{i}.ROIs.results.S(1,:)));
                
            else
                temp_C_raw(ID,:) = ROI_Data{i}.ROIs.results.C_raw(INDEX(ID,i),:);
                temp_C(ID,:) = ROI_Data{i}.ROIs.results.C(INDEX(ID,i),:);
                temp_S(ID,:) = ROI_Data{i}.ROIs.results.S(INDEX(ID,i),:);
                
            end
        end
    end
    
    if  i ==1
        new_C_raw = new_C_raw;
        new_C = new_C;
        new_S = new_S;
        new_timestamps = ROI_Data{i}.Alignment.out.video_times(1:size(new_S,2)); % align timestamps
        day_vector = ones(size(new_S,2),1); % align timestamps

    else
        
        % Timestamp Alignmnet
        temp_new_timestamps = ROI_Data{i}.Alignment.out.video_times(1:size(temp_S,2)); % align timestamps
        new_timestamps = cat(1,new_timestamps,temp_new_timestamps+max(new_timestamps));
        temp_day_vector = ones(size(temp_S,2),1)+(i-1);
        day_vector = cat(1,day_vector,temp_day_vector);

        % ROI Alignemnt
        new_C_raw = cat(2,new_C_raw,temp_C_raw);
        new_C = cat(2,new_C,temp_C);
        new_S = cat(2,new_S,temp_S);
        clear temp_C_raw temp_new_timestamps temp_C temp_S
        
        
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

if plotting ==1;
    figure();
    for i = 1:size(ROI2plot,3);
        hold on;
        for ii=1:size(ROI_Data,2);
            subplot(2,size(ROI_Data,2),ii)
            imagesc(squeeze(ROI_all{i}(:,:,ii)));
        end
        subplot(2,size(ROI_Data,2),size(ROI_Data,2)+1:size(ROI_Data,2)*2)
        plot(smooth(new_C_raw(i,:),30)); ylim([-4 20]);
        pause();
        clf('reset');
        
    end
end



% Save
CombinedROI.C_raw =  new_C_raw;
CombinedROI.C = new_C;
CombinedROI.S = full(new_S);
CombinedROI.timestamps = new_timestamps;
CombinedROI.ROI_all = ROI_all;
CombinedROI.ROI = ROI;
CombinedROI.day_vector = day_vector;








%