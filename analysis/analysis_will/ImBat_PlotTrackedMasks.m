
function ImBat_PlotTrackedMasks(ROI_Data,CombinedROI,cell_registered_struct,days2use,cells2use);


% think it is better to use arrows:
% col = jet(6);
% counter = 1;
% ROI_all = CombinedROI.ROI_all;
% new_C_raw = CombinedROI.C_raw ;
% figure();
% for day2use = 1:5;
%     subplot(1,5,day2use);
%     hold on;
%
% % Plot Max projection
%     E = mat2gray(ROI_Data{day2use}.MaxProj_flights);
% % Plot PNR projection
%     E =  imresize(mat2gray(ROI_Data{day2use}.ROIs.results.Cn),4);
%
% imagesc(E); colormap(gray)
% axis tight ;
% counter = 1;
% for cell2use = [8 22 32 50] %[1 3 11 37 22];
% % get ROI location
% hold on
% I =  mat2gray(squeeze(ROI_all{cell2use}(:,:,day2use)));
% hold on;
% [~,a] = max(mean(I,2));
% [~,b] = max(mean(I,1));
% if a>1
% drawArrow = @(x,y) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0,'LineWidth',3,'MaxHeadSize',50,'color',col(counter,:));
%
% y1 = [a+40 a];
% x1 = [b-40 b];
%
% drawArrow(x1,y1);
% end
%     counter = counter+1;
%
%  clear a b I
% end
%
% end

%%


%ROI2plot = (:,:,zeros(length(Atemp(1,:)));
% get ROI centroids for top 30%;
for iiii = 1:2;
    % think it is better to use arrows:
    if cells2use>0;
        col = jet(size(cells2use,2));
    end
    counter = 1;
    ROI_all = CombinedROI.ROI_all;
    new_C_raw = CombinedROI.C_raw ;
    scaling = 4;
    
    figure();
    for i = 1: size(days2use,2) % = 1:5;
        day2use =  days2use(i);
        subplot(1,size(days2use,2),i);
        title(ROI_Data{day2use}.date);
        hold on;
        
        % Plot Max projection
        % E = mat2gray(ROI_Data{day2use}.MaxProj_flights);
        % Plot PNR projection
        E =  imresize(mat2gray(ROI_Data{day2use}.ROIs.results.Cn),4);
        
        imagesc(E,[0.25 1]); colormap(gray)
        axis tight ;
        counter = 1;
        
        Atemp = full(ROI_Data{day2use}.ROIs.results.A);
        Ysiz = size(ROI_Data{day2use}.ROIs.results.Cn);% only needed for Gal...
        GG = (cell_registered_struct.cell_to_index_map);
        for iii = 1:size(cells2use,2)%[8 22 32 50] %[1 3 11 37 22]; unique ID
            if cells2use>0;
                cell2use = GG(cells2use(iii),day2use);
                if cell2use>0;
                    %create 3d matrix with all ROI heat maps
                    % Atemp = ROI_all{cell2use}(:,:,day2use);
                    
                    %create 3d matrix with all ROI heat maps
                    %binarize the coordinates into mask
                    binaryImage = imbinarize(mat2gray(reshape(Atemp(:,cell2use),Ysiz(1),Ysiz(2))));
                    
                    [y,x]=find(binaryImage);
                    %get ROI coordinates
                    ROI_coords(:,1) = {x*scaling};
                    ROI_coords(:,2) = {y*scaling};
                    %calculate centroids
                    %     centroid(i,1) = mean(ROI_coords{i,1});%*scaling;
                    %     centroid(i,2) = mean(ROI_coords{i,2});%*scaling;%get the centroid of mask
                    
                    % imagesc(imresize(results.Cn,8)); colormap(gray);
                    hold on
                    for ii = 1:length(ROI_coords)
                        try
                            p = plot(ROI_coords{ii,1},ROI_coords{ii,2},'LineWidth',5,'color', col(iii,:));
                            if iiii ==1
                                p.Color(4) = 0;
                            else
                                p.Color(4) = 1;
                            end
                            
                        catch
                        end
                    end
                    clear  ROI_coords
                    hold off
                end
            end
        end
    end
end

end


%
% ROI_all = CombinedROI.ROI_all;
% new_C_raw = CombinedROI.C_raw ;
% figure();
% for cell2use = 1:20;
% for day2use = 1:5;
%     subplot(1,5,day2use);
% E = mat2gray(ROI_Data{day2use}.MaxProj_flights);
% I =  mat2gray(squeeze(ROI_all{cell2use}(:,:,day2use)));
% imshow(I, 'InitialMag', 500)
% imshow(E, 'InitialMag', 500)
% % Make a truecolor all-green image.
% if day2use ==1;
%     c1 = (randi(99,1,1)+1)/100;
%     c2 = (randi(99,1,1)+1)/100;
%     c3 = (randi(99,1,1)+1)/100;
% end
% green = cat(3, ones(size(E))*c1, ones(size(E))*c2, ones(size(E))*c3);
% hold on
% h = imshow(green);
% hold off
% % Use our influence map as the
% % AlphaData for the solid green image.
% set(h, 'AlphaData', I)
% end
% pause();
% clf('reset');
% end





%
% figure();
% for i = 1: 8;
% hold on;
% for ii=1:5;
% subplot(2,5,i)
% imagesc(squeeze(ROI_all{i}(:,:,ii)),[0 10]);
% end
% subplot(2,5,6:10)
% plot(smooth(new_C_raw(i,:),30)); ylim([-4 20]); pause(); end
% hold off;
% end
%
%


%
%
% % Load
%
% ROI_all = CombinedROI.ROI_all;
% new_C_raw = CombinedROI.C_raw ;
%
%     figure();
%     for i = 1:size(ROI_all,2);
%         hold on;
%         for ii=1:size(ROI_Data,2);
%             subplot(2,size(ROI_Data,2),ii)
%             imagesc(squeeze(ROI_all{i}(:,:,ii)));
%         end
%         subplot(2,size(ROI_Data,2),size(ROI_Data,2)+1:size(ROI_Data,2)*2)
%         plot(smooth(new_C_raw(i,:),30)); ylim([-4 20]);
%         pause();
%         clf('reset');
%
%     end
%
%
%     figure();
%    Blue_image = zeros(size(ROI_all{i},1), size(ROI_all{i},2), 3, 'uint8');
%    Blue_image(:, :, 3) = ones(size(ROI_all{i},1), size(ROI_all{i},2), 1, 'uint8')*256;
%    figure();
%    im = image(Blue_image)
%    figure();
%    h = imshow(Blue_image);
%    j = (squeeze(ROI_all{i}(:,:,ii)));
%     set(h, 'AlphaData', mat2gray(j))
%
%     figure();
% mx_proj = ROI_Data{1, 1}.MaxProj_flights;
% figure();
% hh = imshow(mx_proj); colormap(gray);
%     set(hh, 'AlphaData', h)
%
%
%     figure();
%     E = ROI_Data{1, 1}.MaxProj_flights;
%     green = cat(3, zeros(size(E)), ones(size(E)), zeros(size(E)));
% I =  mat2gray(squeeze(ROI_all{i}(:,:,ii)));
% imshow(E, 'InitialMag', 'fit')
% h = imshow(green);
% set(h, 'AlphaData', I)
%
%
%
% hold on
% h = imshow(green);
% hold off
%
%
