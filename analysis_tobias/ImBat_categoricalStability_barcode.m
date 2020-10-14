saveFlag = 1;
batId = 'Gal';
if strcmp(batId,'Gal')
    g = dir('Ga*');
elseif strcmp(batId,'Gen')
    g = dir('Ge*');
end
z = dir('Z1*');
dirTop = vertcat(g,z); %find all folders in top quality directory
nDays = length(dirTop);
%ROI_duplicate = cell(length(dirTop),1); %make cell for indices of duplicated ROIS

clusP = cell(length(dirTop),1);
clusPre = cell(length(dirTop),1);
clusPost = cell(length(dirTop),1);
roiP = cell(length(dirTop),1);
roiPre = cell(length(dirTop),1);
roiPost = cell(length(dirTop),1);
ppCells = cell(length(dirTop),1);
ppreCells = cell(length(dirTop),1);
ppostCells = cell(length(dirTop),1);


for day_i = 1:nDays
    %load results data
    dirExtracted = dir([dirTop(day_i).name filesep '*Extracted*']);
    load([dirTop(day_i).name filesep dirExtracted(end).name]);
    
    batName{day_i} = dirTop(day_i).name(1:3);
    dateSesh{day_i} = dirTop(day_i).name(5:10);
    
    if saveFlag == 1
    saveDir1 = '\\169.229.54.11\server_home\users\tobias\flight\data_processed\topQualityData\analysis_done\plots\';
    % Check if folder exists
    if exist([saveDir1 datestr(now,'yymmdd') filesep 'categoricalStability'])>0;
        disp('Youve been working today..');
    else
        mkdir([saveDir1 datestr(now,'yymmdd') filesep 'categoricalStability'])
    end
    saveDir = [saveDir1 datestr(now,'yymmdd') filesep 'categoricalStability' '\'];
    end
    
    clusP{day_i} = placeCellsAngStable.clusP;
    clusPre{day_i} = placeCellsAngStable.clusPre;
    clusPost{day_i} = placeCellsAngStable.clusPost;
    
    ppCells{day_i} = placeCellsAngStable.pp_cells;
    ppreCells{day_i} = placeCellsAngStable.ppre_cells;
    ppostCells{day_i} = placeCellsAngStable.ppost_cells;
    
    roiP{day_i} = placeCellsAngStable.roiP;
    roiPre{day_i} = placeCellsAngStable.roiPre;
    roiPost{day_i} = placeCellsAngStable.roiPost;
end

%% plot cluster 2 barcode of pre/dur/post cell classification
if strcmp(batId,'Gal')
    clustList = [2 2 2 2 2 2 2 2 2; 3 5 3 3 3 3 3 3 5; 1 3 4 6 4 11 1 1 8]; %cluster conversion so all match for Gal
elseif strcmp(batId,'Gen')
    clustList = [2 2 2 2 2; 3 3 3 3 3]; %cluster conversion so all match for Gen
end

for clust_i = 1:size(clustList,1);

categoricalStabilityBarcode = figure();
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
sgtitle([batId 'Pre (r)/Dur (g)/Post (b) Categorical Stability: Cluster ' num2str(clust_i+1)]);
nROI = size(ppCells{1},2);
ha = tight_subplot(nDays,nROI,[.01 .01],[.03 .1],[.02 .02]);
for day_ii = 1:nDays
     hold on;
      for roi_i = 1:nROI
        if ppreCells{day_ii}(clustList(clust_i,day_ii),roi_i) == 1
           axes(ha((nROI*day_ii)-nROI+roi_i));
           x = [0 1 1 0];
           y = [0 0 1 1];
           plot(x,y,'r','LineWidth',2);
           area(x,y,'FaceColor','red');
           hold on;
           y = [1 2 2 2];
            plot(x,y,'r','LineWidth',2);
           area(x,y,'FaceColor','red');
            axis([0 3 0 2]);
            set(gca,'xticklabel',[],'yticklabel',[]);
         end
        if ppCells{day_ii}(clustList(clust_i,day_ii),roi_i) == 1
           axes(ha((nROI*day_ii)-nROI+roi_i));
           x = [1 2 2 1];
           y = [0 0 1 1];
           plot(x,y,'g','LineWidth',2);
           hold on;
           y = [1 2 2 2];
            plot(x,y,'g','LineWidth',2);
           area(x,y,'FaceColor','green');
            axis([0 3 0 2]);
            set(gca,'xticklabel',[],'yticklabel',[]);
        end
        if ppostCells{day_ii}(clustList(clust_i,day_ii),roi_i) == 1
           axes(ha((nROI*day_ii)-nROI+roi_i));
           x = [2 3 3 2];
           y = [0 0 1 1];
           plot(x,y,'b','LineWidth',2);
           area(x,y,'FaceColor','blue');
           hold on;
           y = [1 2 2 2];
            plot(x,y,'b','LineWidth',2);
           area(x,y,'FaceColor','blue');
            axis([0 3 0 2]);
            set(gca,'xticklabel',[],'yticklabel',[]);
         end
        axes(ha((nROI*day_ii)-nROI+roi_i));
        if roi_i == 1
            ylabel(['Day ' num2str(day_ii)]);
        end
        if day_ii == 1
            title(['ROI ' num2str(roi_i)]);
        end
    end
    drawnow;
end

if saveFlag == 1
    if strcmp(batId,'Gal')
        saveas(categoricalStabilityBarcode,[saveDir filesep 'Gal_200311to20_categoricalStabilityBarcode_clust' num2str(clust_i+1) '_' datestr(now,'yymmdd-HHMMSS') '.tif']);
        savefig(categoricalStabilityBarcode,[saveDir filesep 'Gal_200311to20_categoricalStabilityBarcode_clust' num2str(clust_i+1) '_' datestr(now,'yymmdd-HHMMSS') '.fig']);
    elseif strcmp(batId,'Gen')
        saveas(categoricalStabilityBarcode,[saveDir filesep 'Gen_200319to24_categoricalStabilityBarcode_clust' num2str(clust_i+1) '_' datestr(now,'yymmdd-HHMMSS') '.tif']);
        savefig(categoricalStabilityBarcode,[saveDir filesep 'Gen_200319to24_categoricalStabilityBarcode_clust' num2str(clust_i+1) '_' datestr(now,'yymmdd-HHMMSS') '.fig']);
    end
end

end

% %% plot the functional stability of peak timings from beginning to end of recording session
% 
% day1 = 1;
% dayEnd = nDays;
% 
% day1Cells = cell(3,1); %3 cells for pre, dur, post
% day1Pks = cell(3,1);
% day1Locs = cell(3,1);
% day1Cells{1} = find(ppreCells{day1}(clust_i+1,:)); %preflight active cells
% day1Cells{2} = find(ppCells{day1}(clust_i+1,:)); %durflight active cells
% day1Cells{3} = find(ppostCells{day1}(clust_i+1,:)); %durflight active cells
% for p_i = 1:3
%     day1Pks{p_i} = zeros(length(day1Cells{p_i}),1);    
%     day1Locs{p_i} = zeros(length(day1Cells{p_i}),1);
% for roi_i = 1:length(day1Cells{p_i})
%     if p_i == 1
%     [day1Pks{p_i}(roi_i),day1Locs{p_i}(roi_i)] = findpeaks(placeCellsAngStable.ppre_cells_activity(roi_i,:));
%     elseif p_i == 2
%     [day1Pks{p_i}(roi_i),day1Locs{p_i}(roi_i)] = findpeaks(placeCellsAngStable.pp_cells_activity(roi_i,:));
%     elseif p_i == 3
%     [day1Pks{p_i}(roi_i),day1Locs{p_i}(roi_i)] = findpeaks(placeCellsAngStable.ppost_cells_activity(roi_i,:));
%     end
% end
% end
% 

% for clust_i = 1:size(clustList,1)
%     day1Pre{clust_i} = find(ppreCells{day1}(clust_i+1,:));
%     for 
%     day1Dur{clust_i} = find(ppCells{day1}(clust_i+1,:));
%     day1Post{clust_i} = find(ppostCells{day1}(clust_i+1,:));
%     
%     dayEndPre{clust_i} = find(ppreCells{dayEnd}(clust_i+1,:));
%     dayEndDur{clust_i} = find(ppCells{dayEnd}(clust_i+1,:));
%     dayEndPost{clust_i} = find(ppostCells{dayEnd}(clust_i+1,:));
%     
%     
% end










