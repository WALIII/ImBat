function ImBat_exportProjections(saveFlag)
centroidFlag = 1; %plot centroid ROI# on each cell
binaryMaskFlag = 1; %plot masks on top of max projection
roiHeatFlag = 1; %plot maks on correlation image with/without mask
%macDesk = '/Users/periscope/Desktop/analysis/flight/ROI_manual_selections/'; %local directory to save
%saveFlag = 1; %do you want to save the figures and output structure?
if saveFlag == 1
    saveDir1 = '\\169.229.54.11\server_home\users\tobias\flight\data_processed\topQualityData\analysis_done\plots\';
    % Check if folder exists
    if exist([saveDir1 datestr(now,'yymmdd') filesep 'manualROI_selections'])>0;
        disp('Youve been working today..');
    else
        mkdir([saveDir1 datestr(now,'yymmdd') filesep 'manualROI_selections'])
    end
    saveDir = [saveDir1 datestr(now,'yymmdd')  filesep 'manualROI_selections' '\'];
end

g = dir('G*'); %find all g bats
z = dir('Z*'); %find all z bats
dirTop = vertcat(g,z); %find all folders in top quality directory

%ROI_duplicate = cell(length(dirTop),1); %make cell for indices of duplicated ROIS

%plot_ROI_refined = figure('units','normalized','outerposition',[0 0 0.9 0.9]);
for day_i = length(dirTop)-2:length(dirTop)
    
    try %extract metadata names and enter processed folder
        cd([dirTop(day_i).name filesep 'extracted']);
        flyFolders = dir('*fly*extraction');
        if strcmp(flyFolders(end).name(1),'G')
            batName = flyFolders(end).name(1:3);
            dateSesh = flyFolders(end).name(5:10);
            sessionType = flyFolders(end).name(12:16);
        elseif strcmp(flyFolders(end).name(1:2),'Zo')
            batName = flyFolders(end).name(1:5);
            dateSesh = flyFolders(end).name(7:12);
            sessionType = flyFolders(end).name(14:18);
        else
            batName = flyFolders(end).name(1:4);
            dateSesh = flyFolders(end).name(6:11);
            sessionType = flyFolders(end).name(13:17);
        end
        
        
        cd(flyFolders(end).name);
        dirProcessed = dir('processed_*');
            cd(dirProcessed(end).name);
    catch
        cd(dirTop(day_i).name);
        flyFolders = dir('*fly*extraction');
        if strcmp(flyFolders(end).name(1:2),'Zo')
            batName = flyFolders(end).name(1:5);
            dateSesh = flyFolders(end).name(7:12);
            sessionType = flyFolders(end).name(14:18);
        else
            batName = flyFolders(end).name(1:4);
            dateSesh = flyFolders(end).name(6:11);
            sessionType = flyFolders(end).name(13:17);
        end
        
        cd(flyFolders(end).name);
        dirProcessed = dir('processed_*');
        cd(dirProcessed(end).name);
    end
    
    cellData = load('results.mat');
    videoData = load('Motion_corrected_Data_DS.mat');
    
    %make max projection from video
    [Ymax, Y, maxFig] = ImBat_Dff(videoData.Y,'batname',batName,'datesesh',dateSesh,'sessiontype',sessionType);
    hold on
    %save fig and tif of max projection
    %set(findall(maxFig,'-property','FontSize'),'FontSize',20);
    %savefig(maxFig,[imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/ROI/' fileName '_maxProject.fig']);
    %saveas(maxFig, [imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/ROI/' fileName '_maxProject.tif']);
    if saveFlag == 1
        savefig(maxFig,[saveDir batName '_' dateSesh '_' sessionType '_maxProject.fig']);
        saveas(maxFig,[saveDir batName '_' dateSesh '_' sessionType '_maxProject.tif']);
    end
    %make correlation image and max projection with ROI overlays
    [ROIoverlay,correlationImage,centroidMax] = ImBat_ROIoverlay(cellData.results,'centroid',centroidFlag,'binarymask',binaryMaskFlag,'roiheat',roiHeatFlag,'batname',batName,'datesesh',dateSesh,'sessiontype',sessionType);
    %save fig and tif of max projection
    %set(findall(maxFig,'-property','FontSize'),'FontSize',20);
    %savefig(correlationImage,[imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/ROI/' fileName '_correlationImage.fig']);
    %saveas(correlationImage, [imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/ROI/' fileName '_correlationImage.tif']);
    if saveFlag == 1
        savefig(correlationImage,[saveDir batName '_' dateSesh '_' sessionType '_correlationImage.fig']);
        saveas(correlationImage,[saveDir batName '_' dateSesh '_' sessionType '_correlationImage.tif']);
        savefig(ROIoverlay,[saveDir batName '_' dateSesh '_' sessionType '_ROIoverlay.fig']);
        saveas(ROIoverlay,[saveDir batName '_' dateSesh '_' sessionType '_ROIoverlay.tif']);
        savefig(centroidMax,[saveDir batName '_' dateSesh '_' sessionType '_ROIcentroid.fig']);
        saveas(centroidMax,[saveDir batName '_' dateSesh '_' sessionType '_ROIcentroid.tif']);
    end
    
    maxName = [saveDir batName '_' dateSesh '_' sessionType '_maxProject.tif'];
    corrName = [saveDir batName '_' dateSesh '_' sessionType '_correlationImage.tif'];
    combName = {maxName corrName};
    combName = string(combName);
    combFigure = figure('units','normalized','outerposition',[0 0 1 1]);
    montage(combName);
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    if saveFlag == 1
        savefig(combFigure,[saveDir batName '_' dateSesh '_' sessionType '_maxCorrProject.fig']);
        saveas(combFigure,[saveDir batName '_' dateSesh '_' sessionType '_maxCorrProject.tif']);
    end
    close all;
    cd(dirTop(day_i).folder);
    
end