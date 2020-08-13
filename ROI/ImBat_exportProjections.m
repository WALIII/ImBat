centroidFlag = 1; %plot centroid ROI# on each cell
binaryMaskFlag = 1; %plot masks on top of max projection
roiHeatFlag = 1; %plot maks on correlation image with/without mask
macDesk = '/Users/periscope/Desktop/analysis/flight/ROI_manual_selections/'; %local directory to save


g = dir('G*'); %find all g bats
z = dir('Z*'); %find all z bats
dirTop = z;%vertcat(g,z); %find all folders in top quality directory

%ROI_duplicate = cell(length(dirTop),1); %make cell for indices of duplicated ROIS

%plot_ROI_refined = figure('units','normalized','outerposition',[0 0 0.9 0.9]);
for d = length(dirTop):length(dirTop)-2

    try %extract metadata names and enter processed folder
        cd([dirTop(d).name filesep 'extracted']);
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
        if strcmp(batName(1),'G')
            cd(dirProcessed(1).name);
        else
            cd(dirProcessed(end).name);
        end
    catch
        cd(dirTop(d).name);
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
    savefig(maxFig,[macDesk batName '_' dateSesh '_' sessionType '_maxProject.fig']);
    saveas(maxFig,[macDesk batName '_' dateSesh '_' sessionType '_maxProject.tif']);
    %make correlation image and max projection with ROI overlays
    [ROIoverlay,correlationImage,centroidMax] = ImBat_ROIoverlay(cellData.results,'centroid',centroidFlag,'binarymask',binaryMaskFlag,'roiheat',roiHeatFlag,'batname',batName,'datesesh',dateSesh,'sessiontype',sessionType);
    %save fig and tif of max projection
    %set(findall(maxFig,'-property','FontSize'),'FontSize',20);
    %savefig(correlationImage,[imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/ROI/' fileName '_correlationImage.fig']);
    %saveas(correlationImage, [imageFolders(kk).folder '/' imageFolders(kk).name '/' analysis_Folder '/ROI/' fileName '_correlationImage.tif']);
    savefig(correlationImage,[macDesk batName '_' dateSesh '_' sessionType '_correlationImage.fig']);
    saveas(correlationImage,[macDesk batName '_' dateSesh '_' sessionType '_correlationImage.tif']);    
    savefig(ROIoverlay,[macDesk batName '_' dateSesh '_' sessionType '_ROIoverlay.fig']);
    saveas(ROIoverlay,[macDesk batName '_' dateSesh '_' sessionType '_ROIoverlay.tif']);    
    savefig(centroidMax,[macDesk batName '_' dateSesh '_' sessionType '_ROIcentroid.fig']);
    saveas(centroidMax,[macDesk batName '_' dateSesh '_' sessionType '_ROIcentroid.tif']);    
    
    maxName = [macDesk batName '_' dateSesh '_' sessionType '_maxProject.tif'];
    corrName = [macDesk batName '_' dateSesh '_' sessionType '_correlationImage.tif'];
    combName = {maxName corrName};
    combName = string(combName);
    combFigure = figure('units','normalized','outerposition',[0 0 1 1]);
    montage(combName);
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    savefig(combFigure,[macDesk batName '_' dateSesh '_' sessionType '_maxCorrProject.fig']);
    saveas(combFigure,[macDesk batName '_' dateSesh '_' sessionType '_maxCorrProject.tif']);    
    close all;
    cd(dirTop(d).folder);

end