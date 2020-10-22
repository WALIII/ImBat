batId = 'Gen';
clustNum = 3;
saveFlag = 1;
saveTag = 'sMat';
%load first day placeCellStableROI data
if strcmp(batId,'Gal')
    %cd([dirTop(day1).folder filesep 'plots\200911-preDurPost cells across days']);
        load('200311to200320_Gal_activity_allTrials_allClusts_sMat_dff.mat'); %load activity for pre,dur,post
elseif strcmp(batId,'Gen')
        load('200319to200324_Gen_activity_allTrials_allClusts_sMat_dff.mat'); %load activity for pre,dur,post
end
%make saving directory
if saveFlag == 1
    saveDir1 = '\\169.229.54.11\server_home\users\tobias\flight\data_processed\topQualityData\analysis_done\plots\';
    if ~exist([saveDir1 datestr(now,'yymmdd') filesep 'maxProjFlightAlign'])
        mkdir([saveDir1 datestr(now,'yymmdd') filesep 'maxProjFlightAlign']);
    else
        disp('You have been working today...');
    end
    saveDir = [saveDir1 datestr(now,'yymmdd') filesep 'maxProjFlightAlign' filesep];
end

nDays = size(activity_allTrials.maxMeanFrames_dur{clustNum},2);
maxProjflightAlign = figure('units','normalized','outerposition',[0 0 0.8 1]); 
sgtitle([batId ': Flight Aligned Max Projections Pre,Dur,Post']);
if strcmp(batId,'Gal')
    ha = tight_subplot(3,6,[.02 .01],[.01 .1],[.01 .01]);
elseif strcmp(batId,'Gen') 
    ha = tight_subplot(2,6,[.02 .01],[.01 .1],[.01 .01]);
end
for day_i = 1:nDays
%     if strcmp(batId,'Gal')
%     subplot(3,3,day_i); 
%     elseif strcmp(batId,'Gen')
%         subplot(2,3,day_i)
%     end
%plot the flight aligned max projections    
    axes(ha(2*day_i-1));
    imagesc(activity_allTrials.maxMeanFrames_dur{2}{day_i});  
    hold on;
    colormap(gray); 
    title(['Day ' num2str(day_i) ' Clust ' num2str(clustNum)]);
    set(gca,'xticklabel',[],'yticklabel',[]);
    
    %plot the full session max projections
    axes(ha(2*day_i));
    imagesc(activity_allTrials.YmaxFull{day_i});  
    hold on;
    colormap(gray); 
    title(['Day ' num2str(day_i) ' full sesh']);
    set(gca,'xticklabel',[],'yticklabel',[]);
    
end

if saveFlag == 1
    if strcmp(batId,'Gal')
        saveas(maxProjflightAlign,[saveDir filesep 'Gal_200311and20_maxProjFlightAlign_' saveTag '_' datestr(now,'yymmdd-HHMM') '.tif']);
        savefig(maxProjflightAlign,[saveDir filesep 'Gal_200311and20_maxProjFlightAlign_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
    elseif strcmp(batId,'Gen')
        saveas(maxProjflightAlign,[saveDir filesep 'Gen_200319and24_maxProjFlightAlign_' saveTag '_' datestr(now,'yymmdd-HHMM') '.tif']);
        savefig(maxProjflightAlign,[saveDir filesep 'Gen_200319and24_maxProjFlightAlign_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
    end
end