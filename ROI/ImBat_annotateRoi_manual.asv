saveFlag = 1; %do you want to save the figures and output structure?
reexportFlag = 0; %do you want to reproduce the dff and corr images in new folder
saveTag = 'allTrials_sMat_large_1to4';
if saveFlag == 1
saveDir1 = '\\169.229.54.11\server_home\users\tobias\flight\data_processed\topQualityData\analysis_done\plots\';
% Check if folder exists
if exist([saveDir1 datestr(now,'yymmdd') filesep 'manualROI_selections'])>0;
    disp('Youve been working today..');
else
    mkdir([saveDir1 datestr(now,'yymmdd') filesep 'manualROI_selections'])
end
saveDir = [saveDir1 datestr(now,'yymmdd') filesep 'manualROI_selections' '\'];
end

%run export function to produce all the maxProj files in the output folder
if reexportFlag == 1
ImBat_exportProjections(saveFlag)
end

cd(saveDir); %go to new folder with the maxcorr images
tifDir = dir('*_maxCorrProject.tif');
ROI = struct;
for i = 1:length(tifDir)
    ROI(i).batName = tifDir(i).name(1:3);
    ROI(i).dateSesh = tifDir(i).name(5:10);
    ROI(i).sessionType = tifDir(i).name(12:16);
    ROI(i).fileName = tifDir(i).name;
    ROI(i).folder = tifDir(i).folder;
    ROI(i).data = FS_annotate_image(tifDir(i).name);
end

for ii = length(tifDir):-1:1
        ROI(ii).data = FS_annotate_image(tifDir(ii).name);
end

save([ROI(1).batName '_ROI_' datestr(now,'YYmmDD_hhMM') '.mat'],'ROI');

%%
trackableROI = figure(); 
sgtitle([ROI(1).batName ': stable cells across ' num2str(length(ROI)) ' days']);
for iii = 1:length(ROI)
    subplot(4,3,iii);%,length(ROI.coordinates),iii); 
    hold on; 
    for p = 1:length(ROI(iii).data.coordinates) 
        plot(ROI(iii).data.coordinates{p}(:,1),ROI(iii).data.coordinates{p}(:,2)); 
        txt = text(ROI(iii).data.coordinates{p}(ceil(end/2),1),ROI(iii).data.coordinates{p}(1,2),num2str(p));
    end
    set(gca, 'YDir','reverse')
    title([ROI(iii).batName ' ' ROI(iii).dateSesh]);
    xticklabels([]);
    yticklabels([]);
    
end

savefig(trackableROI,['plot_' ROI(1).batName '_ROI_' datestr(now,'YYmmDD_hhMM') '.fig']);
saveas(trackableROI,['plot_' ROI(1).batName '_ROI_' datestr(now,'YYmmDD_hhMM') '.tif']);

%%

% dirTop = dir('Ga*');
% for day_i = 1:length(nDays) %for each day
%     %load results data
%     try %extract metadata names and enter processed folder
%         cd([dirTop(nDays(day_i)).name filesep 'extracted'])
%         flyFolders = dir('*fly*extraction');
%         batName{day_i} = flyFolders(end).name(1:3);
%         dateSesh{day_i} = flyFolders(end).name(5:10);
%         sessionType{day_i} = flyFolders(end).name(12:16);
%         
%         cd(flyFolders(end).name);
%         %copy _maxCorrProject files to working dir
%         dirAnalysis = dir('analysis_*');
%         cd([dirAnalysis(end).name filesep 'ROI']);
%         maxCorrDir = dir('*_maxCorrProject.tif');
%         
%     fp = dir('*flightPaths.mat');
%     load(fp(end).name); %load flightPaths
%     sd = dir('*snakePlotData.mat');
%     load(sd(end).name);

