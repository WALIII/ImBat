ROI_refined = ROI_refined;% ROI_refined_Gal_11to20_c65;
saveFlag = 1; %do you want to save the figures and output structure?
saveDir1 = '\\169.229.54.11\server_home\users\tobias\flight\data_processed\topQualityData\analysis_done\plots\';
% Check if folder exists
if exist([saveDir1 datestr(now,'yymmdd') filesep 'ROI_covar_scatter'])>0;
    disp('Youve been working today..');
else
    mkdir([saveDir1 datestr(now,'yymmdd') filesep 'ROI_covar_scatter'])
end
saveDir = [saveDir1 datestr(now,'yymmdd') filesep 'ROI_covar_scatter' '\'];

scatterCorrDist = figure('units','normalized','outerposition',[0 0.5 1 0.5]); 
for i = 1:length(ROI_refined.batName)
    subplot(2,length(ROI_refined.batName),i)
    scatter(ROI_refined.timeCorrS{i}(ROI_refined.corrIndAllS{i}),ROI_refined.distance{i}(ROI_refined.corrIndAll{i})) 
    hold on;
    sgtitle('Gal Smatrix: Pairwise ROI distance (um) vs corr coeff (R)')
    title(['Day # ' num2str(i)]);
    if i == 1
        xlabel('CorrCoef');
        ylabel('Distance (um)');
    end
end

for i = 1:length(ROI_refined.batName)
    subplot(2,length(ROI_refined.batName),length(ROI_refined.batName)+i)
    scatter(ROI_refined.timeCorrCloseS{i}(ROI_refined.corrIndCloseS{i}),ROI_refined.distance{i}(ROI_refined.corrIndClose{i}),'r') 
    hold on;
    %sgtitle('SMatrix: Pairwise close ROI distance (um) vs corr coeff (R)')
    %title(['Gal Day# ' num2str(i)]);
    if i == 1
        xlabel('CorrCoef');
        ylabel('Subthresh: Distance (um)');
    end
end

if saveFlag == 1
    %save fig and tif of max projection
    %set(findall(maxFig,'-property','FontSize'),'FontSize',20);
    savefig(scatterCorrDist,[saveDir 'scatterCorrDist_' ROI_refined.batName{1} '_19to24_' datestr(now,'yymmdd-hhMMss') '.fig']);
    saveas(scatterCorrDist, [saveDir 'scatterCorrDist_' ROI_refined.batName{1} '_19to24_' datestr(now,'yymmdd-hhMMss') '.tif']);
    %saveas(scatterCorrDist, [saveDir 'scatterCorrDist' ROI_refined.batName{1} ROI_refined.dateSesh{1} ROI_refined.sessionType{1} '_' datestr(now,'yymmdd-hhMMss') '.svg']); 
end