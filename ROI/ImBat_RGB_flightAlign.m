function ImBat_RGB_flightAlign(batId,fullSeshTag,day1,day2,day3)
%batId = 'Gal';
%fullSeshTag = 1, chooses to look at whole session for 3 day comparison
clustNum = 3;
saveFlag = 1;
saveTag = 'clust3_fulls';
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

if fullSeshTag == 0
    if isempty(day3) %comparing flight aligned max proj
        if day1==day2 %if looking at the cluster flight aligned vs full session
            IM1 = mat2gray(activity_allTrials.maxMeanFrames_dur{clustNum}{day1});
            IM2 = mat2gray(activity_allTrials.YmaxFull{day2});
            plotTitle = [batId ' clust ' num2str(clustNum) ': Day ' num2str(day1) ' (r) v full session (c)'];
        else %if looking at 1 day vs another day only
            IM1 = mat2gray(activity_allTrials.maxMeanFrames_dur{clustNum}{day1});
            IM2 = mat2gray(activity_allTrials.maxMeanFrames_dur{clustNum}{day2});
            plotTitle = [batId ' clust ' num2str(clustNum) ': Day ' num2str(day1) ' (r) v Day ' num2str(day2) ' (c)']
        end
        IM3 = [];
        % Alignment script
        rowTemp = [size(IM1,1),size(IM2,1),size(IM3,1)];
        colTemp = [size(IM1,2),size(IM2,2),size(IM3,2)];
        minRow = min(rowTemp(rowTemp>0));
        minCol = min(colTemp(colTemp>0));
        [IM1_aligned,IM2_aligned] = ImBat_imageAlign(IM1(1:minRow,1:minCol),IM2(1:minRow,1:minCol),IM3);
        IM3_aligned = [];
        % Will RGB script
        [a1,b1] = CaBMI_XMASS(IM1(1:minRow,1:minCol),IM2(1:minRow,1:minCol),IM3);
        [a2,b2] = CaBMI_XMASS(IM1_aligned,IM2_aligned,IM3_aligned);
    else %if looking at comparing 3 days
        IM1 = mat2gray(activity_allTrials.maxMeanFrames_dur{clustNum}{day1});
        IM2 = mat2gray(activity_allTrials.maxMeanFrames_dur{clustNum}{day2});
        IM3 = mat2gray(activity_allTrials.maxMeanFrames_dur{clustNum}{day3});
        minRow = min([size(IM1,1),size(IM2,1),size(IM3,1)]);
        minCol = min([size(IM1,2),size(IM2,2),size(IM3,2)]);
        % Alignment script
        [IM1_aligned,IM2_aligned,IM3_aligned] = ImBat_imageAlign(IM1(1:minRow,1:minCol),IM2(1:minRow,1:minCol),IM3(1:minRow,1:minCol));
        % Will RGB script
        [a1,b1] = CaBMI_XMASS(IM1(1:minRow,1:minCol),IM2(1:minRow,1:minCol),IM3(1:minRow,1:minCol));
        [a2,b2] = CaBMI_XMASS(IM1_aligned,IM2_aligned,IM3_aligned);
        plotTitle = [batId ' clust ' num2str(clustNum) ': Day ' num2str(day1) ' (r) v Day ' num2str(day2) ' (g) v Day ' num2str(day3) ' (b)'];
    end
elseif fullSeshTag == 1
    %if looking at comparing 3 days full sessions
    IM1 = mat2gray(activity_allTrials.YmaxFull{day1});
    IM2 = mat2gray(activity_allTrials.YmaxFull{day2});
    IM3 = mat2gray(activity_allTrials.YmaxFull{day3});
    minRow = min([size(IM1,1),size(IM2,1),size(IM3,1)]);
    minCol = min([size(IM1,2),size(IM2,2),size(IM3,2)]);
    % Alignment script
    [IM1_aligned,IM2_aligned,IM3_aligned] = ImBat_imageAlign(IM1(1:minRow,1:minCol),IM2(1:minRow,1:minCol),IM3(1:minRow,1:minCol));
    % Will RGB script
    [a1,b1] = CaBMI_XMASS(IM1(1:minRow,1:minCol),IM2(1:minRow,1:minCol),IM3(1:minRow,1:minCol));
    [a2,b2] = CaBMI_XMASS(IM1_aligned,IM2_aligned,IM3_aligned);
    plotTitle = [batId ' clust ' num2str(clustNum) ' Full Sesh: Day ' num2str(day1) ' (r) v Day ' num2str(day2) ' (g) v Day ' num2str(day3) ' (b)'];
end
   


plotDay123_unaligned = figure(); 
sgtitle(plotTitle); 
subplot(1,2,1);
image((a1(:,:,:)));
title('un-aligned');
set(gca,'xticklabel',[],'yticklabel',[]);
subplot(1,2,2);
image((a2(:,:,:)));
title('aligned');
set(gca,'xticklabel',[],'yticklabel',[]);

plotDay123_aligned = figure(); 
image((a2(:,:,:)))
title(plotTitle);
set(gca,'xticklabel',[],'yticklabel',[]);
%set(gca,'YDir','normal');
    
if saveFlag == 1
    if strcmp(batId,'Gal')
        saveas(plotDay123_unaligned,[saveDir filesep 'Gal_200311and20_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.tif']);
        savefig(plotDay123_unaligned,[saveDir filesep 'Gal_200311and20_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
        saveas(plotDay123_aligned,[saveDir filesep 'Gal_200311and20_overLap_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.tif']);
        savefig(plotDay123_aligned,[saveDir filesep 'Gal_200311and20_overLap_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
    elseif strcmp(batId,'Gen')
        saveas(plotDay123_unaligned,[saveDir filesep 'Gen_200319and24_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.tif']);
        savefig(plotDay123_unaligned,[saveDir filesep 'Gen_200319and24_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
        saveas(plotDay123_aligned,[saveDir filesep 'Gen_200319and24_overLap_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.tif']);
        savefig(plotDay123_aligned,[saveDir filesep 'Gen_200319and24_overLap_day' num2str(day1) '_day' num2str(day2) '_day' num2str(day3) '_' saveTag '_' datestr(now,'yymmdd-HHMM') '.fig']);
    end
end

%%
% imagesc(imresize(results10.results.Cn,4)); colormap(gray);
% imagesc(imresize(results12.results.Cn,4)); colormap(gray);
% imagesc(imresize(results15.results.Cn,4)); colormap(gray);
% 
% figure();
% [imRGB] = cat(3,results10.results.Cn(1:112,1:153),results12.results.Cn(1:112,1:153),results15.results.Cn(1:112,1:153));
% 
% C = imfuse(results10.results.Cn(1:112,1:153),results12.results.Cn(1:112,1:153),'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
% 
% a = imresize(results10.results.Cn,4);
% b = imresize(results12.results.Cn,4);
% c = imresize(results15.results.Cn,4);
% agray = mat2gray(a);
% bgray = mat2gray(b);
% cgray = mat2gray(c);
% abcRGB = cat(3,agray(1:448,1:612),bgray(1:448,1:612),cgray(1:448,1:612));
% imagesc(abcRGB);

