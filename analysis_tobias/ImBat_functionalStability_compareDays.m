function ImBat_functionalStability_compareDays(batId)

saveFlag = 1;

%load first day placeCellStableROI data
%cd([dirTop(day1).folder filesep 'plots\200911-preDurPost cells across days']);
if strcmp(batId,'Gal')
load('200311to200320_Gal_activity_allTrials_allClusts_sMat_vel.mat'); %load s matrix activity data
dirDates = dir('Gal_*');
elseif strcmp(batId,'Gen')
load('200319to200324_Gen_activity_allTrials_allClusts_sMat_vel.mat'); %load s matrix activity data
dirDates = dir('Gen_*');
end
%load('200311to200320_activity_allTrials_Gal.mat'); %load activity for pre,
%dur, post with craw data
nDays = length(dirDates);
day1 = 1;
dayEnd = 2;
days = [day1 dayEnd];
clusts = 2;
preDur = 90;
postDur = 210;


if saveFlag == 1
    saveTag = [num2str(day1) 'v' num2str(dayEnd)];
    saveDir1 = '\\169.229.54.11\server_home\users\tobias\flight\data_processed\topQualityData\analysis_done\plots\';
    % Check if folder exists
    if exist([saveDir1 datestr(now,'yymmdd') filesep 'functionalStability'])>0;
        disp('Youve been working today..');
    else
        mkdir([saveDir1 datestr(now,'yymmdd') filesep 'functionalStability'])
    end
    saveDir = [saveDir1 datestr(now,'yymmdd') filesep 'functionalStability' '\'];
end

cd(dirDates(day1).name);
dirExtracted = dir('*Extracted_trajectories*');
load(dirExtracted(end).name);
selectiveCells1 = placeCellsAngStable;
cd(dirDates(day1).folder);

%load last day PlaceCellStableROI data
cd(dirDates(dayEnd).name);
dirExtracted = dir('*Extracted_trajectories*');
load(dirExtracted(end).name);
selectiveCellsEnd = placeCellsAngStable;
cd(dirDates(day1).folder);

%pullout the selectively active cells
selectiveCells{1,1} = selectiveCells1.ppre_cells;
selectiveCells{1,2} = selectiveCells1.pp_cells;
selectiveCells{1,3} = selectiveCells1.ppost_cells;
selectiveCells{2,1} = selectiveCellsEnd.ppre_cells;
selectiveCells{2,2} = selectiveCellsEnd.pp_cells;
selectiveCells{2,3} = selectiveCellsEnd.ppost_cells;

% minFlightLength = min(cellfun('size',activity_allTrials.act_dur{clusts}(:,1),2)); %shortest flight length for each day since they may be slightly different lengths
% flightLength = minFlightLength - preDur - postDur;
% mean_act_full = cellfun(@mean,activity_allTrials.act_dur{clusts},'UniformOutput',false); %take the mean of all the trials for each roi/day
% act_aligned = cellfun(@(c) c(:,1:minFlightLength), activity_allTrials.act_dur{clusts},'UniformOutput',false); %cut the traces so they are aligned to the shortest flight trajectory across all the days
% mean_act_aligned = cellfun(@mean,act_aligned{clusts},'UniformOutput',false); %take the mean of the aligned traces
% std_act_aligned = cellfun(@std,act_aligned{clusts},'UniformOutput',false); %take the mean of the aligned traces
% 
% act_pre = cellfun(@(c) c(:,1:preDur), activity_allTrials.act_dur{clusts},'UniformOutput',false); %cut the traces so they are aligned to the shortest flight trajectory across all the days
% mean_act_pre = cellfun(@mean,act_pre,'UniformOutput',false); %take the mean of the aligned traces
% std_act_pre = cellfun(@std,act_pre,'UniformOutput',false); %take the mean of the aligned traces
% act_dur = cellfun(@(c) c(:,preDur+1:minFlightLength-postDur), activity_allTrials.act_dur{clusts},'UniformOutput',false); %cut the traces so they are aligned to the shortest flight trajectory across all the days
% mean_act_dur = cellfun(@mean,act_dur,'UniformOutput',false); %take the mean of the aligned traces
% std_act_dur = cellfun(@std,act_dur,'UniformOutput',false); %take the mean of the aligned traces
% act_post = cellfun(@(c) c(:,minFlightLength-postDur+1:minFlightLength), activity_allTrials.act_dur{clusts},'UniformOutput',false); %cut the traces so they are aligned to the shortest flight trajectory across all the days
% mean_act_post = cellfun(@mean,act_post,'UniformOutput',false); %take the mean of the aligned traces
% std_act_post = cellfun(@std,act_post,'UniformOutput',false); %take the mean of the aligned traces
% 
% %sort the max of each mean
% [~,max_mean_act_aligned] = cellfun(@max,mean_act_aligned);
% [BMax_mean_act_aligned,IMax_mean_act_aligned] = sort(max_mean_act_aligned,2);
% [~,max_mean_act_pre] = cellfun(@max,mean_act_pre);
% [BMax_mean_act_pre,IMax_mean_act_pre] = sort(max_mean_act_pre,2);
% [~,max_mean_act_dur] = cellfun(@max,mean_act_dur);
% [BMax_mean_act_dur,IMax_mean_act_dur] = sort(max_mean_act_dur,2);
% [~,max_mean_act_post] = cellfun(@max,mean_act_post);
% [BMax_mean_act_post,IMax_mean_act_post] = sort(max_mean_act_post,2);

%pull in/cutout all the aligned, pre, dur, and post data
minFlightLength = min(cellfun('size',activity_allTrials.act_dur{clusts}(:,1),2)); %shortest flight length for each day since they may be slightly different lengths
flightLength = minFlightLength - preDur - postDur;
act_aligned = cellfun(@(c) c(:,1:minFlightLength), activity_allTrials.act_dur{clusts},'UniformOutput',false); %cut the traces so they are aligned to the shortest flight trajectory across all the days
act_pre = cellfun(@(c) c(:,1:preDur), activity_allTrials.act_dur{clusts},'UniformOutput',false); %cut the traces so they are aligned to the shortest flight trajectory across all the days
act_dur = cellfun(@(c) c(:,preDur+1:minFlightLength-postDur), activity_allTrials.act_dur{clusts},'UniformOutput',false); %cut the traces so they are aligned to the shortest flight trajectory across all the days
act_post = cellfun(@(c) c(:,minFlightLength-postDur+1:minFlightLength), activity_allTrials.act_dur{clusts},'UniformOutput',false); %cut the traces so they are aligned to the shortest flight trajectory across all the days

%initialize matrices to hold mean and stds of each day/roi
meanMax_act_aligned = zeros(size(days,2),size(act_aligned,2));
stdMax_act_aligned = zeros(size(days,2),size(act_aligned,2));
meanMax_act_pre = zeros(size(days,2),size(act_aligned,2));
stdMax_act_pre = zeros(size(days,2),size(act_aligned,2));
meanMax_act_dur = zeros(size(days,2),size(act_aligned,2));
stdMax_act_dur = zeros(size(days,2),size(act_aligned,2));
meanMax_act_post = zeros(size(days,2),size(act_aligned,2));
stdMax_act_post = zeros(size(days,2),size(act_aligned,2));

% for day_i = 1:size(selectiveCells,1)
%     for phase_i = 1:size(selectiveCells,2)
%         findSelectiveCells{day_i,phase_i} = find(selectiveCells{day_i,phase_i}(2,:));
%         if phase_i == 1
%             [BSelectiveCells{day_i,phase_i},ISelectiveCells{day_i,phase_i}] = sort(max_mean_act_pre(days(day_i),findSelectiveCells{day_i,phase_i}));
%         elseif phase_i == 2
%             [BSelectiveCells{day_i,phase_i},ISelectiveCells{day_i,phase_i}] = sort(max_mean_act_dur(days(day_i),findSelectiveCells{day_i,phase_i}));
%         elseif phase_i == 3
%             [BSelectiveCells{day_i,phase_i},ISelectiveCells{day_i,phase_i}] = sort(max_mean_act_post(days(day_i),findSelectiveCells{day_i,phase_i}));
%         end
%     end
%     sumSelectiveCells(day_i) = length(findSelectiveCells{day_i,1}) + length(findSelectiveCells{day_i,2}) + length(findSelectiveCells{day_i,3});
% end

%determine the max,mean, and std for the pre, dur, post data
for day_i = 1:length(days)
    for roi_i = 1:size(act_aligned,2)
        for flight_i = 1:size(act_aligned{days(day_i),roi_i},1)
            [~,max_act_aligned{day_i,roi_i}(flight_i,1)] = max(act_aligned{days(day_i),roi_i}(flight_i,:));
            [~,max_act_pre{day_i,roi_i}(flight_i,1)] = max(act_pre{days(day_i),roi_i}(flight_i,:));
            [~,max_act_dur{day_i,roi_i}(flight_i,1)] = max(act_dur{days(day_i),roi_i}(flight_i,:));
            [~,max_act_post{day_i,roi_i}(flight_i,1)] = max(act_post{days(day_i),roi_i}(flight_i,:));
        end
        meanMax_act_aligned(day_i,roi_i) = round(mean(max_act_aligned{day_i,roi_i}));
        stdMax_act_aligned(day_i,roi_i) = round(std(max_act_aligned{day_i,roi_i}));
        meanMax_act_pre(day_i,roi_i) = round(mean(max_act_pre{day_i,roi_i}));
        stdMax_act_pre(day_i,roi_i) = round(std(max_act_pre{day_i,roi_i}));
        meanMax_act_dur(day_i,roi_i) = round(mean(max_act_dur{day_i,roi_i}));
        stdMax_act_dur(day_i,roi_i) = round(std(max_act_dur{day_i,roi_i}));
        meanMax_act_post(day_i,roi_i) = round(mean(max_act_post{day_i,roi_i}));
        stdMax_act_post(day_i,roi_i) = round(std(max_act_post{day_i,roi_i}));
    end
    %sort the means of maxes for each full, pre, dur, and post
    [BmeanMax_act_aligned(day_i,:),ImeanMax_act_aligned(day_i,:)] = sort(meanMax_act_aligned(day_i,:));
    [BmeanMax_act_pre(day_i,:),ImeanMax_act_pre(day_i,:)] = sort(meanMax_act_pre(day_i,:));
    [BmeanMax_act_dur(day_i,:),ImeanMax_act_dur(day_i,:)] = sort(meanMax_act_dur(day_i,:));
    [BmeanMax_act_post(day_i,:),ImeanMax_act_post(day_i,:)] = sort(meanMax_act_post(day_i,:));
%sort the means of the maxes so that they are in chronological order
    %across ROIs
    for phase_i = 1:size(selectiveCells,2)
        findSelectiveCells{day_i,phase_i} = find(selectiveCells{day_i,phase_i}(clusts,:));
        if phase_i == 1
            [BSelectiveCells{day_i,phase_i},ISelectiveCells{day_i,phase_i}] = sort(meanMax_act_pre(day_i,findSelectiveCells{day_i,phase_i}));
        elseif phase_i == 2
            [BSelectiveCells{day_i,phase_i},ISelectiveCells{day_i,phase_i}] = sort(meanMax_act_dur(day_i,findSelectiveCells{day_i,phase_i}));
        elseif phase_i == 3
            [BSelectiveCells{day_i,phase_i},ISelectiveCells{day_i,phase_i}] = sort(meanMax_act_post(day_i,findSelectiveCells{day_i,phase_i}));
        end
    end
    sumSelectiveCells(day_i) = length(findSelectiveCells{day_i,1}) + length(findSelectiveCells{day_i,2}) + length(findSelectiveCells{day_i,3});

end
%% plot the max peak of each selective pre/dur/post cell
functionalStabilityDotPlot = figure();
set(gcf, 'units','normalized','outerposition',[0 0 0.8 0.6]);
sgtitle('Functional Stability of Selectively Active Cells: Cluster 2');
colorDays = ["r","b"];

for day_i = 1:size(selectiveCells,1)
    roiCounter = sumSelectiveCells(day_i);
    for phase_i = 1:size(selectiveCells,2)
        %findSelectiveCells{day_i,phase_i} = find(selectiveCells{day_i,phase_i}(2,:));
        for roi_i = 1:size(findSelectiveCells{day_i,phase_i},2)
            subplot(1,3,day_i); %plot the individual day 1 and day 9
            if phase_i == 1
                p1 = plot(BSelectiveCells{day_i,phase_i}(roi_i),roiCounter,[char(colorDays(day_i)) '.'],'MarkerSize',20);
                %plot(BSelectiveCells{day_i,phase_i}(roi_i),roiCounter,[char(colorDays(day_i)) '.'],'MarkerSize',20);
                %plot(max_mean_act_pre(days(day_i),findSelectiveCells{day_i,phase_i}(roi_i)),roiCounter,[char(colorDays(day_i)) '.'],'MarkerSize',20);
            elseif phase_i == 2
                p1 = plot(preDur+BSelectiveCells{day_i,phase_i}(roi_i),roiCounter,[char(colorDays(day_i)) '.'],'MarkerSize',20);
                %plot(preDur + max_mean_act_dur(days(day_i),findSelectiveCells{day_i,phase_i}(roi_i)), roiCounter,[char(colorDays(day_i)) '.'],'MarkerSize',20);
            elseif phase_i == 3
                p1 = plot(preDur + flightLength + BSelectiveCells{day_i,phase_i}(roi_i),roiCounter,[char(colorDays(day_i)) '.'],'MarkerSize',20);
                %plot(preDur + flightLength + max_mean_act_post(days(day_i),findSelectiveCells{day_i,phase_i}(roi_i)), roiCounter,[char(colorDays(day_i)) '.'],'MarkerSize',20);
            end
            hold on;
            xPre=preDur*ones(1,max(sumSelectiveCells)*2+1);
            xPost=(preDur+flightLength)*ones(1,max(sumSelectiveCells)*2+1);
            plot(xPre,1:0.5:max(sumSelectiveCells)+1,'k.');
            plot(xPost,1:0.5:max(sumSelectiveCells)+1,'k.');
            title(['Day ' num2str(days(day_i))]);
            axis([0 minFlightLength 0 max(sumSelectiveCells)+1]);
            set(gca,'YTick',[1:sumSelectiveCells(day_i)],'yticklabel',flip([ISelectiveCells{day_i,:}]),'ycolor',char(colorDays(day_i)));
            xt = get(gca, 'XTick');
            set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
            xlabel('Time (s)');
            ylabel('ROI #');
            
            subplot(1,3,3); %plot the overlap
            colororder({'r','b'})
            if day_i == 1
                yyaxis left
            else
                yyaxis right
            end
            if phase_i == 1
                p1 = plot(BSelectiveCells{day_i,phase_i}(roi_i),roiCounter,[char(colorDays(day_i)) '.'],'MarkerSize',20);
                %plot(max_mean_act_pre(days(day_i),findSelectiveCells{day_i,phase_i}(roi_i)), roiCounter,[char(colorDays(day_i)) '.'],'MarkerSize',20);
            elseif phase_i == 2
                p1 = plot(preDur+BSelectiveCells{day_i,phase_i}(roi_i),roiCounter,[char(colorDays(day_i)) '.'],'MarkerSize',20);
                %plot(preDur + max_mean_act_dur(days(day_i),findSelectiveCells{day_i,phase_i}(roi_i)), roiCounter,[char(colorDays(day_i)) '.'],'MarkerSize',20);
            elseif phase_i == 3
                plot(preDur + flightLength + BSelectiveCells{day_i,phase_i}(roi_i),roiCounter,[char(colorDays(day_i)) '.'],'MarkerSize',20);
                %plot(preDur + flightLength + max_mean_act_post(days(day_i),findSelectiveCells{day_i,phase_i}(roi_i)), roiCounter,[char(colorDays(day_i)) '.'],'MarkerSize',20);
            end
            hold on;
            plot(xPre,1:0.5:max(sumSelectiveCells)+1,'k.');
            plot(xPost,1:0.5:max(sumSelectiveCells)+1,'k.');
            title(['Overlap']);
            axis([0 minFlightLength 0 max(sumSelectiveCells)+1]);
            set(gca,'yticklabel',[]);
            xlabel('Time (s)');
            ylabel('ROI #');
            set(gca,'YTick',[1:sumSelectiveCells(day_i)],'yticklabel',flip([ISelectiveCells{day_i,:}]));
            xt = get(gca, 'XTick');
            set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));

            roiCounter = roiCounter - 1;
        end
    end
end

if saveFlag == 1
    if strcmp(batId,'Gal')
    saveas(functionalStabilityDotPlot,[saveDir filesep 'Gal_200311and20_functionalStability_dotPlot_selectiveCells_clust' num2str(clusts) '_' saveTag '_' datestr(now,'yymmdd-HHMMSS') '.tif']);
    savefig(functionalStabilityDotPlot,[saveDir filesep 'Gal_200311and20_functionalStability_dotPlot_selectiveCells_clust' num2str(clusts) '_' saveTag '_' datestr(now,'yymmdd-HHMMSS') '.fig']);
    elseif strcmp(batId,'Gen')
    saveas(functionalStabilityDotPlot,[saveDir filesep 'Gen_200319and24_functionalStability_dotPlot_selectiveCells_clust' num2str(clusts) '_' saveTag '_' datestr(now,'yymmdd-HHMMSS') '.tif']);
    savefig(functionalStabilityDotPlot,[saveDir filesep 'Gen_200319and24_functionalStability_dotPlot_selectiveCells_clust' num2str(clusts) '_' saveTag '_' datestr(now,'yymmdd-HHMMSS') '.fig']);
    end
end


%% plot only consistent cells

functionalStabilityDotPlot_consistent = figure();
set(gcf, 'units','normalized','outerposition',[0 0 0.8 0.6]);
sgtitle('Functional Stability of Consistently Selectively Active Cells: Cluster 2');
colorDays = ["r","b"];

%find the overlapping cells and phase
for day_i = 1:size(selectiveCells,1)
    for phase_i = 1:size(selectiveCells,2)
        if day_i == 1
            [overlapVal{day_i,phase_i},overlapPos{day_i,phase_i}]=intersect(ISelectiveCells{1,phase_i},ISelectiveCells{2,phase_i});
        elseif day_i == 2
            [overlapVal{day_i,phase_i},overlapPos{day_i,phase_i}]=intersect(ISelectiveCells{2,phase_i},ISelectiveCells{1,phase_i}); %reverse if it's for day2
        end
        [BoverlapPosSorted{day_i,phase_i},IoverlapPosSorted{day_i,phase_i}] = sort(overlapPos{day_i,phase_i});
        ISelectiveCellsOverlap{day_i,phase_i} = ISelectiveCells{day_i,phase_i}(BoverlapPosSorted{day_i,phase_i});
    end
    overlapCat(day_i,:) = vertcat(overlapPos{day_i,:})';
end

%sort out and find day2 order based on day1 sorting
for phase_i = 1:size(selectiveCells,2)
    for roi_i = 1:length(ISelectiveCellsOverlap{1,phase_i})
        overlapPosSorted_dayEnd{phase_i}(roi_i) = find(ismember(ISelectiveCells{2,phase_i},ISelectiveCellsOverlap{1,phase_i}(roi_i)));
    end
end

%plot day 1
sumOverlaps = length(overlapPos{1,1})+length(overlapPos{1,2})+length(overlapPos{1,3});
for plot_i = [1 3]
roiStableCounter = sumOverlaps;
for phase_i = 1:size(selectiveCells,2)
    for roi_i = 1:length(BoverlapPosSorted{1,phase_i})
        subplot(1,3,plot_i); %plot the individual day 1 and overlap
        if phase_i == 1
            p1 = plot(BSelectiveCells{1,phase_i}(BoverlapPosSorted{1,phase_i}(roi_i)),roiStableCounter,[char(colorDays(1)) '.'],'MarkerSize',20);
        elseif phase_i == 2
            p1 = plot(preDur+BSelectiveCells{1,phase_i}(BoverlapPosSorted{1,phase_i}(roi_i)),roiStableCounter,[char(colorDays(1)) '.'],'MarkerSize',20);
        elseif phase_i == 3
            p1 = plot(preDur + flightLength + BSelectiveCells{1,phase_i}(BoverlapPosSorted{1,phase_i}(roi_i)),roiStableCounter,[char(colorDays(1)) '.'],'MarkerSize',20);
        end
        hold on;
        p1.Color(4) = 0.5;
        xPre=preDur*ones(1,sumOverlaps*2+1);
        xPost=(preDur+flightLength)*ones(1,sumOverlaps*2+1);
        plot(xPre,1:0.5:sumOverlaps+1,'k.');
        plot(xPost,1:0.5:sumOverlaps+1,'k.');
        title(['Day ' num2str(days(1))]);
        axis([0 minFlightLength 0 sumOverlaps+1]);
        set(gca,'YTick',[1:sumOverlaps],'yticklabel',flip(vertcat(IoverlapPosSorted{1,:})'));%,'ycolor',char(colorDays(1)));
        xt = get(gca, 'XTick');
        set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
        xlabel('Time (s)');
        ylabel('ROI #');
    roiStableCounter = roiStableCounter - 1;
    end
end
end
%plot the day2
for plot_ii = [2 3]
roiStableCounter = sumOverlaps;
for phase_i = 1:size(selectiveCells,2)
    for roi_i = 1:length(BoverlapPosSorted{1,phase_i})
        subplot(1,3,plot_ii); %plot the individual day 1 and day 9
        if phase_i == 1
            p2 = plot(BSelectiveCells{2,phase_i}(overlapPosSorted_dayEnd{1,phase_i}(roi_i)),roiStableCounter,[char(colorDays(2)) '.'],'MarkerSize',20);
        elseif phase_i == 2
            p2 = plot(preDur+BSelectiveCells{2,phase_i}(overlapPosSorted_dayEnd{1,phase_i}(roi_i)),roiStableCounter,[char(colorDays(2)) '.'],'MarkerSize',20);
        elseif phase_i == 3
            p2 = plot(preDur + flightLength + BSelectiveCells{2,phase_i}(overlapPosSorted_dayEnd{1,phase_i}(roi_i)),roiStableCounter,[char(colorDays(2)) '.'],'MarkerSize',20);
        end
        hold on;
        p2.Color(4) = 0.5;
        xPre=preDur*ones(1,sumOverlaps*2+1);
        xPost=(preDur+flightLength)*ones(1,sumOverlaps*2+1);
        plot(xPre,1:0.5:sumOverlaps+1,'k.');
        plot(xPost,1:0.5:sumOverlaps+1,'k.');
        if plot_ii == 2
            title(['Day ' num2str(days(end))]);
        else
        title('Overlap');
        end
        axis([0 minFlightLength 0 sumOverlaps+1]);
        set(gca,'YTick',[1:sumOverlaps],'yticklabel',flip(vertcat(IoverlapPosSorted{1,:})'));%,'ycolor',char(colorDays(2)));
        xt = get(gca, 'XTick');
        set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
        xlabel('Time (s)');
        ylabel('ROI #');
        roiStableCounter = roiStableCounter - 1;
    end
end
end

if saveFlag == 1
    if strcmp(batId,'Gal')
        saveas(functionalStabilityDotPlot_consistent,[saveDir filesep 'Gal_200311and20_functionalStability_dotPlot_selectiveCells_clust2_consistent_' saveTag '_' datestr(now,'yymmdd-HHMMSS') '.tif']);
        savefig(functionalStabilityDotPlot_consistent,[saveDir filesep 'Gal_200311and20_functionalStability_dotPlot_selectiveCells_clust2_consistent_' saveTag '_' datestr(now,'yymmdd-HHMMSS') '.fig']);
    elseif strcmp(batId,'Gen')
        saveas(functionalStabilityDotPlot_consistent,[saveDir filesep 'Gen_200319and24_functionalStability_dotPlot_selectiveCells_clust2_consistent_' saveTag '_' datestr(now,'yymmdd-HHMMSS') '.tif']);
        savefig(functionalStabilityDotPlot_consistent,[saveDir filesep 'Gen_200319and24_functionalStability_dotPlot_selectiveCells_clust2_consistent_' saveTag '_' datestr(now,'yymmdd-HHMMSS') '.fig']);
    end
end