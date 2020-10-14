function ImBat_functionalStability_allDays(dataPreDurPost)
batId = 'Gen';
saveFlag = 1;
stdFlag = 0;
cRaw = 0;
%load first day placeCellStableROI data
%cd([dirTop(day1).folder filesep 'plots\200911-preDurPost cells across days']);
%if cRaw == 1
%load('200311to200320_activity_allTrials_Gal.mat'); %load activity for pre,dur,post
%load('200311to200320_Gal_activity_allTrials_allClusts_sMat.mat'); %load s matrix activity data
saveTag = ['selectiveCells_consistent_allDays_sMat_std'];
if strcmp(batId,'Gal')
    dirDates = dir('Gal_*');
elseif strcmp(batId,'Gen')
    dirDates = dir('Gen_*');
end
nDays = length(dirDates);
day1 = 1; %first day to start comparing
dayEnd = nDays; %last day to start comparing
days = [1:5]; %all days concat
clusts = 2;
preDur = 90; %pre takeoff
postDur = 210; %post takeoff
if strcmp(batId,'Gal')
    clustList = [2 2 2 2 2 2 2 2 2; 3 5 3 3 3 3 3 3 5; 1 3 4 6 4 11 1 1 8];
elseif strcmp(batId,'Gen')
    clustList = [2 2 2 2 2; 3 3 3 3 3];
end


if saveFlag == 1
    saveDir1 = '\\169.229.54.11\server_home\users\tobias\flight\data_processed\topQualityData\analysis_done\plots\';
    % Check if folder exists
    if exist([saveDir1 datestr(now,'yymmdd') filesep 'functionalStability'])>0;
        disp('Youve been working today..');
    else
        mkdir([saveDir1 datestr(now,'yymmdd') filesep 'functionalStability'])
    end
    saveDir = [saveDir1 datestr(now,'yymmdd') filesep 'functionalStability' '\'];
end
%
% %pull in/cutout all the aligned, pre, dur, and post data
% minFlightLength = min(cellfun('size',activity_allTrials.act_dur{clusts}(:,1),2)); %shortest flight length for each day since they may be slightly different lengths
% flightLength = minFlightLength - preDur - postDur;
% act_aligned = cellfun(@(c) c(:,1:minFlightLength), activity_allTrials.act_dur{clusts},'UniformOutput',false); %cut the traces so they are aligned to the shortest flight trajectory across all the days
% act_pre = cellfun(@(c) c(:,1:preDur), activity_allTrials.act_dur{clusts},'UniformOutput',false); %cut the traces so they are aligned to the shortest flight trajectory across all the days
% act_dur = cellfun(@(c) c(:,preDur+1:minFlightLength-postDur), activity_allTrials.act_dur{clusts},'UniformOutput',false); %cut the traces so they are aligned to the shortest flight trajectory across all the days
% act_post = cellfun(@(c) c(:,minFlightLength-postDur+1:minFlightLength), activity_allTrials.act_dur{clusts},'UniformOutput',false); %cut the traces so they are aligned to the shortest flight trajectory across all the days
%
% %initialize matrices to hold mean and stds of each day/roi
% meanMax_act_aligned = zeros(size(days,2),size(act_aligned,2));
% stdMax_act_aligned = zeros(size(days,2),size(act_aligned,2));
% meanMax_act_pre = zeros(size(days,2),size(act_aligned,2));
% stdMax_act_pre = zeros(size(days,2),size(act_aligned,2));
% meanMax_act_dur = zeros(size(days,2),size(act_aligned,2));
% stdMax_act_dur = zeros(size(days,2),size(act_aligned,2));
% meanMax_act_post = zeros(size(days,2),size(act_aligned,2));
% stdMax_act_post = zeros(size(days,2),size(act_aligned,2));
% ylab = cell(1,3);
% overlapSorted = cell(size(days,2),3);
% intersectDays = cell(1,3);
%
% %determine the max,mean, and std for the pre, dur, post data
% for day_i = 1:length(days)
%     for roi_i = 1:size(act_aligned,2)
%         for flight_i = 1:size(act_aligned{days(day_i),roi_i},1)
%             [~,max_act_aligned{day_i,roi_i}(flight_i,1)] = max(act_aligned{days(day_i),roi_i}(flight_i,:));
%             [~,max_act_pre{day_i,roi_i}(flight_i,1)] = max(act_pre{days(day_i),roi_i}(flight_i,:));
%             [~,max_act_dur{day_i,roi_i}(flight_i,1)] = max(act_dur{days(day_i),roi_i}(flight_i,:));
%             [~,max_act_post{day_i,roi_i}(flight_i,1)] = max(act_post{days(day_i),roi_i}(flight_i,:));
%         end
%         meanMax_act_aligned(day_i,roi_i) = round(mean(max_act_aligned{day_i,roi_i}));
%         stdMax_act_aligned(day_i,roi_i) = round(std(max_act_aligned{day_i,roi_i}));
%         meanMax_act_pre(day_i,roi_i) = round(mean(max_act_pre{day_i,roi_i}));
%         stdMax_act_pre(day_i,roi_i) = round(std(max_act_pre{day_i,roi_i}));
%         meanMax_act_dur(day_i,roi_i) = round(mean(max_act_dur{day_i,roi_i}));
%         stdMax_act_dur(day_i,roi_i) = round(std(max_act_dur{day_i,roi_i}));
%         meanMax_act_post(day_i,roi_i) = round(mean(max_act_post{day_i,roi_i}));
%         stdMax_act_post(day_i,roi_i) = round(std(max_act_post{day_i,roi_i}));
%     end
%     %sort the means of maxes for each full, pre, dur, and post
%     [BmeanMax_act_aligned(day_i,:),ImeanMax_act_aligned(day_i,:)] = sort(meanMax_act_aligned(day_i,:));
%     [BmeanMax_act_pre(day_i,:),ImeanMax_act_pre(day_i,:)] = sort(meanMax_act_pre(day_i,:));
%     [BmeanMax_act_dur(day_i,:),ImeanMax_act_dur(day_i,:)] = sort(meanMax_act_dur(day_i,:));
%     [BmeanMax_act_post(day_i,:),ImeanMax_act_post(day_i,:)] = sort(meanMax_act_post(day_i,:));
%
%     %load the stable place cell data for each day
%     cd(dirGalDates(days(day_i)).name);
%     dirExtracted = dir('*Extracted_trajectories*');
%     load(dirExtracted(end).name);
%     cd(dirGalDates(days(day_i)).folder);
%     %pullout the selectively active cells
%     selectiveCells{day_i,1} = placeCellsAngStable.ppre_cells;
%     selectiveCells{day_i,2} = placeCellsAngStable.pp_cells;
%     selectiveCells{day_i,3} = placeCellsAngStable.ppost_cells;
%
%     %sort the means of the maxes so that they are in chronological order
%     %across ROIs
%     for phase_i = 1:size(selectiveCells,2)
%         findSelectiveCells{day_i,phase_i} = find(selectiveCells{day_i,phase_i}(clusts,:));
%         if phase_i == 1
%             [BSelectiveCells{day_i,phase_i},ISelectiveCells{day_i,phase_i}] = sort(meanMax_act_pre(day_i,findSelectiveCells{day_i,phase_i}));
%         elseif phase_i == 2
%             [BSelectiveCells{day_i,phase_i},ISelectiveCells{day_i,phase_i}] = sort(meanMax_act_dur(day_i,findSelectiveCells{day_i,phase_i}));
%         elseif phase_i == 3
%             [BSelectiveCells{day_i,phase_i},ISelectiveCells{day_i,phase_i}] = sort(meanMax_act_post(day_i,findSelectiveCells{day_i,phase_i}));
%         end
%     end
% end
%
% %find ROIs that intersect across all days
% for day_i = 1:size(selectiveCells,1)
%     for phase_i = 1:size(selectiveCells,2)
%         intersectDays{phase_i} = mintersect(ISelectiveCells{:,phase_i}); %all intersecting ROIs across all days
%         roiNum = 1; %counter
%         for roi_i = 1:length(ISelectiveCells{day_i,phase_i})
%             if ismember(ISelectiveCells{day_i,phase_i}(roi_i),intersectDays{phase_i})
%                 overlapSorted{day_i,phase_i}(roiNum) = ISelectiveCells{day_i,phase_i}(roi_i); %sort each intersecting ROI based on the sorting from the original data
%                 if day_i ==1
%                     ylab{phase_i} = [ylab{phase_i} overlapSorted{1,phase_i}(roiNum)];%make label for the yaxis for the plots
%                 end
%                 roiNum = roiNum +1;
%             end
%         end
%     end
% end
%% plot only consistent cells
if strcmp(batId,'Gal')
    clustList = [2 2 2 2 2 2 2 2 2; 3 5 3 3 3 3 3 3 5; 1 3 4 6 4 11 1 1 8];
elseif strcmp(batId,'Gen')
    clustList = [2 2 2 2 2; 3 3 3 3 3];
end
functionalStabilityDotPlot_consistent_all = figure();
set(gcf, 'units','normalized','outerposition',[0 0 0.8 0.6]);
for clust_i = 1:size(clustList,1);
    sgtitle([batId ': Functional Stability of Consistently Selectively Active Cells: All Days: Cluster ' num2str(clust_i+1)]);
    colorDays = jet(length(days)); %make the color vector for each day
    for day_i = 1:size(dataPreDurPost.selectiveCells,1)
        roiStableCounter = length([dataPreDurPost.intersectDays{clust_i+1}{:}]); %roi counter
        for phase_i = 1:size(dataPreDurPost.selectiveCells,2)
            for roi_i = 1:length([dataPreDurPost.intersectDays{clust_i+1}{phase_i}])
                hold on;
                %plot the mean of the max and std for each ROI
                if phase_i == 1
                    if stdFlag ==1
                        p1(day_i) = plot(dataPreDurPost.BSelectiveCells{clustList(clust_i,day_i)}{day_i,phase_i}(dataPreDurPost.intersectDays{clust_i+1}{phase_i}(roi_i)),roiStableCounter-str2num(['0.' num2str(day_i)]),'.','Color',[colorDays(day_i,:) 0.5],'MarkerSize',20);
                        p2 = plot(dataPreDurPost.BSelectiveCells{clustList(clust_i,day_i)}{day_i,phase_i}(dataPreDurPost.intersectDays{clust_i+1}{phase_i}(roi_i))-round(dataPreDurPost.stdMax_act_pre{clustList(clust_i,day_i)}(day_i,roi_i)/2):1:dataPreDurPost.BSelectiveCells{clustList(clust_i,day_i)}{day_i,phase_i}(dataPreDurPost.intersectDays{clust_i+1}{phase_i}(roi_i))+round(dataPreDurPost.stdMax_act_pre{clustList(clust_i,day_i)}(day_i,roi_i)/2),...
                            (ones(round(dataPreDurPost.stdMax_act_pre{clustList(clust_i,day_i)}(day_i,roi_i)/2)*2)*roiStableCounter)-(ones(round(dataPreDurPost.stdMax_act_pre{clustList(clust_i,day_i)}(day_i,roi_i)/2)*2)*str2num(['0.' num2str(day_i)])) ,'-','Color',[colorDays(day_i,:) 0.5],'LineWidth',1);
                    else
                        p1(day_i) = plot(dataPreDurPost.BSelectiveCells{clustList(clust_i,day_i)}{day_i,phase_i}(dataPreDurPost.intersectDays{clust_i+1}{phase_i}(roi_i)),roiStableCounter,'.','Color',[colorDays(day_i,:) 0.5],'MarkerSize',20);
                    end
                elseif phase_i == 2
                    if stdFlag ==1
                        p1(day_i) = plot(preDur + dataPreDurPost.BSelectiveCells{clustList(clust_i,day_i)}{day_i,phase_i}(dataPreDurPost.intersectDays{clust_i+1}{phase_i}(roi_i)),roiStableCounter-str2num(['0.' num2str(day_i)]),'.','Color',[colorDays(day_i,:) 0.5],'MarkerSize',20);
                        p2 = plot(1+preDur+dataPreDurPost.BSelectiveCells{clustList(clust_i,day_i)}{day_i,phase_i}(dataPreDurPost.intersectDays{clust_i+1}{phase_i}(roi_i))-round(stdMax_act_dur{clustList(clust_i,day_i)}(day_i,roi_i)/2):1:preDur+dataPreDurPost.BSelectiveCells{clustList(clust_i,day_i)}{day_i,phase_i}(dataPreDurPost.intersectDays{clust_i+1}{phase_i}(roi_i))+round(dataPreDurPost.stdMax_act_dur{clustList(clust_i,day_i)}(day_i,roi_i)/2),...
                            (ones(round(dataPreDurPost.stdMax_act_dur{clustList(clust_i,day_i)}(day_i,roi_i)/2)*2)*roiStableCounter)-(ones(round(dataPreDurPost.stdMax_act_dur{clustList(clust_i,day_i)}(day_i,roi_i)/2)*2)*str2num(['0.' num2str(day_i)])) ,'-','Color',[colorDays(day_i,:) 0.5],'LineWidth',1);
                    else
                        p1(day_i) = plot(preDur + dataPreDurPost.BSelectiveCells{clustList(clust_i,day_i)}{day_i,phase_i}(dataPreDurPost.intersectDays{clust_i+1}{phase_i}(roi_i)),roiStableCounter,'.','Color',[colorDays(day_i,:) 0.5],'MarkerSize',20);
                    end
                elseif phase_i == 3
                    if stdFlag ==1
                        p1(day_i) = plot(preDur + dataPreDurPost.flightLength{clustList(clust_i,day_i)} + dataPreDurPost.BSelectiveCells{clustList(clust_i,day_i)}{day_i,phase_i}(dataPreDurPost.intersectDays{clust_i+1}{phase_i}(roi_i)),roiStableCounter-str2num(['0.' num2str(day_i)]),'.','Color',[colorDays(day_i,:) 0.5],'MarkerSize',20);
                        p2 = plot(1+preDur+dataPreDurPost.flightLength{clustList(clust_i,day_i)}+dataPreDurPost.BSelectiveCells{clustList(clust_i,day_i)}{day_i,phase_i}(dataPreDurPost.intersectDays{clust_i+1}{phase_i}(roi_i))-round(dataPreDurPost.stdMax_act_post{clustList(clust_i,day_i)}(day_i,roi_i)/2):1:preDur+dataPreDurPost.flightLength{clustList(clust_i,day_i)}+dataPreDurPost.BSelectiveCells{clustList(clust_i,day_i)}{day_i,phase_i}(dataPreDurPost.intersectDays{clust_i+1}{phase_i}(roi_i))+round(dataPreDurPost.stdMax_act_post{clustList(clust_i,day_i)}(day_i,roi_i)/2),...
                            (ones(round(dataPreDurPost.stdMax_act_post{clustList(clust_i,day_i)}(day_i,roi_i)/2)*2)*roiStableCounter)-(ones(round(dataPreDurPost.stdMax_act_post{clustList(clust_i,day_i)}(day_i,roi_i)/2)*2)*str2num(['0.' num2str(day_i)])) ,'-','Color',[colorDays(day_i,:) 0.5],'LineWidth',1);
                    else
                        p1(day_i) = plot(preDur + dataPreDurPost.flightLength{clustList(clust_i,day_i)} + dataPreDurPost.BSelectiveCells{clustList(clust_i,day_i)}{day_i,phase_i}(dataPreDurPost.intersectDays{clust_i+1}{phase_i}(roi_i)),roiStableCounter,'.','Color',[colorDays(day_i,:) 0.5],'MarkerSize',20);
                    end
                end
                roiStableCounter = roiStableCounter - 1;
                hold on;
            end
            flipYlab{phase_i} = flip(dataPreDurPost.ylab{clustList(clust_i,day_i)}{phase_i}); %flip the ylabels so they are top to bottom in order
        end
    end
    sumOverlaps = length([dataPreDurPost.intersectDays{clust_i+1}{:}]); %total number of ROIs plotted
    %make dotted line for pre and post marker
    xPre=preDur*ones(1,sumOverlaps*2+1);
    xPost=(preDur+dataPreDurPost.flightLength{clustList(clust_i,day_i)})*ones(1,sumOverlaps*2+1);
    plot(xPre,1:0.5:sumOverlaps+1,'k.');
    plot(xPost,1:0.5:sumOverlaps+1,'k.');
    axis([0 dataPreDurPost.minFlightLength{clustList(clust_i,day_i)} 0 sumOverlaps+1]);
    set(gca,'YTick',[1:sumOverlaps],'yticklabel',[dataPreDurPost.ylab{clustList(clust_i,day_i)}{:}]);
    xt = get(gca, 'XTick');
    set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
    xlabel('Time (s)');
    ylabel('ROI #');
    try
        if strcmp(batId,'Gen')
        legend([p1(1),p1(2),p1(3),p1(4),p1(5)],'Day 1','Day 2','Day 3','Day 4','Day 5');
        elseif strcmp(batId,'Gal')
            legend([p1(1),p1(2),p1(3),p1(4),p1(5),p1(6),p1(7),p1(8),p1(9)],'Day 1','Day 2','Day 3','Day 4','Day 5','Day 6','Day 7','Day 8','Day 9');
        end
    catch
    end
    
    
    if saveFlag == 1
        if strcmp(batId,'Gal')
            saveas(functionalStabilityDotPlot_consistent_all,[saveDir filesep 'Gal_200311and20_functionalStability_dotPlot_' saveTag '_clust' num2str(clust_i+1) '_' datestr(now,'yymmdd-HHMM') '.tif']);
            savefig(functionalStabilityDotPlot_consistent_all,[saveDir filesep 'Gal_200311and20_functionalStability_dotPlot_' saveTag '_clust' num2str(clust_i+1) '_' datestr(now,'yymmdd-HHMM') '.fig']);
        elseif strcmp(batId,'Gen')
            saveas(functionalStabilityDotPlot_consistent_all,[saveDir filesep 'Gen_200319and24_functionalStability_dotPlot_' saveTag '_clust' num2str(clust_i+1) '_' datestr(now,'yymmdd-HHMM') '.tif']);
            savefig(functionalStabilityDotPlot_consistent_all,[saveDir filesep 'Gen_200319and24_functionalStability_dotPlot_' saveTag '_clust' num2str(clust_i+1) '_' datestr(now,'yymmdd-HHMM') '.fig']);
            
        end
        
    end
    
    %zoom in and resave
    axis([preDur+dataPreDurPost.flightLength{clustList(clust_i,day_i)} 360 0 sumOverlaps+1]); %dataPreDurPost.minFlightLength{clustList(clust_i,day_i)} 0 sumOverlaps+1]);
    drawnow;
    if saveFlag == 1
        if strcmp(batId,'Gal')
            saveas(functionalStabilityDotPlot_consistent_all,[saveDir filesep 'Gal_200311and20_functionalStability_dotPlot_' saveTag '_zoom_clust' num2str(clust_i+1) '_' datestr(now,'yymmdd-HHMM') '.tif']);
            savefig(functionalStabilityDotPlot_consistent_all,[saveDir filesep 'Gal_200311and20_functionalStability_dotPlot_' saveTag '_zoom_clust' num2str(clust_i+1) '_' datestr(now,'yymmdd-HHMM') '.fig']);
        elseif strcmp(batId,'Gen')
            saveas(functionalStabilityDotPlot_consistent_all,[saveDir filesep 'Gen_200319and24_functionalStability_dotPlot_' saveTag '_clust' num2str(clust_i+1) '_' datestr(now,'yymmdd-HHMM') '.tif']);
            savefig(functionalStabilityDotPlot_consistent_all,[saveDir filesep 'Gen_200319and24_functionalStability_dotPlot_' saveTag '_clust' num2str(clust_i+1) '_' datestr(now,'yymmdd-HHMM') '.fig']);
            
        end
    end
    clf;
end
