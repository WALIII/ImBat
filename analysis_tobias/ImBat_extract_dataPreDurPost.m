function dataPreDurPost = ImBat_extract_dataPreDurPost(batId)
saveFlag = 1;
saveTag = 'sMat_vel';
cRaw = 0;
%load first day placeCellStableROI data
if strcmp(batId,'Gal')
    %cd([dirTop(day1).folder filesep 'plots\200911-preDurPost cells across days']);
    if cRaw ==1
        load('200311to200320_Gal_activity_allTrials_allClusts_cRaw_vel.mat'); %load activity for pre,dur,post
    else
        load('200311to200320_Gal_activity_allTrials_allClusts_sMat_vel.mat'); %load s matrix activity data
    end
    dirDates = dir('Gal_*');
elseif strcmp(batId,'Gen')
    if cRaw ==1
        load('200319to200324_Gen_activity_allTrials_allClusts_cRaw_vel.mat'); %load activity for pre,dur,post
    else
        load('200319to200324_Gen_activity_allTrials_allClusts_sMat_vel.mat'); %load s matrix activity data
    end
    dirDates = dir('Gen_*');
    
end
nDays = length(dirDates);
day1 = 1; %first day to start comparing
dayEnd = nDays; %last day to start comparing
days = [1:dayEnd]; %all days concat
clusts = 2;
CNMFe_Fs = 30;
track_Fs = 120;
preDur = 3; %pre takeoff
postDur = 7; %post takeoff
if strcmp(batId,'Gal')
    clustList = [2 2 2 2 2 2 2 2 2; 3 5 3 3 3 3 3 3 5; 1 3 4 6 4 11 1 1 8];
elseif strcmp(batId,'Gen')
    clustList = [2 2 2 2 2; 3 3 3 3 3];
end

minFlightLength = cell(1,max(max(clustList)));
minBehavFlightLength= cell(1,max(max(clustList)));
flightLength= cell(1,max(max(clustList)));
behavLength= cell(1,max(max(clustList)));
act_aligned= cell(1,max(max(clustList)));
act_pre= cell(1,max(max(clustList)));
act_dur= cell(1,max(max(clustList)));
act_post= cell(1,max(max(clustList)));
vel_aligned= cell(1,max(max(clustList)));
vel_pre= cell(1,max(max(clustList)));
vel_dur= cell(1,max(max(clustList)));
vel_post= cell(1,max(max(clustList)));
%pull in/cutout all the aligned, pre, dur, and post data
for clust_i = 1:15%max(max(clustList))
    minFlightSizesVect = cellfun('size',activity_allTrials.act_dur{clust_i}(:,1),2);
    minFlightLength{clust_i} = min(minFlightSizesVect(minFlightSizesVect>0));
    %minFlightLength{clust_i} = min(cellfun('size',activity_allTrials.act_dur{clust_i}(:,1),2)(cellfun('size',activity_allTrials.act_dur{clust_i}(:,1),2)>0)); %shortest flight length for each day since they may be slightly different lengths
    minBehavSizesVect = cellfun('size',activity_allTrials.vel_dur{clust_i}(1,:),2);
    minBehavFlightLength{clust_i} = min(minBehavSizesVect(minBehavSizesVect>0));
    %minBehavFlightLength{clust_i} = min(cellfun('size',activity_allTrials.vel_dur{clust_i}(1,:),2)); %for behavior
    flightLength{clust_i} = minFlightLength{clust_i} - preDur*CNMFe_Fs - postDur*CNMFe_Fs;
    behavLength{clust_i} = minBehavFlightLength{clust_i} - preDur*track_Fs - postDur*track_Fs;
    %try
    act_aligned{clust_i} = cellfun(@(c) c(:,1:minFlightLength{clust_i}), activity_allTrials.act_dur{clust_i},'UniformOutput',false); %cut the traces so they are aligned to the shortest flight trajectory across all the days
    act_pre{clust_i} = cellfun(@(c) c(:,1:preDur*CNMFe_Fs), activity_allTrials.act_dur{clust_i},'UniformOutput',false); %cut the traces so they are aligned to the shortest flight trajectory across all the days
    act_dur{clust_i} = cellfun(@(c) c(:,1+preDur*CNMFe_Fs:minFlightLength{clust_i}-postDur*CNMFe_Fs), activity_allTrials.act_dur{clust_i},'UniformOutput',false); %cut the traces so they are aligned to the shortest flight trajectory across all the days
    act_post{clust_i} = cellfun(@(c) c(:,1+ minFlightLength{clust_i}-postDur*CNMFe_Fs:minFlightLength{clust_i}), activity_allTrials.act_dur{clust_i},'UniformOutput',false); %cut the traces so they are aligned to the shortest flight trajectory across all the days
    vel_aligned{clust_i} = cellfun(@(c) c(:,1:minBehavFlightLength{clust_i}), activity_allTrials.vel_dur{clust_i},'UniformOutput',false); %cut the traces so they are aligned to the shortest flight trajectory across all the days
    vel_pre{clust_i} = cellfun(@(c) c(:,1:preDur*track_Fs), activity_allTrials.vel_dur{clust_i},'UniformOutput',false); %cut the traces so they are aligned to the shortest flight trajectory across all the days
    vel_dur{clust_i} = cellfun(@(c) c(:,1+preDur*track_Fs:minBehavFlightLength{clust_i}-postDur*track_Fs), activity_allTrials.vel_dur{clust_i},'UniformOutput',false); %cut the traces so they are aligned to the shortest flight trajectory across all the days
    vel_post{clust_i} = cellfun(@(c) c(:,1+minBehavFlightLength{clust_i}-postDur*track_Fs:minBehavFlightLength{clust_i}), activity_allTrials.vel_dur{clust_i},'UniformOutput',false); %cut the traces so they are aligned to the shortest flight trajectory across all the days
    %catch
    %end
    
    %mean and std of the trials as a whole for correlations in future
    mean_act_aligned{clust_i} = cellfun(@nanmean,act_aligned{clust_i},'UniformOutput',false); %take the mean of the aligned traces
    std_act_aligned{clust_i} = cellfun(@nanstd,act_aligned{clust_i},'UniformOutput',false); %take the mean of the aligned traces
    mean_act_pre{clust_i} = cellfun(@nanmean,act_pre{clust_i},'UniformOutput',false); %take the mean of the aligned traces
    std_act_pre{clust_i} = cellfun(@nanstd,act_pre{clust_i},'UniformOutput',false); %take the mean of the aligned traces
    mean_act_dur{clust_i} = cellfun(@nanmean,act_dur{clust_i},'UniformOutput',false); %take the mean of the aligned traces
    std_act_dur{clust_i} = cellfun(@nanstd,act_dur{clust_i},'UniformOutput',false); %take the mean of the aligned traces
    mean_act_post{clust_i} = cellfun(@nanmean,act_post{clust_i},'UniformOutput',false); %take the mean of the aligned traces
    std_act_post{clust_i} = cellfun(@nanstd,act_post{clust_i},'UniformOutput',false); %take the mean of the aligned traces
    mean_vel_aligned{clust_i} = cellfun(@nanmean,vel_aligned{clust_i},'UniformOutput',false); %take the mean of the aligned traces
    std_vel_aligned{clust_i} = cellfun(@nanstd,vel_aligned{clust_i},'UniformOutput',false); %take the mean of the aligned traces
    mean_vel_pre{clust_i} = cellfun(@nanmean,vel_pre{clust_i},'UniformOutput',false); %take the mean of the aligned traces
    std_vel_pre{clust_i} = cellfun(@nanstd,vel_pre{clust_i},'UniformOutput',false); %take the mean of the aligned traces
    mean_vel_dur{clust_i} = cellfun(@nanmean,vel_dur{clust_i},'UniformOutput',false); %take the mean of the aligned traces
    std_vel_dur{clust_i} = cellfun(@nanstd,vel_dur{clust_i},'UniformOutput',false); %take the mean of the aligned traces
    mean_vel_post{clust_i} = cellfun(@nanmean,vel_post{clust_i},'UniformOutput',false); %take the mean of the aligned traces
    std_vel_post{clust_i} = cellfun(@nanstd,vel_post{clust_i},'UniformOutput',false); %take the mean of the aligned traces
    
    
    %initialize matrices to hold mean and stds of each day/roi
    meanMax_act_aligned{clust_i} = zeros(size(days,2),size(act_aligned{clust_i},2));
    stdMax_act_aligned{clust_i} = zeros(size(days,2),size(act_aligned{clust_i},2));
    meanMax_act_pre{clust_i} = zeros(size(days,2),size(act_aligned{clust_i},2));
    stdMax_act_pre{clust_i} = zeros(size(days,2),size(act_aligned{clust_i},2));
    meanMax_act_dur{clust_i} = zeros(size(days,2),size(act_aligned{clust_i},2));
    stdMax_act_dur{clust_i} = zeros(size(days,2),size(act_aligned{clust_i},2));
    meanMax_act_post{clust_i} = zeros(size(days,2),size(act_aligned{clust_i},2));
    stdMax_act_post{clust_i} = zeros(size(days,2),size(act_aligned{clust_i},2));
    ylab{clust_i} = cell(1,3);
    overlapSorted{clust_i} = cell(size(days,2),3);
    intersectDays{clust_i} = cell(1,3);
    
    %determine the max,mean, and std for the pre, dur, post data
    for day_i = 1:length(days)
        for roi_i = 1:size(act_aligned{clust_i},2)
            for flight_i = 1:size(act_aligned{clust_i}{days(day_i),roi_i},1)
                [~,max_act_aligned{clust_i}{day_i,roi_i}(flight_i,1)] = max(act_aligned{clust_i}{days(day_i),roi_i}(flight_i,:));
                [~,max_act_pre{clust_i}{day_i,roi_i}(flight_i,1)] = max(act_pre{clust_i}{days(day_i),roi_i}(flight_i,:));
                [~,max_act_dur{clust_i}{day_i,roi_i}(flight_i,1)] = max(act_dur{clust_i}{days(day_i),roi_i}(flight_i,:));
                [~,max_act_post{clust_i}{day_i,roi_i}(flight_i,1)] = max(act_post{clust_i}{days(day_i),roi_i}(flight_i,:));
            end
            meanMax_act_aligned{clust_i}(day_i,roi_i) = round(mean(max_act_aligned{clust_i}{day_i,roi_i}));
            stdMax_act_aligned{clust_i}(day_i,roi_i) = round(std(max_act_aligned{clust_i}{day_i,roi_i}));
            meanMax_act_pre{clust_i}(day_i,roi_i) = round(mean(max_act_pre{clust_i}{day_i,roi_i}));
            stdMax_act_pre{clust_i}(day_i,roi_i) = round(std(max_act_pre{clust_i}{day_i,roi_i}));
            meanMax_act_dur{clust_i}(day_i,roi_i) = round(mean(max_act_dur{clust_i}{day_i,roi_i}));
            stdMax_act_dur{clust_i}(day_i,roi_i) = round(std(max_act_dur{clust_i}{day_i,roi_i}));
            meanMax_act_post{clust_i}(day_i,roi_i) = round(mean(max_act_post{clust_i}{day_i,roi_i}));
            stdMax_act_post{clust_i}(day_i,roi_i) = round(std(max_act_post{clust_i}{day_i,roi_i}));
        end
        %sort the means of maxes for each full, pre, dur, and post
        [BmeanMax_act_aligned{clust_i}(day_i,:),ImeanMax_act_aligned{clust_i}(day_i,:)] = sort(meanMax_act_aligned{clust_i}(day_i,:));
        [BmeanMax_act_pre{clust_i}(day_i,:),ImeanMax_act_pre{clust_i}(day_i,:)] = sort(meanMax_act_pre{clust_i}(day_i,:));
        [BmeanMax_act_dur{clust_i}(day_i,:),ImeanMax_act_dur{clust_i}(day_i,:)] = sort(meanMax_act_dur{clust_i}(day_i,:));
        [BmeanMax_act_post{clust_i}(day_i,:),ImeanMax_act_post{clust_i}(day_i,:)] = sort(meanMax_act_post{clust_i}(day_i,:));
        
        %load the stable place cell data for each day
        cd(dirDates(days(day_i)).name);
        dirExtracted = dir('*Extracted_trajectories*');
        load(dirExtracted(end).name);
        cd(dirDates(days(day_i)).folder);
        %pullout the selectively active cells
        selectiveCells{day_i,1} = placeCellsAngStable.ppre_cells;
        selectiveCells{day_i,2} = placeCellsAngStable.pp_cells;
        selectiveCells{day_i,3} = placeCellsAngStable.ppost_cells;
    end
end
%sort the means of the maxes so that they are in chronological order
%across ROIs
for day_i = 1:length(days)
    for clust_ii = 1:size(selectiveCells{day_i,1},1)
        for phase_i = 1:size(selectiveCells,2)
            findSelectiveCells{clust_ii}{day_i,phase_i} = find(selectiveCells{day_i,phase_i}(clust_ii,:));
            if phase_i == 1
                [BSelectiveCells{clust_ii}{day_i,phase_i},ISelectiveCells{clust_ii}{day_i,phase_i}] = sort(meanMax_act_pre{clust_ii}(day_i,findSelectiveCells{clust_ii}{day_i,phase_i}));
            elseif phase_i == 2
                [BSelectiveCells{clust_ii}{day_i,phase_i},ISelectiveCells{clust_ii}{day_i,phase_i}] = sort(meanMax_act_dur{clust_ii}(day_i,findSelectiveCells{clust_ii}{day_i,phase_i}));
            elseif phase_i == 3
                [BSelectiveCells{clust_ii}{day_i,phase_i},ISelectiveCells{clust_ii}{day_i,phase_i}] = sort(meanMax_act_post{clust_ii}(day_i,findSelectiveCells{clust_ii}{day_i,phase_i}));
            end
        end
    end
end

%find ROIs that intersect across all days
for clust_ii = 1:size(clustList,1)
    for day_i = 1:size(selectiveCells,1)
        for phase_i = 1:size(selectiveCells,2)
            try %make a list of the intersecting cells so it can be concatenated and sorted later
                intersectList{clust_ii,phase_i}{day_i} = ISelectiveCells{clustList(clust_ii,day_i)}{day_i,phase_i};
            catch
            end
        end
    end
end
for clust_ii = 1:size(clustList,1)
    for day_i = 1:size(selectiveCells,1)
        for phase_i = 1:size(selectiveCells,2)
            intersectDays{clust_ii+1}{phase_i} = mintersect(intersectList{clust_ii,phase_i}{:});
            roiNum = 1; %counter
            for roi_i = 1:length(ISelectiveCells{clust_ii+1}{day_i,phase_i})
                if ismember(ISelectiveCells{clust_ii+1}{day_i,phase_i}(roi_i),intersectDays{clust_ii+1}{phase_i})
                    overlapSorted{clust_ii+1}{day_i,phase_i}(roiNum) = ISelectiveCells{clust_ii+1}{day_i,phase_i}(roi_i); %sort each intersecting ROI based on the sorting from the original data
                    if day_i ==1
                        ylab{clust_ii+1}{phase_i} = [ylab{clust_ii+1}{phase_i} overlapSorted{clust_ii+1}{1,phase_i}(roiNum)];%make label for the yaxis for the plots
                    end
                    roiNum = roiNum +1;
                end
            end
        end
    end
end
% for day_i = 1:size(selectiveCells,1)
%     for clust_ii = 1:size(selectiveCells{day_i,1},1)
%         for phase_i = 1:size(selectiveCells,2)
%             intersectDays2{clust_ii}{phase_i} = mintersect(ISelectiveCells{clust_ii}{:,phase_i}); %all intersecting ROIs across all days
%             roiNum = 1; %counter
%             for roi_i = 1:length(ISelectiveCells{clust_ii}{day_i,phase_i})
%                 if ismember(ISelectiveCells{clust_ii}{day_i,phase_i}(roi_i),intersectDays2{clust_ii}{phase_i})
%                     overlapSorted2{clust_ii}{day_i,phase_i}(roiNum) = ISelectiveCells{clust_ii}{day_i,phase_i}(roi_i); %sort each intersecting ROI based on the sorting from the original data
%                     if day_i ==1
%                         ylab2{clust_ii}{phase_i} = [ylab{clust_ii}{phase_i} overlapSorted2{clust_ii}{1,phase_i}(roiNum)];%make label for the yaxis for the plots
%                     end
%                     roiNum = roiNum +1;
%                 end
%             end
%         end
%     end
% end



dataPreDurPost.batId = batId;
dataPreDurPost.nDays = nDays;
dataPreDurPost.day1 = day1;
dataPreDurPost.dayEnd = dayEnd;
dataPreDurPost.days = days;
dataPreDurPost.clusts = clusts;
dataPreDurPost.preDur = preDur;
dataPreDurPost.postDur = postDur;
dataPreDurPost.CNMFe_Fs = CNMFe_Fs;
dataPreDurPost.track_Fs = track_Fs;
dataPreDurPost.minFlightLength = minFlightLength;
dataPreDurPost.flightLength = flightLength;
dataPreDurPost.act_aligned = act_aligned;
dataPreDurPost.act_pre = act_pre;
dataPreDurPost.act_dur = act_dur;
dataPreDurPost.act_post = act_post;
dataPreDurPost.mean_act_aligned = mean_act_aligned;
dataPreDurPost.std_act_aligned = std_act_aligned;
dataPreDurPost.mean_act_pre = mean_act_pre;
dataPreDurPost.std_act_pre = std_act_pre;
dataPreDurPost.mean_act_dur = mean_act_dur;
dataPreDurPost.std_act_dur = std_act_dur;
dataPreDurPost.mean_act_post = mean_act_post;
dataPreDurPost.std_act_post = std_act_post;

dataPreDurPost.minBehavFlightLength = minBehavFlightLength;
dataPreDurPost.behavLength = behavLength;
dataPreDurPost.vel_aligned = vel_aligned;
dataPreDurPost.vel_pre = vel_pre;
dataPreDurPost.vel_dur = vel_dur;
dataPreDurPost.vel_post = vel_post;
dataPreDurPost.mean_vel_aligned = mean_vel_aligned;
dataPreDurPost.std_vel_aligned = std_vel_aligned;
dataPreDurPost.mean_vel_pre = mean_vel_pre;
dataPreDurPost.std_vel_pre = std_vel_pre;
dataPreDurPost.mean_vel_dur = mean_vel_dur;
dataPreDurPost.std_vel_dur = std_vel_dur;
dataPreDurPost.mean_vel_post = mean_vel_post;
dataPreDurPost.std_vel_post = std_vel_post;

dataPreDurPost.meanMax_act_aligned = meanMax_act_aligned;
dataPreDurPost.stdMax_act_aligned = stdMax_act_aligned;
dataPreDurPost.meanMax_act_pre = meanMax_act_pre;
dataPreDurPost.stdMax_act_pre = stdMax_act_aligned;
dataPreDurPost.meanMax_act_dur = meanMax_act_dur;
dataPreDurPost.stdMax_act_dur = stdMax_act_aligned;
dataPreDurPost.meanMax_act_post = meanMax_act_post;
dataPreDurPost.stdMax_act_post = stdMax_act_aligned;
dataPreDurPost.ylab = ylab;
dataPreDurPost.overlapSorted = overlapSorted;
dataPreDurPost.intersectDays = intersectDays;
dataPreDurPost.BmeanMax_act_aligned = BmeanMax_act_aligned;
dataPreDurPost.ImeanMax_act_aligned = ImeanMax_act_aligned;
dataPreDurPost.BmeanMax_act_pre = BmeanMax_act_pre;
dataPreDurPost.ImeanMax_act_pre = ImeanMax_act_pre;
dataPreDurPost.BmeanMax_act_dur = BmeanMax_act_dur;
dataPreDurPost.ImeanMax_act_dur = ImeanMax_act_dur;
dataPreDurPost.BmeanMax_act_post = BmeanMax_act_post;
dataPreDurPost.ImeanMax_act_post = ImeanMax_act_post;
dataPreDurPost.selectiveCells = selectiveCells;
dataPreDurPost.findSelectiveCells = findSelectiveCells;
dataPreDurPost.BSelectiveCells = BSelectiveCells;
dataPreDurPost.ISelectiveCells = ISelectiveCells;

if saveFlag == 1
    if strcmp(batId,'Gal')
    save([pwd filesep '200311to200320_Gal_dataPreDurPost_' saveTag '.mat'],'dataPreDurPost');
    elseif strcmp(batId,'Gen')
      save([pwd filesep '200319to200324_Gen_dataPreDurPost_' saveTag '.mat'],'dataPreDurPost');
  end
end

