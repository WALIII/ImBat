function [dataPreDurPost] = ImBat_extract_dataPreDurPost
saveFlag = 1;
saveTag = 'cRaw_vel';
cRaw = 1;
%load first day placeCellStableROI data
%cd([dirTop(day1).folder filesep 'plots\200911-preDurPost cells across days']);
if cRaw ==1
    load('200311to200320_Gal_activity_allTrials_allClusts_cRaw_vel.mat'); %load activity for pre,dur,post
else
    load('200311to200320_Gal_activity_allTrials_allClusts_sMat.mat'); %load s matrix activity data
end
dirGalDates = dir('Gal_*');
nDays = length(dirGalDates);
day1 = 1; %first day to start comparing
dayEnd = nDays; %last day to start comparing
days = [1:9]; %all days concat
clusts = 2;
CNMFe_Fs = 30;
track_Fs = 120;
preDur = 3; %pre takeoff
postDur = 7; %post takeoff

%pull in/cutout all the aligned, pre, dur, and post data
minFlightLength = min(cellfun('size',activity_allTrials.act_dur{clusts}(:,1),2)); %shortest flight length for each day since they may be slightly different lengths
minBehavFlightLength = min(cellfun('size',activity_allTrials.vel_dur{clusts}(1,:),2)); %for behavior
flightLength = minFlightLength - preDur*CNMFe_Fs - postDur*CNMFe_Fs;
behavLength = minBehavFlightLength - preDur*track_Fs - postDur*track_Fs;
act_aligned = cellfun(@(c) c(:,1:minFlightLength), activity_allTrials.act_dur{clusts},'UniformOutput',false); %cut the traces so they are aligned to the shortest flight trajectory across all the days
act_pre = cellfun(@(c) c(:,1:preDur*CNMFe_Fs), activity_allTrials.act_dur{clusts},'UniformOutput',false); %cut the traces so they are aligned to the shortest flight trajectory across all the days
act_dur = cellfun(@(c) c(:,1+preDur*CNMFe_Fs:minFlightLength-postDur*CNMFe_Fs), activity_allTrials.act_dur{clusts},'UniformOutput',false); %cut the traces so they are aligned to the shortest flight trajectory across all the days
act_post = cellfun(@(c) c(:,1+ minFlightLength-postDur*CNMFe_Fs:minFlightLength), activity_allTrials.act_dur{clusts},'UniformOutput',false); %cut the traces so they are aligned to the shortest flight trajectory across all the days
vel_aligned = cellfun(@(c) c(:,1:minBehavFlightLength), activity_allTrials.vel_dur{clusts},'UniformOutput',false); %cut the traces so they are aligned to the shortest flight trajectory across all the days
vel_pre = cellfun(@(c) c(:,1:preDur*track_Fs), activity_allTrials.vel_dur{clusts},'UniformOutput',false); %cut the traces so they are aligned to the shortest flight trajectory across all the days
vel_dur = cellfun(@(c) c(:,1+preDur*track_Fs:minBehavFlightLength-postDur*track_Fs), activity_allTrials.vel_dur{clusts},'UniformOutput',false); %cut the traces so they are aligned to the shortest flight trajectory across all the days
vel_post = cellfun(@(c) c(:,1+minBehavFlightLength-postDur*track_Fs:minBehavFlightLength), activity_allTrials.vel_dur{clusts},'UniformOutput',false); %cut the traces so they are aligned to the shortest flight trajectory across all the days

%mean and std of the trials as a whole for correlations in future
mean_act_aligned = cellfun(@mean,act_aligned,'UniformOutput',false); %take the mean of the aligned traces
std_act_aligned = cellfun(@std,act_aligned,'UniformOutput',false); %take the mean of the aligned traces
mean_act_pre = cellfun(@mean,act_pre,'UniformOutput',false); %take the mean of the aligned traces
std_act_pre = cellfun(@std,act_pre,'UniformOutput',false); %take the mean of the aligned traces
mean_act_dur = cellfun(@mean,act_dur,'UniformOutput',false); %take the mean of the aligned traces
std_act_dur = cellfun(@std,act_dur,'UniformOutput',false); %take the mean of the aligned traces
mean_act_post = cellfun(@mean,act_post,'UniformOutput',false); %take the mean of the aligned traces
std_act_post = cellfun(@std,act_post,'UniformOutput',false); %take the mean of the aligned traces
mean_vel_aligned = cellfun(@mean,vel_aligned,'UniformOutput',false); %take the mean of the aligned traces
std_vel_aligned = cellfun(@std,vel_aligned,'UniformOutput',false); %take the mean of the aligned traces
mean_vel_pre = cellfun(@mean,vel_pre,'UniformOutput',false); %take the mean of the aligned traces
std_vel_pre = cellfun(@std,vel_pre,'UniformOutput',false); %take the mean of the aligned traces
mean_vel_dur = cellfun(@mean,vel_dur,'UniformOutput',false); %take the mean of the aligned traces
std_vel_dur = cellfun(@std,vel_dur,'UniformOutput',false); %take the mean of the aligned traces
mean_vel_post = cellfun(@mean,vel_post,'UniformOutput',false); %take the mean of the aligned traces
std_vel_post = cellfun(@std,vel_post,'UniformOutput',false); %take the mean of the aligned traces
 

%initialize matrices to hold mean and stds of each day/roi
meanMax_act_aligned = zeros(size(days,2),size(act_aligned,2));
stdMax_act_aligned = zeros(size(days,2),size(act_aligned,2));
meanMax_act_pre = zeros(size(days,2),size(act_aligned,2));
stdMax_act_pre = zeros(size(days,2),size(act_aligned,2));
meanMax_act_dur = zeros(size(days,2),size(act_aligned,2));
stdMax_act_dur = zeros(size(days,2),size(act_aligned,2));
meanMax_act_post = zeros(size(days,2),size(act_aligned,2));
stdMax_act_post = zeros(size(days,2),size(act_aligned,2));
ylab = cell(1,3);
overlapSorted = cell(size(days,2),3);
intersectDays = cell(1,3);

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
    
    %load the stable place cell data for each day
    cd(dirGalDates(days(day_i)).name);
    dirExtracted = dir('*Extracted_trajectories*');
    load(dirExtracted(end).name);
    cd(dirGalDates(days(day_i)).folder);
    %pullout the selectively active cells
    selectiveCells{day_i,1} = placeCellsAngStable.ppre_cells;
    selectiveCells{day_i,2} = placeCellsAngStable.pp_cells;
    selectiveCells{day_i,3} = placeCellsAngStable.ppost_cells;
    
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
end

%find ROIs that intersect across all days
for day_i = 1:size(selectiveCells,1)
    for phase_i = 1:size(selectiveCells,2)
        intersectDays{phase_i} = mintersect(ISelectiveCells{:,phase_i}); %all intersecting ROIs across all days
        roiNum = 1; %counter
        for roi_i = 1:length(ISelectiveCells{day_i,phase_i})
            if ismember(ISelectiveCells{day_i,phase_i}(roi_i),intersectDays{phase_i})
                overlapSorted{day_i,phase_i}(roiNum) = ISelectiveCells{day_i,phase_i}(roi_i); %sort each intersecting ROI based on the sorting from the original data
                if day_i ==1
                    ylab{phase_i} = [ylab{phase_i} overlapSorted{1,phase_i}(roiNum)];%make label for the yaxis for the plots
                end
                roiNum = roiNum +1;
            end
        end
    end
end

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
dataPreDurPost.findSelectiveCells = findSelectiveCells;
dataPreDurPost.BSelectiveCells = BSelectiveCells;
dataPreDurPost.ISelectiveCells = ISelectiveCells;

if saveFlag == 1
save([pwd filesep '200311to200320_Gal_dataPreDurPost_' saveTag '.mat'],'dataPreDurPost');
end

