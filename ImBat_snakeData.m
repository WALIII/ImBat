function [snakeTrace] = ImBat_snakeData(cellData,flightPaths,alignment,varargin)

batName = [];
dateSesh = [];
sessionType = [];
loadFlag = 0; %do you want to load and save the data individually outside of ImBatAnalyze
saveFlag = 0;
%offset = 0.1; % account for slow calcium estimation ~move locations back 100ms in time... This is the knob to turn for 'prospective' coding...
nClusters = 5; %number of flight trajectories to look at and for k means clustering of whole time series by peak

% User inputs overrides
nparams=length(varargin);
for i=1:2:nparams
    switch lower(varargin{i})
        case 'batname'
            batName=varargin{i+1};
        case 'datesesh'
            dateSesh = varargin{i+1};
        case 'sessiontype'
            sessionType = varargin{i+1};
        case 'loadflag'
            loadFlag = varargin{i+1};
        case 'saveflag'
            saveFlag = varargin{i+1};
    end
end

%labels for loading and saving data if running independent fromImBat_analyze
if loadFlag == 1
    date = strcat(lower(batName(1:2)),dateSesh);
    label = [batName '_' dateSesh '_' sessionType];
    
    cellData = load([pwd '/processed/Motion_corrected_Data_DS_results.mat']);
    alignment = load([pwd '/processed/Alignment.mat']);
    load([pwd '/analysis/' label '_flightPaths_6clusters.mat']);
end
%padding for during flight snake plot to include some time before and after flight
prePad = 2; %number of seconds to plot before alignment point
postPad = 6; %number of seconds to plot after alignment point
prePadCalcium = prePad*cellData.results.Fs; %number of frames (seconds*freq) to include in the trace extraction
postPadCalcium = postPad*cellData.results.Fs; %add 2 seconds to the end of the plots to include delay in peak time
prePadSpeed = prePad*120; %add 2 seconds * FS of tracking data (120)
postPadSpeed = postPad*120;%add 6 seconds * FS of tracking data (120)
%padding for pre and post flight snakePlots
preFlightPad = 10; %number of seconds to include before flight starts
postFlightPad = 10; %of of seconds to include after flight ends
preFlightPadCalcium = preFlightPad*cellData.results.Fs; %number of frames (seconds*freq) to include in the trace extraction
postFlightPadCalcium = postFlightPad*cellData.results.Fs; %add 2 seconds to the end of the plots to include delay in peak time
preFlightPadSpeed = preFlightPad*120; %add 2 seconds * FS of tracking data (120)
postFlightPadSpeed = postFlightPad*120;%add 6 seconds * FS of tracking data (120)

%meanTrace = cell(1,length(flightPaths.clusterIndex));
meanTraceFlightAll = []; %initialize the variable to concatenate all traces for uniform zscore
meanTracePreAll = []; %initialize the variable to concatenate all traces for uniform zscore
meanTracePostAll = []; %initialize the variable to concatenate all traces for uniform zscore
meanTracePreFlightPostAll = []; %initialize the variable to concatenate all pre/post/flight traces for uniform zscore


%this is to merge specific clusters if 2 or more of them are the same but
%were separated by the flight k-means
flightPaths.clusterIndex{2}=cat(2,flightPaths.clusterIndex{2},flightPaths.clusterIndex{3});
flightPaths.clusterIndex{3} =[];
flightPaths.clusterIndex= flightPaths.clusterIndex(~cellfun('isempty',flightPaths.clusterIndex));

%for each cluster type
for clust_i = 1:nClusters %length(flightPaths.clusterIndex)
    %for each cell
    for cell_i = 1:length(cellData.results.C(:,1))
        %for each trial in cluster, calculate the start/stop times and duration of the calcium videos
        for dur_i = 1:length(flightPaths.clusterIndex{clust_i})
            %get imaging times of start and stop index converting from tracking to video times
            [minValueStart(dur_i),closestIndexStart(dur_i)] = min(abs(alignment.out.video_times-alignment.out.Location_time(flightPaths.flight_starts_idx(flightPaths.clusterIndex{clust_i}(dur_i)))));
            [minValueEnd(dur_i),closestIndexEnd(dur_i)] = min(abs(alignment.out.video_times-alignment.out.Location_time(flightPaths.flight_ends_idx(flightPaths.clusterIndex{clust_i}(dur_i)))));
            %calculate duration of each flight in a particular cluster so
            %you can pad all flights to the longest flight in that cluster
            dur(dur_i) = closestIndexEnd(dur_i)-closestIndexStart(dur_i);
            durSpeed(dur_i)= flightPaths.flight_ends_idx(flightPaths.clusterIndex{clust_i}(dur_i))-flightPaths.flight_starts_idx(flightPaths.clusterIndex{clust_i}(dur_i));
        end
        %calculate max duration for each cluster of trajectories
        maxDur = max(dur);
        maxDurSpeed = max(durSpeed);
        %meanTrace{clust_i}=zeros(length(cellData.results.C(:,1)),maxDur+preWindow+1);
        %initialize the vector to store the neural activity of each flight
        traceFlight = zeros(length(flightPaths.clusterIndex{clust_i}),maxDur+prePadCalcium+postPadCalcium+1);
        speedFlight{clust_i} = zeros(length(flightPaths.clusterIndex{clust_i}),maxDurSpeed+prePadSpeed+postPadSpeed+1);
        tracePre = zeros(length(flightPaths.clusterIndex{clust_i}),preFlightPadCalcium+1);
        speedPre{clust_i} = zeros(length(flightPaths.clusterIndex{clust_i}),preFlightPadSpeed+1);
        tracePost = zeros(length(flightPaths.clusterIndex{clust_i}),postFlightPadCalcium+1);
        speedPost{clust_i} = zeros(length(flightPaths.clusterIndex{clust_i}),postFlightPadSpeed+1);
        %build the calcium and speed vectors for each flight within each cluster
        for trace_i = 1:length(flightPaths.clusterIndex{clust_i})
            try
                traceFlight(trace_i,:) = cellData.results.C(cell_i,closestIndexStart(trace_i) - prePadCalcium:closestIndexEnd(trace_i) + (maxDur-dur(trace_i)) + postPadCalcium);
                speedFlight{clust_i}(trace_i,:) = flightPaths.batSpeed(flightPaths.flight_starts_idx(flightPaths.clusterIndex{clust_i}(trace_i)) - prePadSpeed:flightPaths.flight_ends_idx(flightPaths.clusterIndex{clust_i}(trace_i)) + (maxDurSpeed-durSpeed(trace_i)) + postPadSpeed);
                smoothSpeedRawFlight{clust_i}(trace_i,:) = smooth(speedFlight{clust_i}(trace_i,:),100);
                tracePre(trace_i,:) = cellData.results.C(cell_i,closestIndexStart(trace_i) - preFlightPadCalcium:closestIndexStart(trace_i));
                speedPre{clust_i}(trace_i,:) = flightPaths.batSpeed(flightPaths.flight_starts_idx(flightPaths.clusterIndex{clust_i}(trace_i)) - preFlightPadSpeed:flightPaths.flight_starts_idx(flightPaths.clusterIndex{clust_i}(trace_i)));
                smoothSpeedRawPre{clust_i}(trace_i,:) = smooth(speedPre{clust_i}(trace_i,:),100);
                tracePost(trace_i,:) = cellData.results.C(cell_i,closestIndexEnd(trace_i):closestIndexEnd(trace_i) + postFlightPadCalcium);
                speedPost{clust_i}(trace_i,:) = flightPaths.batSpeed(flightPaths.flight_ends_idx(flightPaths.clusterIndex{clust_i}(trace_i)):flightPaths.flight_ends_idx(flightPaths.clusterIndex{clust_i}(trace_i)) + postFlightPadSpeed);
                smoothSpeedRawPost{clust_i}(trace_i,:) = smooth(speedPost{clust_i}(trace_i,:),100);
            catch
                sizeToRecordingEnd = size(cellData.results.C(cell_i,closestIndexStart(trace_i) - preFlightPadCalcium:end),2);
                sizeToTraceEnd = size(traceFlight(trace_i,:),2);
                traceFlight(trace_i,:) = (cellData.results.C(cell_i,closestIndexStart(trace_i) - preFlightPadCalcium:end + postFlightPadCalcium)+(zeros(1,sizeToTraceEnd - sizeToRecordingEnd)));
                speedFlight{clust_i}(trace_i,:) = (flightPaths.batSpeed(closestIndexStart(trace_i) - preFlightPadSpeed:end + postFlightPadSpeed)+(zeros(1,sizeToTraceEnd - sizeToRecordingEnd)));
                disp('End of rec')
            end
        end
        
        %calculate the mean neural activity across all flights in a cluster for each cell
        meanTraceFlight{clust_i}(cell_i,:) = mean(traceFlight);
        meanTraceOdd{clust_i}(cell_i,:) = mean(traceFlight(1:2:end,:));
        meanTraceEven{clust_i}(cell_i,:) = mean(traceFlight(2:2:end,:));
        meanTracePre{clust_i}(cell_i,:) = mean(tracePre);
        meanTracePost{clust_i}(cell_i,:) = mean(tracePost);
        meanSpeedFlight{clust_i} = mean(speedFlight{clust_i});
        meanSpeedPre{clust_i} = mean(speedPre{clust_i});
        meanSpeedPost{clust_i} = mean(speedPost{clust_i});

        
        %smooth and zscore the neural data. subtract the min of the zscore so the
        %min is 0 rather than mean 0
        normTraceFlight{clust_i}(cell_i,:) = zscore(smooth(meanTraceFlight{clust_i}(cell_i,:),3));
        normTraceFlight{clust_i}(cell_i,:) = normTraceFlight{clust_i}(cell_i,:) - min(normTraceFlight{clust_i}(cell_i,:));
        smoothSpeedFlight{clust_i} = smooth(meanSpeedFlight{clust_i},40);
        normTraceOdd{clust_i}(cell_i,:) = zscore(smooth(meanTraceOdd{clust_i}(cell_i,:),3));
        normTraceOdd{clust_i}(cell_i,:) = normTraceOdd{clust_i}(cell_i,:) - min(normTraceOdd{clust_i}(cell_i,:));
        normTraceEven{clust_i}(cell_i,:) = zscore(smooth(meanTraceEven{clust_i}(cell_i,:),3));
        normTraceEven{clust_i}(cell_i,:) = normTraceEven{clust_i}(cell_i,:) - min(normTraceEven{clust_i}(cell_i,:));
        normTracePre{clust_i}(cell_i,:) = zscore(smooth(meanTracePre{clust_i}(cell_i,:),3));
        normTracePre{clust_i}(cell_i,:) = normTracePre{clust_i}(cell_i,:) - min(normTracePre{clust_i}(cell_i,:));       
        smoothSpeedPre{clust_i} = smooth(meanSpeedPre{clust_i},40);
        normTracePost{clust_i}(cell_i,:) = zscore(smooth(meanTracePost{clust_i}(cell_i,:),3));
        normTracePost{clust_i}(cell_i,:) = normTracePost{clust_i}(cell_i,:) - min(normTracePost{clust_i}(cell_i,:));       
        smoothSpeedPost{clust_i} = smooth(meanSpeedPost{clust_i},40);
        %find time index of max peaks
        [~,maxnormTraceFlight{clust_i}(cell_i,1)] = max(normTraceFlight{clust_i}(cell_i,:));
        [~,maxnormTracePre{clust_i}(cell_i,1)] = max(normTracePre{clust_i}(cell_i,:));
        [~,maxnormTracePost{clust_i}(cell_i,1)] = max(normTracePost{clust_i}(cell_i,:));
        [~,maxNormTraceOdd{clust_i}(cell_i,1)] = max(normTraceOdd{clust_i}(cell_i,:));
        
    end
    %sort each cell by the timing of its peak firing
    [BFlight{clust_i},IFlight{clust_i}] = sort(maxnormTraceFlight{clust_i});
    [BPre{clust_i},IPre{clust_i}] = sort(maxnormTracePre{clust_i});
    [BPost{clust_i},IPost{clust_i}] = sort(maxnormTracePost{clust_i});
    [B1Flight{clust_i},I1Flight{clust_i}] = sort(maxnormTraceFlight{1}); %sort by cluster 1
    [B1Pre{clust_i},I1Pre{clust_i}] = sort(maxnormTracePre{1}); %sort by cluster 1
    [B1Post{clust_i},I1Post{clust_i}] = sort(maxnormTracePost{1}); %sort by cluster 1    
    [Bodd{clust_i},Iodd{clust_i}] = sort(maxNormTraceOdd{clust_i});
    %split dataset into even and odd for comparing 2 halves
    %     normTraceOdd{clust_i} = normTraceFlight{clust_i}(1:2:end,:);
    %      normTraceEven{clust_i} = normTraceFlight{clust_i}(2:2:end,:);
    %     %sort odd clusters in ascending order of the peak of the odds
    %     if length(normTraceOdd{clust_i}(:,1))>length(normTraceEven{clust_i}(:,1)) %if number of odd elements is greater than even
    %     [Bodd{clust_i},Iodd{clust_i}] = sort(maxNormTrace{clust_i}(1:2:end-1));
    %     else
    %     [Bodd{clust_i},Iodd{clust_i}] = sort(maxnormTraceFlight{clust_i}(1:2:end));
    %     end
    
    %concatenate all clusters together and then zscore to equalize the
    %peaks across all trials
    meanTraceFlightAll = [meanTraceFlightAll meanTraceFlight{clust_i}]; %concatenate all traces
    meanTracePreAll = [meanTracePreAll meanTracePre{clust_i}]; %concatenate all traces
    meanTracePostAll = [meanTracePostAll meanTracePost{clust_i}]; %concatenate all traces
    meanTracePreFlightPostAll = [meanTracePreFlightPostAll meanTracePreFlightPost{clust_i} 
end
%% this is to smooth, zscore, and sort the entire cell data by their preferred flight according to a homemade k-means (max dff across flights)
%zscore the full data set and subtract min to start at 0
for cell_ii = 1:length(cellData.results.C(:,1))
    normMeanTraceFlightAll(cell_ii,:) = zscore(smooth(meanTraceFlightAll(cell_ii,:),10));
    normMeanTraceFlightAll(cell_ii,:) = normMeanTraceFlightAll(cell_ii,:) - min(normMeanTraceFlightAll(cell_ii,:));
    normMeanTracePreAll(cell_ii,:) = zscore(smooth(meanTracePreAll(cell_ii,:),10));
    normMeanTracePreAll(cell_ii,:) = normMeanTracePreAll(cell_ii,:) - min(normMeanTracePreAll(cell_ii,:));
    normMeanTracePostAll(cell_ii,:) = zscore(smooth(meanTracePostAll(cell_ii,:),10));
    normMeanTracePostAll(cell_ii,:) = normMeanTracePostAll(cell_ii,:) - min(normMeanTracePostAll(cell_ii,:));
end
%split the data back into its clusters
traceFlightIndConcat = [1 length(normTraceFlight{1}(1,:))]; %initialize index variable to hold the indices for start/stop of each cluster
tracePreIndConcat = [1 length(normTracePre{1}(1,:))]; %initialize index variable to hold the indices for start/stop of each cluster
tracePostIndConcat = [1 length(normTracePost{1}(1,:))]; %initialize index variable to hold the indices for start/stop of each cluster
%find the start and stop indices for 2nd through n cluster
for clust_ii = 2:clust_i
    traceFlightIndConcat = vertcat(traceFlightIndConcat,[traceFlightIndConcat(clust_ii-1,2)+1 traceFlightIndConcat(clust_ii-1,2)+length(normTraceFlight{clust_ii}(1,:))]);
    tracePreIndConcat = vertcat(tracePreIndConcat,[tracePreIndConcat(clust_ii-1,2)+1 tracePreIndConcat(clust_ii-1,2)+length(normTracePre{clust_ii}(1,:))]);
    tracePostIndConcat = vertcat(tracePostIndConcat,[tracePostIndConcat(clust_ii-1,2)+1 tracePostIndConcat(clust_ii-1,2)+length(normTracePost{clust_ii}(1,:))]);

end
%regroup by each flight cluster
for clust_iii = 1:clust_i
    normTraceFlightAll{clust_iii} = normMeanTraceFlightAll(:,traceFlightIndConcat(clust_iii,1):traceFlightIndConcat(clust_iii,2));
    normTracePreAll{clust_iii} = normMeanTracePreAll(:,tracePreIndConcat(clust_iii,1):tracePreIndConcat(clust_iii,2));
    normTracePostAll{clust_iii} = normMeanTracePostAll(:,tracePostIndConcat(clust_iii,1):tracePostIndConcat(clust_iii,2));
    %find the order of the maximum for each cluster within the regrouping
    for cell_iii = 1:length(cellData.results.C(:,1))
        [~,maxNormFlightAll{clust_iii}(cell_iii,1)] = max(normTraceFlightAll{clust_iii}(cell_iii,:));
        [~,maxNormPreAll{clust_iii}(cell_iii,1)] = max(normTracePreAll{clust_iii}(cell_iii,:));
        [~,maxNormPostAll{clust_iii}(cell_iii,1)] = max(normTracePostAll{clust_iii}(cell_iii,:));
    end
    %sort each cell by the timing of its peak firing
    [BnormFlightAll{clust_iii},InormFlightAll{clust_iii}] = sort(maxNormFlightAll{clust_iii});
    [B1normFlightAll{clust_iii},I1normFlightAll{clust_iii}] = sort(maxNormFlightAll{1}); %sort by cluster 1
    [BnormPreAll{clust_iii},InormPreAll{clust_iii}] = sort(maxNormPreAll{clust_iii});
    [B1normPreAll{clust_iii},I1normPreAll{clust_iii}] = sort(maxNormPreAll{1}); %sort by cluster 1
    [BnormPosttAll{clust_iii},InormPostAll{clust_iii}] = sort(maxNormPostAll{clust_iii});
    [B1normPosttAll{clust_iii},I1normPostAll{clust_iii}] = sort(maxNormPostAll{1}); %sort by cluster 1
end


%find max and sort based on the peaks of the cell across the whole
%timeseries
[~,maxAllFlight] = max(normMeanTraceFlightAll,[],2);
[~,maxAllPre] = max(normMeanTracePreAll,[],2);
[~,maxAllPost] = max(normMeanTracePostAll,[],2);

%kPeaks = kmeans(maxAll,nClusters); %take k-means cluster of the peaks of each cell whole time
%sort based on each clustered peaks for flight, pre and post
rng(2);
for n = 1:size(maxAllFlight,1);
    for ii = 1:nClusters
        if maxAllFlight(n) >= traceFlightIndConcat(ii,1) & maxAllFlight(n) < traceFlightIndConcat(ii,2); %if the peak is within each cluster, sort it into that particular cluster
            kPeaksFlight(:,n) = 1;
        end
    end
end
for n = 1:size(maxAllPre,1);
    for ii = 1:nClusters
        if maxAllPre(n) >= tracePreIndConcat(ii,1) & maxAllPre(n) < tracePreIndConcat(ii,2); %if the peak is within each cluster, sort it into that particular cluster
            kPeaksPre(:,n) = 1;
        end
    end
end
for n = 1:size(maxAllPost,1);
    for ii = 1:nClusters
        if maxAllPost(n) >= tracePostIndConcat(ii,1) & maxAllPost(n) < tracePostIndConcat(ii,2); %if the peak is within each cluster, sort it into that particular cluster
            kPeaksPost(:,n) = 1;
        end
    end
end

%sort the clusters according to their peak times
[BkPeaksFlight,IkPeaksFlight] = sort(kPeaksFlight);
normMeanTraceSortFlight = normMeanTraceFlightAll(IkPeaksFlight,:);
[BkPeaksPre,IkPeaksPre] = sort(kPeaksPre);
normMeanTraceSortPre = normMeanTracePreAll(IkPeaksPre,:);
[BkPeaksPost,IkPeaksPost] = sort(kPeaksPost);
normMeanTraceSortPost = normMeanTracePostAll(IkPeaksPost,:);
%for each cluster, find the peaks of the cluster (aMat aka aTemp), tke the max
%(maxInd), and sort that max (aSort) from the cluster, add this to the
%normMeanTraceSort variable
for c = 1:nClusters
    aMatFlight = find(BkPeaksFlight == c);
    aTempFlight = normMeanTraceSortFlight(aMatFlight,:);
    [~,maxIndFlight] = max(aTempFlight,[],2);
    [BmaxFlight,ImaxFlight] = sort(maxIndFlight);
    aSortFlight = aTempFlight(ImaxFlight,:);
    normMeanTraceSortFlight(aMatFlight,:) = aSortFlight;
    clear aMatFlight aTempFlight maxIndFlight ImaxFlight BaxFlight aSortFlight
    aMatPre = find(BkPeaksPre == c);
    aTempPre = normMeanTraceSortPre(aMatPre,:);
    [~,maxIndPre] = max(aTempPre,[],2);
    [BmaxPre,ImaxPre] = sort(maxIndPre);
    aSortPre = aTempPre(ImaxPre,:);
    normMeanTraceSortPre(aMatPre,:) = aSortPre;
    clear aMatPre aTempPre maxIndPre ImaxPre BaxPre aSortPre
    aMatPost = find(BkPeaksPost == c);
    aTempPost = normMeanTraceSortPost(aMatPost,:);
    [~,maxIndPost] = max(aTempPost,[],2);
    [BmaxPost,ImaxPost] = sort(maxIndPost);
    aSortPost = aTempPost(ImaxPost,:);
    normMeanTraceSortPost(aMatPost,:) = aSortPost;
    clear aMatPost aTempPost maxIndPost ImaxPost BaxPost aSortPost
end
%divide the normMeanTraceSort overall variable into the respective clusters
for c = 1:nClusters
    normMeanTraceEachFlight{c} = normMeanTraceSortFlight(:,traceFlightIndConcat(c,1):traceFlightIndConcat(c,2));
    normMeanTraceEachPre{c} = normMeanTraceSortPre(:,tracePreIndConcat(c,1):tracePreIndConcat(c,2));
    normMeanTraceEachPost{c} = normMeanTraceSortPost(:,tracePostIndConcat(c,1):tracePostIndConcat(c,2));
end

%% save to snakeTrace variable
snakeTrace.meanTraceFlight = meanTraceFlight;
snakeTrace.normTraceFlight = normTraceFlight;
snakeTrace.maxnormTraceFlight = maxnormTraceFlight;
snakeTrace.meanTracePre = meanTracePre;
snakeTrace.normTracePre = normTracePre;
snakeTrace.maxnormTracePre = maxnormTracePre;
snakeTrace.meanTracePost = meanTracePost;
snakeTrace.normTracePost = normTracePost;
snakeTrace.maxnormTracePost = maxnormTracePost;
snakeTrace.BFlight = BFlight;
snakeTrace.IFlight = IFlight;
snakeTrace.BPre = BPre;
snakeTrace.IPre = IPre;
snakeTrace.BPost = BPost;
snakeTrace.IPost = IPost;
snakeTrace.BOdd = Bodd;
snakeTrace.IOdd = Iodd;
snakeTrace.B1Flight = B1Flight;
snakeTrace.I1Flight = I1Flight;
snakeTrace.B1Pre = B1Pre;
snakeTrace.I1Pre = I1Pre;
snakeTrace.B1Post = B1Post;
snakeTrace.I1Post = I1Post;
snakeTrace.Iodd = Iodd;
snakeTrace.Bodd = Bodd;
snakeTrace.normTraceEven = normTraceEven;
snakeTrace.normTraceOdd = normTraceOdd;
snakeTrace.meanSpeedFlight = meanSpeedFlight;
snakeTrace.smoothSpeedFlight = smoothSpeedFlight;
snakeTrace.meanSpeedPre = meanSpeedPre;
snakeTrace.smoothSpeedPre = smoothSpeedPre;
snakeTrace.meanSpeedPost = meanSpeedPost;
snakeTrace.smoothSpeedPost = smoothSpeedPost;
snakeTrace.smoothSpeedRawPre = smoothSpeedRawPre;
snakeTrace.smoothSpeedRawPost = smoothSpeedRawPost;
snakeTrace.smoothSpeedRawFlight = smoothSpeedRawFlight;
snakeTrace.normMeanTraceEachFlight = normMeanTraceEachFlight;
snakeTrace.normMeanTraceSortFlight = normMeanTraceSortFlight;
snakeTrace.normMeanTraceFlightAll = normMeanTraceFlightAll;
snakeTrace.traceFlightIndConcat = traceFlightIndConcat;
snakeTrace.normMeanTraceEachPre = normMeanTraceEachPre;
snakeTrace.normMeanTraceSortPre = normMeanTraceSortPre;
snakeTrace.normMeanTracePreAll = normMeanTracePreAll;
snakeTrace.tracePreIndConcat = tracePreIndConcat;
snakeTrace.normMeanTraceEachPost = normMeanTraceEachPost;
snakeTrace.normMeanTraceSortPost = normMeanTraceSortPost;
snakeTrace.normMeanTracePostAll = normMeanTracePostAll;
snakeTrace.tracePostIndConcat = tracePostIndConcat;
snakeTrace.nClusters = nClusters;
snakeTrace.batName = batName;
snakeTrace.dateSesh = dateSesh;
snakeTrace.sessionType = sessionType;
snakeTrace.InormFlightAll = InormFlightAll;
snakeTrace.normTraceFlightAll = normTraceFlightAll;
snakeTrace.I1normFlightAll = I1normFlightAll;

if saveFlag == 1
    snakeTrace.label = label;
   save([pwd '/analysis/' label '_snakePlotData.mat'],'snakeTrace');
end

%snakeTrace.normMeanTraceAllSmooth = normMeanTraceAllSmooth;

end