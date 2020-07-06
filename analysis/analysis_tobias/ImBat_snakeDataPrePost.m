function [snakeTracePrePost] = ImBat_snakeDataPrePost(cellData,flightPaths,alignment,varargin)

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
prePad = 0; %number of seconds to plot before alignment point for during flight
postPad = 0; %number of seconds to plot after alignment point for during flight
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
meanTracePreFlightPostAll = []; %initialize the variable to concatenate all pre/post/flight traces for uniform zscore

%for each cell in each sort type (pre,during,post)

%for each cell
for cell_i = 1:length(cellData.results.C(:,1))
    %for each trial in cluster, calculate the start/stop times and duration of the calcium videos
    for dur_i = 1:length(flightPaths.flight_starts_idx)
        %get imaging times of start and stop index converting from tracking to video times
        [minValueStart(dur_i),closestIndexStart(dur_i)] = min(abs(alignment.out.video_times-alignment.out.Location_time(flightPaths.flight_starts_idx(dur_i))));
        [minValueEnd(dur_i),closestIndexEnd(dur_i)] = min(abs(alignment.out.video_times-alignment.out.Location_time(flightPaths.flight_ends_idx(dur_i))));
        %calculate duration of each flight in a particular cluster so
        %you can pad all flights to the longest flight in that cluster
        dur(dur_i) = closestIndexEnd(dur_i)-closestIndexStart(dur_i);
        durSpeed(dur_i)= flightPaths.flight_ends_idx(dur_i)-flightPaths.flight_starts_idx(dur_i);
    end
    %calculate max duration for each cluster of trajectories
    maxDur = max(dur);
    maxDurSpeed = max(durSpeed);
    
    %initialize the vector to store the neural activity of each flight
    traceFlight = zeros(length(flightPaths.flight_starts_idx),maxDur+prePadCalcium+postPadCalcium+1);
    speedFlight = zeros(length(flightPaths.flight_starts_idx),maxDurSpeed+prePadSpeed+postPadSpeed+1);
    tracePre = zeros(length(flightPaths.flight_starts_idx),preFlightPadCalcium+1);
    speedPre = zeros(length(flightPaths.flight_starts_idx),preFlightPadSpeed+1);
    tracePost = zeros(length(flightPaths.flight_starts_idx),postFlightPadCalcium+1);
    speedPost = zeros(length(flightPaths.flight_starts_idx),postFlightPadSpeed+1);
    
    %build the calcium and speed vectors for each flight within each cluster
        for trace_i = 1:length(flightPaths.flight_starts_idx)
            try
                traceFlight(trace_i,:) = cellData.results.C(cell_i,closestIndexStart(trace_i) - prePadCalcium:closestIndexEnd(trace_i) + (maxDur-dur(trace_i)) + postPadCalcium);
                speedFlight(trace_i,:) = flightPaths.batSpeed(flightPaths.flight_starts_idx(trace_i) - prePadSpeed:flightPaths.flight_ends_idx(trace_i) + (maxDurSpeed-durSpeed(trace_i)) + postPadSpeed);
                smoothSpeedRawFlight(trace_i,:) = smooth(speedFlight(trace_i,:),100);
                tracePre(trace_i,:) = cellData.results.C(cell_i,closestIndexStart(trace_i) - preFlightPadCalcium:closestIndexStart(trace_i));
                speedPre(trace_i,:) = flightPaths.batSpeed(flightPaths.flight_starts_idx(trace_i) - preFlightPadSpeed:flightPaths.flight_starts_idx(trace_i));
                smoothSpeedRawPre(trace_i,:) = smooth(speedPre(trace_i,:),100);
                tracePost(trace_i,:) = cellData.results.C(cell_i,closestIndexEnd(trace_i)+postPadCalcium:closestIndexEnd(trace_i)+postPadCalcium + postFlightPadCalcium);
                speedPost(trace_i,:) = flightPaths.batSpeed(flightPaths.flight_ends_idx(trace_i)+postPadCalcium:flightPaths.flight_ends_idx(trace_i)+postPadCalcium + postFlightPadSpeed);
                smoothSpeedRawPost(trace_i,:) = smooth(speedPost(trace_i,:),100);
            catch
                sizeToRecordingEnd = size(cellData.results.C(cell_i,closestIndexStart(trace_i) - preFlightPadCalcium:end),2);
                sizeToTraceEnd = size(traceFlight(trace_i,:),2);
                traceFlight(trace_i,:) = (cellData.results.C(cell_i,closestIndexStart(trace_i) - preFlightPadCalcium:end + postFlightPadCalcium)+(zeros(1,sizeToTraceEnd - sizeToRecordingEnd)));
                speedFlight(trace_i,:) = (flightPaths.batSpeed(closestIndexStart(trace_i) - preFlightPadSpeed:end + postFlightPadSpeed)+(zeros(1,sizeToTraceEnd - sizeToRecordingEnd)));
                disp('End of rec')
            end
        end
     %calculate the mean neural activity across all flights in a cluster for each cell
        meanTraceFlight(cell_i,:) = mean(traceFlight);
        meanTracePre(cell_i,:) = mean(tracePre);
        meanTracePost(cell_i,:) = mean(tracePost);
        meanSpeedFlight = mean(speedFlight);
        meanSpeedPre = mean(speedPre);
        meanSpeedPost = mean(speedPost);

        
        %smooth and zscore the neural data. subtract the min of the zscore so the
        %min is 0 rather than mean 0
        normTraceFlight(cell_i,:) = zscore(smooth(meanTraceFlight(cell_i,:),3));
        normTraceFlight(cell_i,:) = normTraceFlight(cell_i,:) - min(normTraceFlight(cell_i,:));
        smoothSpeedFlight = smooth(meanSpeedFlight,40);
        normTracePre(cell_i,:) = zscore(smooth(meanTracePre(cell_i,:),3));
        normTracePre(cell_i,:) = normTracePre(cell_i,:) - min(normTracePre(cell_i,:));       
        smoothSpeedPre = smooth(meanSpeedPre,40);
        normTracePost(cell_i,:) = zscore(smooth(meanTracePost(cell_i,:),3));
        normTracePost(cell_i,:) = normTracePost(cell_i,:) - min(normTracePost(cell_i,:));       
        smoothSpeedPost = smooth(meanSpeedPost,40);
        %find time index of max peaks
        [~,maxnormTraceFlight(cell_i,1)] = max(normTraceFlight(cell_i,:));
        [~,maxnormTracePre(cell_i,1)] = max(normTracePre(cell_i,:));
        [~,maxnormTracePost(cell_i,1)] = max(normTracePost(cell_i,:));
end    
 %sort each cell by the timing of its peak firing
    [BFlight,IFlight] = sort(maxnormTraceFlight);
    [BPre,IPre] = sort(maxnormTracePre);
    [BPost,IPost] = sort(maxnormTracePost);    
    
    %concatenate all clusters together and then zscore to equalize the
    %peaks across all trials
    meanTracePreFlightPostAll = [meanTracePre meanTraceFlight meanTracePost]; %concatenate all traces
   

%% this is to smooth, zscore, and sort the entire cell data by their preferred flight according to a homemade k-means (max dff across flights)
%zscore the full data set and subtract min to start at 0
for cell_ii = 1:length(cellData.results.C(:,1))
    normMeanTracePreFlightPostAll(cell_ii,:) = zscore(smooth(meanTracePreFlightPostAll(cell_ii,:),10));
    normMeanTracePreFlightPostAll(cell_ii,:) = normMeanTracePreFlightPostAll(cell_ii,:) - min(normMeanTracePreFlightPostAll(cell_ii,:));
end
%split the data back into its pre/during/post flight groups
%regroup by each flight group
normTraceFlightAll = normMeanTracePreFlightPostAll(:,preFlightPadCalcium+1:preFlightPadCalcium+maxDur+prePadCalcium+postPadCalcium);
normTracePreAll = normMeanTracePreFlightPostAll(:,1:preFlightPadCalcium);
normTracePostAll = normMeanTracePreFlightPostAll(:,preFlightPadCalcium+maxDur+1+prePadCalcium+postPadCalcium:preFlightPadCalcium++prePadCalcium+postPadCalcium+maxDur+postFlightPadCalcium);
%find the order of the maximum for each flight group within the regrouping
for cell_iii = 1:length(cellData.results.C(:,1))
    [~,maxNormFlightAll(cell_iii,1)] = max(normTraceFlightAll(cell_iii,:));
    [~,maxNormPreAll(cell_iii,1)] = max(normTracePreAll(cell_iii,:));
    [~,maxNormPostAll(cell_iii,1)] = max(normTracePostAll(cell_iii,:));
end
%sort each cell by the timing of its peak firing
[BnormFlightAll,InormFlightAll] = sort(maxNormFlightAll);
[BnormPreAll,InormPreAll] = sort(maxNormPreAll);
[BnormPosttAll,InormPostAll] = sort(maxNormPostAll);


%find max and sort based on the peaks of the cell across the whole
%timeseries
[~,maxAllPreFlightPost] = max(normMeanTracePreFlightPostAll,[],2);


%kPeaks = kmeans(maxAll,nClusters); %take k-means cluster of the peaks of each cell whole time
%sort based on each clustered peaks for flight, pre and post
rng(2);
for n = 1:size(maxAllPreFlightPost,1);    
    if maxAllPreFlightPost(n) >= 1 & maxAllPreFlightPost(n) < prePadCalcium +maxDur+preFlightPadCalcium+postFlightPadCalcium+postPadCalcium+1; %if the peak is within each cluster, sort it into that particular cluster
        kPeaksPreFlightPost(:,n) = 1;
    end
end

%sort the clusters according to their peak times
[BkPeaksPreFlightPost,IkPeaksPreFlightPost] = sort(kPeaksPreFlightPost);
normMeanTraceSortPreFlightPost = normMeanTracePreFlightPostAll(IkPeaksPreFlightPost,:);

%for each cluster, find the peaks of the cluster (aMat aka aTemp), tke the max
%(maxInd), and sort that max (aSort) from the cluster, add this to the normMeanTraceSort variable
aMatFlight = find(BkPeaksPreFlightPost == 1);
aTempFlight = normMeanTraceSortPreFlightPost(aMatFlight,:);
[~,maxIndFlight] = max(aTempFlight,[],2);
[BmaxFlight,ImaxFlight] = sort(maxIndFlight);
aSortFlight = aTempFlight(ImaxFlight,:);
normMeanTraceSortPreFlightPost(aMatFlight,:) = aSortFlight;
clear aMatFlight aTempFlight maxIndFlight ImaxFlight BaxFlight aSortFlight

%divide the normMeanTraceSort overall variable into the respective clusters
normMeanTraceEachFlight = normMeanTraceSortPreFlightPost(:,preFlightPadCalcium+1:preFlightPadCalcium+prePadCalcium+postPadCalcium+maxDur);
normMeanTraceEachPre = normMeanTraceSortPreFlightPost(:,1:preFlightPadCalcium);
normMeanTraceEachPost = normMeanTraceSortPreFlightPost(:,preFlightPadCalcium+prePadCalcium+postPadCalcium+maxDur+1:preFlightPadCalcium+prePadCalcium+postPadCalcium+maxDur+postFlightPadCalcium);

%% save to snakeTrace variable
snakeTracePrePost.normMeanTraceSortPreFlightPost = normMeanTraceSortPreFlightPost;
snakeTracePrePost.meanTracePreFlightPostAll=meanTracePreFlightPostAll;
snakeTracePrePost.maxAllPreFlightPost=maxAllPreFlightPost;
snakeTracePrePost.normMeanTracePreFlightPostAll=normMeanTracePreFlightPostAll;
snakeTracePrePost.kPeaksPreFlightPost=kPeaksPreFlightPost;
snakeTracePrePost.IkPeaksPreFlightPost=IkPeaksPreFlightPost;
snakeTracePrePost.normMeanTraceSortPreFlightPost = normMeanTraceSortPreFlightPost;
snakeTracePrePost.meanTraceFlight = meanTraceFlight;
snakeTracePrePost.normTraceFlight = normTraceFlight;
snakeTracePrePost.maxnormTraceFlight = maxnormTraceFlight;
snakeTracePrePost.meanTracePre = meanTracePre;
snakeTracePrePost.normTracePre = normTracePre;
snakeTracePrePost.maxnormTracePre = maxnormTracePre;
snakeTracePrePost.meanTracePost = meanTracePost;
snakeTracePrePost.normTracePost = normTracePost;
snakeTracePrePost.maxnormTracePost = maxnormTracePost;
snakeTracePrePost.BFlight = BFlight;
snakeTracePrePost.IFlight = IFlight;
snakeTracePrePost.BPre = BPre;
snakeTracePrePost.IPre = IPre;
snakeTracePrePost.BPost = BPost;
snakeTracePrePost.IPost = IPost;
snakeTracePrePost.meanSpeedFlight = meanSpeedFlight;
snakeTracePrePost.smoothSpeedFlight = smoothSpeedFlight;
snakeTracePrePost.meanSpeedPre = meanSpeedPre;
snakeTracePrePost.smoothSpeedPre = smoothSpeedPre;
snakeTracePrePost.meanSpeedPost = meanSpeedPost;
snakeTracePrePost.smoothSpeedPost = smoothSpeedPost;
snakeTracePrePost.smoothSpeedRawPre = smoothSpeedRawPre;
snakeTracePrePost.smoothSpeedRawPost = smoothSpeedRawPost;
snakeTracePrePost.smoothSpeedRawFlight = smoothSpeedRawFlight;
snakeTracePrePost.normMeanTraceEachFlight = normMeanTraceEachFlight;
snakeTracePrePost.normMeanTraceEachPre = normMeanTraceEachPre;
snakeTracePrePost.normMeanTraceEachPost = normMeanTraceEachPost;
snakeTracePrePost.batName = batName;
snakeTracePrePost.dateSesh = dateSesh;
snakeTracePrePost.sessionType = sessionType;
snakeTracePrePost.InormFlightAll = InormFlightAll;
snakeTracePrePost.normTraceFlightAll = normTraceFlightAll;

if saveFlag == 1
    snakeTracePrePost.label = label;
   save([pwd '/analysis/' label '_snakePlotData.mat'],'snakeTracePrePost');
end