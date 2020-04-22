function [snakeTrace_cRaw,snakeTrace_c,snakeTrace_s] = ImBat_snakeData(cellData,flightPaths,alignment,varargin)

batName = [];
dateSesh = [];
sessionType = [];
loadFlag = 0; %do you want to load and save the data individually outside of ImBatAnalyze
saveFlag = 0;
%offset = 0.1; % account for slow calcium estimation ~move locations back  0ms in time... This is the knob to turn for 'prospective' coding...
%clusters = [1 2 4 6];
%nClusters = 4;%length(clusters); %number of flight trajectories to look at and for k means clustering of whole time series by peak



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
        case 'nclusters'
            nClusters = varargin{i+1};
        case 'preflightpad'
            preFlightPad = varargin{i+1};
        case 'postflightpad'
            postFlightPad = varargin{i+1};
        case 'analysisfolder'
            analysis_Folder = varargin{i+1};
    end
end

%labels for loading and saving data if running independent fromImBat_analyze
if loadFlag == 1
    date = strcat(lower(batName(1:2)),dateSesh);
    label = [batName '_' dateSesh '_' sessionType];
    
    cellData = load([pwd '/processed/Motion_corrected_Data_DS_results.mat']);
    alignment = load([pwd '/processed/Alignment.mat']);
    load([pwd '/' analysis_Folder '/' label '_flightPaths.mat']);
end
%padding for during flight snake plot to include some time before and after flight
prePad = 3; %number of seconds to plot before alignment point
postPad = 5; %number of seconds to plot after alignment point
prePadCalcium = prePad*cellData.results.Fs; %number of frames (seconds*freq) to include in the trace extraction
postPadCalcium = postPad*cellData.results.Fs; %add 2 seconds to the end of the plots to include delay in peak time
prePadSpeed = prePad*120; %add 2 seconds * FS of tracking data (120)
postPadSpeed = postPad*120;%add 6 seconds * FS of tracking data (120)
%padding for pre and post flight snakePlots
preFlightPad = 3; %number of seconds to include before flight starts
postFlightPad = 5; %of of seconds to include after flight ends
preFlightPadCalcium = preFlightPad*cellData.results.Fs; %number of frames (seconds*freq) to include in the trace extraction
postFlightPadCalcium = postFlightPad*cellData.results.Fs; %add 2 seconds to the end of the plots to include delay in peak time
preFlightPadSpeed = preFlightPad*120; %add 2 seconds * FS of tracking data (120)
postFlightPadSpeed = postFlightPad*120;%add 6 seconds * FS of tracking data (120)
smoothSpeed = 100; %number of tracking frames to smooth the data by using mean smooth
smoothTrace = 2; %number of calcium video frames to smooth the data
%number of clusters to plot data against
if isfield(flightPaths,'nClusters')
    nClusters = flightPaths.nClusters;
else
    nClusters = 3;
end

% if ~isempty(clustDelete) & clustDelete <6
%     flightPaths.clusterIndex{clustDelete} = [];
% end
%flightPaths.clusterIndex cat(2,flightPaths.clusterIndex{clustComb(1)},flightPaths.clusterIndex{clustComb(2)});

 %flightPaths.clusterIndex{3} =flightPaths.clusterIndex{6}; 
% flightPaths.clusterIndex{3} =[];
%flightPaths.clusterIndex= flightPaths.clusterIndex(~cellfun('isempty',flightPaths.clusterIndex));

%for each C_raw and C matrix
data_i = 1; %start with the cRaw data
for data_i = 1:3
    if data_i ==1
        traceData = cellData.results.C_raw; %assign the data from cnmfe output for C_raw matrix
        traceSmooth = 3; %number of imaging frames to smooth the data by using mean smooth
    elseif data_i ==2
        traceData = cellData.results.C; %assign the data from cnmfe output for C matrix
        traceSmooth = 3; %number of imaging frames to smooth the data by using mean smooth
    elseif data_i ==3
        traceData = cellData.results.S; %assign the data from cnmfe output for C matrix
        traceSmooth = 1; %number of imaging frames to smooth the data by using mean smooth
    end
    
    %initialize variables
    %meanTrace = cell(1,length(flightPaths.clusterIndex));
    meanTraceFlightAll = []; %initialize the variable to concatenate all traces for uniform zscore
    meanTracePreAll = []; %initialize the variable to concatenate all traces for uniform zscore
    meanTracePostAll = []; %initialize the variable to concatenate all traces for uniform zscore
    meanTracePreFlightPostAll = []; %initialize the variable to concatenate all pre/post/flight traces for uniform zscore
    traceFlight = cell(1,nClusters);
    tracePre = cell(1,nClusters);
    tracePost = cell(1,nClusters);
    speedFlight = cell(1,nClusters);
    speedPre = cell(1,nClusters);
    speedPost = cell(1,nClusters);
    meanTraceFlight = cell(1,nClusters);
    sdTraceFlight = cell(1,nClusters);
    meanTraceOdd = cell(1,nClusters);
    meanTraceEven = cell(1,nClusters);
    meanTracePreFlightPost = cell(1,nClusters);
    sdTracePreFlightPost = cell(1,nClusters);
    normTraceFlight = cell(1,nClusters);
    normTraceOdd= cell(1,nClusters);
    normTraceEven = cell(1,nClusters);
    tracePreFlightPost = cell(1,nClusters);
    smoothTraceRawFlight = cell(1,nClusters);
    normTraceRawFlight = cell(1,nClusters);
    smoothTraceRawPre = cell(1,nClusters);
    normTraceRawPre = cell(1,nClusters);
    smoothTraceRawPost = cell(1,nClusters);
    normTraceRawPost = cell(1,nClusters);
    smoothTraceRawPreFlightPost = cell(1,nClusters);
    normTraceRawPreFlightPost = cell(1,nClusters);
    speedFlight = cell(1,nClusters);
    smoothSpeedRawFlight = cell(1,nClusters);
    smoothSpeedRawPre = cell(1,nClusters);
    speedPre = cell(1,nClusters);
    smoothSpeedRawPost = cell(1,nClusters);
    speedPost = cell(1,nClusters);
    speedPreFlightPost = cell(1,nClusters);
    smoothSpeedRawPreFlightPost = cell(1,nClusters);
    normTracePre = cell(1,nClusters);
    smoothSpeedPre = cell(1,nClusters);
    normTracePost = cell(1,nClusters);
    normTracePost = cell(1,nClusters);
    smoothSpeedPost = cell(1,nClusters);
    semTracePreFlightPost = cell(1,nClusters);
    normMeanTraceFlightAll = [];
    normMeanTracePreAll= [];
    normMeanTracePreAll= [];
    normMeanTracePostAll = [];
    normMeanTracePreFlightPostAll = cell(1,nClusters);
    %for each cluster type
    for clust_i = 1:nClusters %length(flightPaths.clusterIndex)
        %for each trial in cluster, calculate the start/stop times and duration of the calcium videos
        for dur_i = 1:length(flightPaths.clusterIndex{clust_i})
            %get imaging times of start and stop index converting from tracking to video times
            [minValueStart(dur_i),closestIndexStart(dur_i)] = min(abs(alignment.out.video_timesDS-alignment.out.Location_time(flightPaths.flight_starts_idx(flightPaths.clusterIndex{clust_i}(dur_i)))));
            [minValueEnd(dur_i),closestIndexEnd(dur_i)] = min(abs(alignment.out.video_timesDS-alignment.out.Location_time(flightPaths.flight_ends_idx(flightPaths.clusterIndex{clust_i}(dur_i)))));
            %calculate duration of each flight in a particular cluster so
            %you can pad all flights to the longest flight in that cluster
            dur(dur_i) = closestIndexEnd(dur_i)-closestIndexStart(dur_i);
            durSpeed(dur_i)= flightPaths.flight_ends_idx(flightPaths.clusterIndex{clust_i}(dur_i))-flightPaths.flight_starts_idx(flightPaths.clusterIndex{clust_i}(dur_i));
        end
        %calculate max duration for each cluster of trajectories
        maxDur = max(dur);
        maxDurSpeed = max(durSpeed);
        %meanTrace{clust_i}=zeros(length(traceData(:,1)),maxDur+preWindow+1);
        %initialize the vector to store the neural activity of each flight
        traceFlight{clust_i} = zeros(length(flightPaths.clusterIndex{clust_i}),maxDur+prePadCalcium+postPadCalcium+1,length(traceData(:,1)));
        speedFlight{clust_i} = zeros(length(flightPaths.clusterIndex{clust_i}),maxDurSpeed+prePadSpeed+postPadSpeed+1);
        tracePre{clust_i} = zeros(length(flightPaths.clusterIndex{clust_i}),preFlightPadCalcium+1,length(traceData(:,1)));
        speedPre{clust_i} = zeros(length(flightPaths.clusterIndex{clust_i}),preFlightPadSpeed+1);
        tracePost{clust_i} = zeros(length(flightPaths.clusterIndex{clust_i}),postFlightPadCalcium+1,length(traceData(:,1)));
        speedPost{clust_i} = zeros(length(flightPaths.clusterIndex{clust_i}),postFlightPadSpeed+1);
        tracePreFlightPost{clust_i} = zeros(length(flightPaths.clusterIndex{clust_i}),length(tracePre{clust_i}(1,:,1)) + length(traceFlight{clust_i}(1,:,1)) + length(tracePost{clust_i}(1,:,1)));
        speedPreFlightPost{clust_i} = zeros(length(flightPaths.clusterIndex{clust_i}),length(speedPre{clust_i}(1,:)) + length(speedFlight{clust_i}(1,:)) + length(speedPost{clust_i}(1,:)));
        %for each cell
        for cell_i = 1:length(traceData(:,1))
            %build the calcium and speed vectors for each flight within each cluster
            for trace_i = 1:length(flightPaths.clusterIndex{clust_i})
                try
                    traceFlightIdx{clust_i}(trace_i) = flightPaths.clusterIndex{clust_i}(trace_i);
                    traceFlight{clust_i}(trace_i,:,cell_i) = traceData(cell_i,closestIndexStart(trace_i) - prePadCalcium:closestIndexEnd(trace_i) + (maxDur-dur(trace_i)) + postPadCalcium);
                    speedFlight{clust_i}(trace_i,:) = flightPaths.batSpeed(flightPaths.flight_starts_idx(flightPaths.clusterIndex{clust_i}(trace_i)) - prePadSpeed:flightPaths.flight_ends_idx(flightPaths.clusterIndex{clust_i}(trace_i)) + (maxDurSpeed-durSpeed(trace_i)) + postPadSpeed);
                    smoothTraceRawFlight{clust_i}(trace_i,:,cell_i) = smooth(traceFlight{clust_i}(trace_i,:,cell_i),traceSmooth);
                    normTraceRawFlight{clust_i}(trace_i,:,cell_i) = zscore(smoothTraceRawFlight{clust_i}(trace_i,:,cell_i),0,smoothTrace);
                    normTraceRawFlight{clust_i}(trace_i,:,cell_i) = normTraceRawFlight{clust_i}(trace_i,:,cell_i) - min(normTraceRawFlight{clust_i}(trace_i,:,cell_i));
                    smoothSpeedRawFlight{clust_i}(trace_i,:) = smooth(speedFlight{clust_i}(trace_i,:),smoothSpeed);
                    tracePre{clust_i}(trace_i,:,cell_i) = traceData(cell_i,closestIndexStart(trace_i) - preFlightPadCalcium:closestIndexStart(trace_i));
                    speedPre{clust_i}(trace_i,:) = flightPaths.batSpeed(flightPaths.flight_starts_idx(flightPaths.clusterIndex{clust_i}(trace_i)) - preFlightPadSpeed:flightPaths.flight_starts_idx(flightPaths.clusterIndex{clust_i}(trace_i)));
                    smoothTraceRawPre{clust_i}(trace_i,:,cell_i) = smooth(tracePre{clust_i}(trace_i,:,cell_i),traceSmooth);
                    normTraceRawPre{clust_i}(trace_i,:,cell_i) = zscore(smoothTraceRawPre{clust_i}(trace_i,:,cell_i),0,smoothTrace);
                    normTraceRawPre{clust_i}(trace_i,:,cell_i) = normTraceRawPre{clust_i}(trace_i,:,cell_i) - min(normTraceRawPre{clust_i}(trace_i,:,cell_i));
                    smoothSpeedRawPre{clust_i}(trace_i,:) = smooth(speedPre{clust_i}(trace_i,:),smoothSpeed);
                    tracePost{clust_i}(trace_i,:,cell_i) = traceData(cell_i,closestIndexEnd(trace_i):closestIndexEnd(trace_i) + postFlightPadCalcium);
                    speedPost{clust_i}(trace_i,:) = flightPaths.batSpeed(flightPaths.flight_ends_idx(flightPaths.clusterIndex{clust_i}(trace_i))+(maxDurSpeed-durSpeed(trace_i)):flightPaths.flight_ends_idx(flightPaths.clusterIndex{clust_i}(trace_i)) + (maxDurSpeed-durSpeed(trace_i)) + postFlightPadSpeed);
                    smoothTraceRawPost{clust_i}(trace_i,:,cell_i) = smooth(tracePost{clust_i}(trace_i,:,cell_i),traceSmooth);
                    normTraceRawPost{clust_i}(trace_i,:,cell_i) = zscore(smoothTraceRawPost{clust_i}(trace_i,:,cell_i),0,smoothTrace);
                    normTraceRawPost{clust_i}(trace_i,:,cell_i) = normTraceRawPost{clust_i}(trace_i,:,cell_i) - min(normTraceRawPost{clust_i}(trace_i,:,cell_i));
                    smoothSpeedRawPost{clust_i}(trace_i,:) = smooth(speedPost{clust_i}(trace_i,:),smoothSpeed);
                    tracePreFlightPost{clust_i}(trace_i,:,cell_i) = cat(2,tracePre{clust_i}(trace_i,:,cell_i),traceFlight{clust_i}(trace_i,:,cell_i),tracePost{clust_i}(trace_i,:,cell_i));
                    speedPreFlightPost{clust_i}(trace_i,:) = horzcat(speedPre{clust_i}(trace_i,:),speedFlight{clust_i}(trace_i,:),speedPost{clust_i}(trace_i,:));
                    smoothTraceRawPreFlightPost{clust_i}(trace_i,:,cell_i) = smooth(tracePreFlightPost{clust_i}(trace_i,:,cell_i),traceSmooth);
                    normTraceRawPreFlightPost{clust_i}(trace_i,:,cell_i) = zscore(smoothTraceRawPreFlightPost{clust_i}(trace_i,:,cell_i),0,smoothTrace);
                    normTraceRawPreFlightPost{clust_i}(trace_i,:,cell_i) = normTraceRawPreFlightPost{clust_i}(trace_i,:,cell_i) - min(normTraceRawPreFlightPost{clust_i}(trace_i,:,cell_i));
                    smoothSpeedRawPreFlightPost{clust_i}(trace_i,:) = smooth(speedPreFlightPost{clust_i}(trace_i,:),smoothSpeed);
                    
                    
                catch
                    %                     tracePost{clust_i}(trace_i,:,cell_i) = [traceData(cell_i,closestIndexEnd(trace_i):end),zeros(1,size(tracePost{clust_i}(trace_i,:,cell_i),2) - size(traceData(cell_i,closestIndexEnd(trace_i):end),2))];
                    %                     speedPost{clust_i}(trace_i,:) = [flightPaths.batSpeed(flightPaths.flight_ends_idx(flightPaths.clusterIndex{clust_i}(trace_i))+(maxDurSpeed-durSpeed(trace_i)):end),zeros(1,size(speedPost{clust_i}(trace_i,:),2)-size(flightPaths.batSpeed(flightPaths.flight_ends_idx(flightPaths.clusterIndex{clust_i}(trace_i))+(maxDurSpeed-durSpeed(trace_i)):end,2)))];
                    %                     smoothTraceRawPost{clust_i}(trace_i,:,cell_i) = smooth(tracePost{clust_i}(trace_i,:,cell_i),traceSmooth);
                    %                     normTraceRawPost{clust_i}(trace_i,:,cell_i) = zscore(smoothTraceRawPost{clust_i}(trace_i,:,cell_i),0,2);
                    %                     normTraceRawPost{clust_i}(trace_i,:,cell_i) = normTraceRawPost{clust_i}(trace_i,:,cell_i) - min(normTraceRawPost{clust_i}(trace_i,:,cell_i));
                    %                     smoothSpeedRawPost{clust_i}(trace_i,:) = smooth(speedPost{clust_i}(trace_i,:),speedSmooth);
                    disp('End of rec')
                end
                
            end
            
            %calculate the mean neural activity across all flights in a cluster for each cell
            if length(traceFlight{clust_i}(:,1,cell_i)) <2
                meanTraceFlight{clust_i}(cell_i,:) = traceFlight{clust_i}(:,:,cell_i);
                meanTracePre{clust_i}(cell_i,:) = tracePre{clust_i}(:,:,cell_i);
                meanTracePost{clust_i}(cell_i,:) = tracePost{clust_i}(:,:,cell_i);
            else 
                meanTraceFlight{clust_i}(cell_i,:) = mean(traceFlight{clust_i}(:,:,cell_i));
                meanTracePre{clust_i}(cell_i,:) = mean(tracePre{clust_i}(:,:,cell_i));
                meanTracePost{clust_i}(cell_i,:) = mean(tracePost{clust_i}(:,:,cell_i));
                
            end
            sdTraceFlight{clust_i}(cell_i,:) = std(traceFlight{clust_i}(:,:,cell_i));
            meanTraceOdd{clust_i}(cell_i,:) = mean(traceFlight{clust_i}(1:2:end,:,cell_i));
            meanTraceEven{clust_i}(cell_i,:) = mean(traceFlight{clust_i}(2:2:end,:,cell_i));
            meanSpeedFlight{clust_i} = mean(speedFlight{clust_i});
            meanSpeedPre{clust_i} = mean(speedPre{clust_i});
            meanSpeedPost{clust_i} = mean(speedPost{clust_i});
            meanTracePreFlightPost{clust_i}(cell_i,:) = mean(tracePreFlightPost{clust_i}(:,:,cell_i));
            sdTracePreFlightPost{clust_i}(cell_i,:) = std(tracePreFlightPost{clust_i}(:,:,cell_i));
            semTracePreFlightPost{clust_i}(cell_i,:) = std(tracePreFlightPost{clust_i}(:,:,cell_i))/sqrt(length(tracePreFlightPost{clust_i}(:,:,cell_i)));
            meanSpeedPreFlightPost{clust_i} = mean(speedPreFlightPost{clust_i});
            
            %smooth and zscore the neural data. subtract the min of the zscore so the
            %min is 0 rather than mean 0
            normTraceFlight{clust_i}(cell_i,:) = zscore(smooth(meanTraceFlight{clust_i}(cell_i,:),traceSmooth));
            normTraceFlight{clust_i}(cell_i,:) = normTraceFlight{clust_i}(cell_i,:) - min(normTraceFlight{clust_i}(cell_i,:));
            smoothSpeedFlight{clust_i} = smooth(meanSpeedFlight{clust_i},smoothSpeed);
            normTraceOdd{clust_i}(cell_i,:) = zscore(smooth(meanTraceOdd{clust_i}(cell_i,:),traceSmooth));
            normTraceOdd{clust_i}(cell_i,:) = normTraceOdd{clust_i}(cell_i,:) - min(normTraceOdd{clust_i}(cell_i,:));
            normTraceEven{clust_i}(cell_i,:) = zscore(smooth(meanTraceEven{clust_i}(cell_i,:),traceSmooth));
            normTraceEven{clust_i}(cell_i,:) = normTraceEven{clust_i}(cell_i,:) - min(normTraceEven{clust_i}(cell_i,:));
            normTracePre{clust_i}(cell_i,:) = zscore(smooth(meanTracePre{clust_i}(cell_i,:),traceSmooth));
            normTracePre{clust_i}(cell_i,:) = normTracePre{clust_i}(cell_i,:) - min(normTracePre{clust_i}(cell_i,:));
            smoothSpeedPre{clust_i} = smooth(meanSpeedPre{clust_i},smoothSpeed);
            normTracePost{clust_i}(cell_i,:) = zscore(smooth(meanTracePost{clust_i}(cell_i,:),traceSmooth));
            normTracePost{clust_i}(cell_i,:) = normTracePost{clust_i}(cell_i,:) - min(normTracePost{clust_i}(cell_i,:));
            smoothSpeedPost{clust_i} = smooth(meanSpeedPost{clust_i},smoothSpeed);
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
        meanTracePreFlightPostAll{clust_i} = [meanTracePre{clust_i} meanTraceFlight{clust_i} meanTracePost{clust_i}];
        
        
    end
    
    
    %% this is to smooth, zscore, and sort the entire cell data by their preferred flight according to a homemade k-means (max dff across flights)
    %zscore the full data set and subtract min to start at 0
    for cell_ii = 1:length(traceData(:,1))
        normMeanTraceFlightAll(cell_ii,:) = zscore(smooth(meanTraceFlightAll(cell_ii,:),traceSmooth));
        normMeanTraceFlightAll(cell_ii,:) = normMeanTraceFlightAll(cell_ii,:) - min(normMeanTraceFlightAll(cell_ii,:));
        normMeanTracePreAll(cell_ii,:) = zscore(smooth(meanTracePreAll(cell_ii,:),traceSmooth));
        normMeanTracePreAll(cell_ii,:) = normMeanTracePreAll(cell_ii,:) - min(normMeanTracePreAll(cell_ii,:));
        normMeanTracePostAll(cell_ii,:) = zscore(smooth(meanTracePostAll(cell_ii,:),traceSmooth));
        normMeanTracePostAll(cell_ii,:) = normMeanTracePostAll(cell_ii,:) - min(normMeanTracePostAll(cell_ii,:));
        for clust_i = 1:nClusters
            normMeanTracePreFlightPostAll{clust_i}(cell_ii,:) = zscore(smooth(meanTracePreFlightPostAll{clust_i}(cell_ii,:),traceSmooth));
            normMeanTracePreFlightPostAll{clust_i}(cell_ii,:) = normMeanTracePreFlightPostAll{clust_i}(cell_ii,:) - min(normMeanTracePreFlightPostAll{clust_i}(cell_ii,:));
        end
    end
    %split the data back into its clusters
    traceFlightIndConcat = [1 length(normTraceFlight{1}(1,:))]; %initialize index variable to hold the indices for start/stop of each cluster
    tracePreIndConcat = [1 length(normTracePre{1}(1,:))]; %initialize index variable to hold the indices for start/stop of each cluster
    tracePostIndConcat = [1 length(normTracePost{1}(1,:))]; %initialize index variable to hold the indices for start/stop of each cluster
    
    %find the start and stop indices for 2nd through n cluster
    for clust_ii = 2:nClusters
        traceFlightIndConcat = vertcat(traceFlightIndConcat,[traceFlightIndConcat(clust_ii-1,2)+1 traceFlightIndConcat(clust_ii-1,2)+length(normTraceFlight{clust_ii}(1,:))]);
        tracePreIndConcat = vertcat(tracePreIndConcat,[tracePreIndConcat(clust_ii-1,2)+1 tracePreIndConcat(clust_ii-1,2)+length(normTracePre{clust_ii}(1,:))]);
        tracePostIndConcat = vertcat(tracePostIndConcat,[tracePostIndConcat(clust_ii-1,2)+1 tracePostIndConcat(clust_ii-1,2)+length(normTracePost{clust_ii}(1,:))]);
        
    end
    %regroup by each flight cluster
    for clust_iii = 1:nClusters
        normTraceFlightAll{clust_iii} = normMeanTraceFlightAll(:,traceFlightIndConcat(clust_iii,1):traceFlightIndConcat(clust_iii,2));
        normTracePreAll{clust_iii} = normMeanTracePreAll(:,tracePreIndConcat(clust_iii,1):tracePreIndConcat(clust_iii,2));
        normTracePostAll{clust_iii} = normMeanTracePostAll(:,tracePostIndConcat(clust_iii,1):tracePostIndConcat(clust_iii,2));
        normTracePreClustAll{clust_iii} = normMeanTracePreFlightPostAll{clust_iii}(:,1:preFlightPadCalcium);
        normTraceFlightClustAll{clust_iii} = normMeanTracePreFlightPostAll{clust_iii}(:,preFlightPadCalcium+1:preFlightPadCalcium+length(normTraceFlightAll{clust_iii}(1,:)));
        normTracePostClustAll{clust_iii} = normMeanTracePreFlightPostAll{clust_iii}(:,preFlightPadCalcium+length(normTraceFlightAll{clust_iii}(1,:))+1:preFlightPadCalcium+length(normTraceFlightAll{clust_iii}(1,:))+postFlightPadCalcium);
        %find the order of the maximum for each cluster within the regrouping
        for cell_iii = 1:length(traceData(:,1))
            [~,maxNormFlightAll{clust_iii}(cell_iii,1)] = max(normTraceFlightAll{clust_iii}(cell_iii,:));
            [~,maxNormPreAll{clust_iii}(cell_iii,1)] = max(normTracePreAll{clust_iii}(cell_iii,:));
            [~,maxNormPostAll{clust_iii}(cell_iii,1)] = max(normTracePostAll{clust_iii}(cell_iii,:));
            [~,maxNormPreClustAll{clust_iii}(cell_iii,1)] = max(normTracePreClustAll{clust_iii}(cell_iii,:));
            [~,maxNormFlightClustAll{clust_iii}(cell_iii,1)] = max(normTraceFlightClustAll{clust_iii}(cell_iii,:));
            [~,maxNormPostClustAll{clust_iii}(cell_iii,1)] = max(normTracePostClustAll{clust_iii}(cell_iii,:));
        end
        %sort each cell by the timing of its peak firing
        [BnormFlightAll{clust_iii},InormFlightAll{clust_iii}] = sort(maxNormFlightAll{clust_iii});
        [B1normFlightAll{clust_iii},I1normFlightAll{clust_iii}] = sort(maxNormFlightAll{1}); %sort by cluster 1
        [BnormPreAll{clust_iii},InormPreAll{clust_iii}] = sort(maxNormPreAll{clust_iii});
        [B1normPreAll{clust_iii},I1normPreAll{clust_iii}] = sort(maxNormPreAll{1}); %sort by cluster 1
        [BnormPostAll{clust_iii},InormPostAll{clust_iii}] = sort(maxNormPostAll{clust_iii});
        [B1normPostAll{clust_iii},I1normPostAll{clust_iii}] = sort(maxNormPostAll{1}); %sort by cluster 1
        
        [BnormPreClustAll{clust_iii},InormPreClustAll{clust_iii}] = sort(maxNormPreClustAll{clust_iii});
        [B1normPreClustAll{clust_iii},I1normPreClustAll{clust_iii}] = sort(maxNormPreClustAll{1}); %sort by cluster 1
        [BnormFlightClustAll{clust_iii},InormFlightClustAll{clust_iii}] = sort(maxNormFlightClustAll{clust_iii});
        [B1FlightClustAll{clust_iii},I1FlightClustAll{clust_iii}] = sort(maxNormFlightClustAll{1}); %sort by cluster 1
        [BnormPostClustAll{clust_iii},InormPostClustAll{clust_iii}] = sort(maxNormPostClustAll{clust_iii});
        [B1normPostClustAll{clust_iii},I1normPostClustAll{clust_iii}] = sort(maxNormPostClustAll{1}); %sort by cluster 1
        
        
        
        
    end
    
    
    %find max and sort based on the peaks of the cell across the whole
    %timeseries
    [~,maxAllFlight] = max(normMeanTraceFlightAll,[],2);
    [~,maxAllPre] = max(normMeanTracePreAll,[],2);
    [~,maxAllPost] = max(normMeanTracePostAll,[],2);
    for clust_i = 1:nClusters
        [~,maxAllPreFlightPost{clust_i}] = max(normMeanTracePreFlightPostAll{clust_iii},[],2);
    end
    
    
    %kPeaks = kmeans(maxAll,nClusters); %take k-means cluster of the peaks of each cell whole time
    %sort based on each clustered peaks for flight, pre and post
    rng(2);
    for n = 1:size(maxAllFlight,1)
        for ii = 1:nClusters
            if maxAllFlight(n) >= traceFlightIndConcat(ii,1) & maxAllFlight(n) < traceFlightIndConcat(ii,2); %if the peak is within each cluster, sort it into that particular cluster
                kPeaksFlight(:,n) = 1;
            end
            
        end
    end
    for n = 1:size(maxAllPre,1)
        for ii = 1:nClusters
            if maxAllPre(n) >= tracePreIndConcat(ii,1) & maxAllPre(n) < tracePreIndConcat(ii,2); %if the peak is within each cluster, sort it into that particular cluster
                kPeaksPre(:,n) = 1;
            end
        end
    end
    for n = 1:size(maxAllPost,1)
        for ii = 1:nClusters
            if maxAllPost(n) >= tracePostIndConcat(ii,1) & maxAllPost(n) < tracePostIndConcat(ii,2); %if the peak is within each cluster, sort it into that particular cluster
                kPeaksPost(:,n) = 1;
            end
        end
    end
    for clust_i = 1:nClusters
        for n = 1:size(maxAllPreFlightPost{clust_i},1)
            if maxAllPreFlightPost{clust_i}(n) >= 1 & maxAllPreFlightPost{clust_i}(n) < (length(normTraceFlightAll{clust_iii}(1,:))+preFlightPadCalcium+postFlightPadCalcium+1) %if the peak is within each cluster, sort it into that particular cluster
                kPeaksPreFlightPost{clust_i}(:,n) = 1;
            end
        end
        [BkPeaksPreFlightPost{clust_i},IkPeaksPreFlightPost{clust_i}] = sort(kPeaksPreFlightPost{clust_i});
        normMeanTraceSortPreFlightPost{clust_i} = normMeanTracePreFlightPostAll{clust_i}(IkPeaksPreFlightPost{clust_i},:);
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
        aMatPreFlightPost = find(BkPeaksPreFlightPost{c} == 1);
        aTempPreFlightPost = normMeanTraceSortPreFlightPost{c}(aMatPreFlightPost,:);
        [~,maxIndPreFlightPost] = max(aTempPreFlightPost,[],2);
        [BmaxPreFlightPost,ImaxPreFlightPost] = sort(maxIndPreFlightPost);
        aSortPreFlightPost = aTempPreFlightPost(ImaxPreFlightPost,:);
        normMeanTraceSortPreFlightPost{c}(aMatPreFlightPost,:) = aSortPreFlightPost;
        clear aMatPreFlightPost aTempPreFlightPost maxIndPreFlightPost ImaxPreFlightPost BaxPreFlightPost aSortPreFlightPost
        
        
    end
    %divide the normMeanTraceSort overall variable into the respective clusters
    for c = 1:nClusters
        normMeanTraceEachFlight{c} = normMeanTraceSortFlight(:,traceFlightIndConcat(c,1):traceFlightIndConcat(c,2));
        normMeanTraceEachPre{c} = normMeanTraceSortPre(:,tracePreIndConcat(c,1):tracePreIndConcat(c,2));
        normMeanTraceEachPost{c} = normMeanTraceSortPost(:,tracePostIndConcat(c,1):tracePostIndConcat(c,2));
        normMeanTraceEachPFlightP{c} = normMeanTraceSortPreFlightPost{c}(:,preFlightPadCalcium+1:length(normTraceFlightAll{c}(1,:))+preFlightPadCalcium);
        normMeanTraceEachPreFP{c} = normMeanTraceSortPreFlightPost{c}(:,1:preFlightPadCalcium);
        normMeanTraceEachPFPost{c} = normMeanTraceSortPreFlightPost{c}(:,length(normTraceFlightAll{c}(1,:))+preFlightPadCalcium+1:length(normTraceFlightAll{c}(1,:))+preFlightPadCalcium+postFlightPadCalcium+1);
        
    end
    
    %% save to snakeTrace variable
        snakeTraceData.batName = batName;
        snakeTraceData.dateSesh = dateSesh;
        snakeTraceData.sessionType = sessionType;
        snakeTraceData.dur = dur;
        snakeTraceData.BPre = BPre;
        snakeTraceData.IPre = IPre;
        snakeTraceData.BFlight = BFlight;
        snakeTraceData.IFlight = IFlight;
        snakeTraceData.BPost = BPost;
        snakeTraceData.IPost = IPost;
        snakeTraceData.Bodd = Bodd;
        snakeTraceData.Iodd = Iodd;
        snakeTraceData.B1Pre = B1Pre;
        snakeTraceData.I1Pre = I1Pre;
        snakeTraceData.B1Flight = B1Flight;
        snakeTraceData.I1Flight = I1Flight;
        snakeTraceData.B1Post = B1Post;
        snakeTraceData.I1Post = I1Post;
        snakeTraceData.InormFlightAll = InormFlightAll;
        snakeTraceData.I1normFlightAll = I1normFlightAll;
        snakeTraceData.IkPeaksPreFlightPost=IkPeaksPreFlightPost;
        snakeTraceData.kPeaksPreFlightPost=kPeaksPreFlightPost;
        snakeTraceData.maxAllPreFlightPost=maxAllPreFlightPost;
        snakeTraceData.maxnormTracePre = maxnormTracePre;
        snakeTraceData.maxnormTraceFlight = maxnormTraceFlight;
        snakeTraceData.maxnormTracePost = maxnormTracePost;snakeTraceData.meanTracePre = meanTracePre;
        snakeTraceData.meanTraceFlight = meanTraceFlight;
        snakeTraceData.meanTracePost = meanTracePost;
        snakeTraceData.meanSpeedPre = meanSpeedPre;
        snakeTraceData.meanSpeedFlight = meanSpeedFlight;
        snakeTraceData.meanSpeedPost = meanSpeedPost;
        snakeTraceData.meanSpeedPreFlightPost = meanSpeedPreFlightPost;
        snakeTraceData.meanTracePreFlightPost = meanTracePreFlightPost;
        snakeTraceData.meanTracePreFlightPostAll=meanTracePreFlightPostAll;
        snakeTraceData.meanTraceEven = meanTraceEven;
        snakeTraceData.meanTraceOdd = meanTraceOdd;
        snakeTraceData.nClusters = nClusters;
        snakeTraceData.normTracePre = normTracePre;
        snakeTraceData.normTraceFlight = normTraceFlight;
        snakeTraceData.normTracePost = normTracePost;
        snakeTraceData.normTraceEven = normTraceEven;
        snakeTraceData.normTraceOdd = normTraceOdd;
        snakeTraceData.normTraceFlightAll = normTraceFlightAll;
        snakeTraceData.normTraceRawPre = normTraceRawPre;
        snakeTraceData.normTraceRawFlight = normTraceRawFlight;
        snakeTraceData.normTraceRawPost = normTraceRawPost;
        snakeTraceData.normTraceRawPreFlightPost = normTraceRawPreFlightPost;
        snakeTraceData.normMeanTraceSortPreFlightPost = normMeanTraceSortPreFlightPost;
        snakeTraceData.normMeanTracePreFlightPostAll=normMeanTracePreFlightPostAll;
        snakeTraceData.normMeanTraceEachPFlightP = normMeanTraceEachPFlightP;
        snakeTraceData.normMeanTraceEachPreFP = normMeanTraceEachPreFP;
        snakeTraceData.normMeanTraceEachPFPost = normMeanTraceEachPFPost;
        snakeTraceData.normMeanTraceEachPre = normMeanTraceEachPre;
        snakeTraceData.normMeanTraceSortPre = normMeanTraceSortPre;
        snakeTraceData.normMeanTracePreAll = normMeanTracePreAll;
        snakeTraceData.normMeanTraceEachFlight = normMeanTraceEachFlight;
        snakeTraceData.normMeanTraceSortFlight = normMeanTraceSortFlight;
        snakeTraceData.normMeanTraceFlightAll = normMeanTraceFlightAll;
        snakeTraceData.normMeanTraceEachPost = normMeanTraceEachPost;
        snakeTraceData.normMeanTraceSortPost = normMeanTraceSortPost;
        snakeTraceData.normMeanTracePostAll = normMeanTracePostAll;       
        snakeTraceData.preFlightPadCalcium = preFlightPadCalcium;
        snakeTraceData.postFlightPadCalcium = postFlightPadCalcium;
        snakeTraceData.preFlightPadSpeed = preFlightPadSpeed;
        snakeTraceData.postFlightPadSpeed = postFlightPadSpeed;
        snakeTraceData.preFlightPad = preFlightPad;
        snakeTraceData.postFlightPad = postFlightPad;
        snakeTraceData.sdTracePreFlightPost = sdTracePreFlightPost;
        snakeTraceData.semTracePreFlightPost = semTracePreFlightPost;
        snakeTraceData.smoothTraceRawPreFlightPost = smoothTraceRawPreFlightPost;
        snakeTraceData.smoothSpeedPre = smoothSpeedPre;
        snakeTraceData.smoothSpeedFlight = smoothSpeedFlight;
        snakeTraceData.smoothSpeedPost = smoothSpeedPost;
        snakeTraceData.smoothSpeedRawPre = smoothSpeedRawPre;
        snakeTraceData.smoothSpeedRawFlight = smoothSpeedRawFlight;
        snakeTraceData.smoothSpeedRawPost = smoothSpeedRawPost;
        snakeTraceData.smoothTraceRawPre = smoothTraceRawPre;
        snakeTraceData.smoothTraceRawFlight = smoothTraceRawFlight;
        snakeTraceData.smoothTraceRawPost = smoothTraceRawPost;
        snakeTraceData.smoothSpeedRawPreFlightPost = smoothSpeedRawPreFlightPost;
        snakeTraceData.speedPreFlightPost = speedPreFlightPost;
        snakeTraceData.traceFlightIndConcat = traceFlightIndConcat;
        snakeTraceData.tracePreIndConcat = tracePreIndConcat;
        snakeTraceData.tracePostIndConcat = tracePostIndConcat;
        snakeTraceData.traceFlightIdx = traceFlightIdx;
        snakeTraceData.tracePre = tracePre;
        snakeTraceData.traceFlight = traceFlight;
        snakeTraceData.tracePost = tracePost;
        snakeTraceData.tracePreFlightPost = tracePreFlightPost;
        if data_i ==1
            snakeTrace_cRaw = snakeTraceData;
        elseif data_i == 2
            snakeTrace_c = snakeTraceData;
        elseif data_i == 3
            snakeTrace_s = snakeTraceData;
        end   
end
if loadFlag ==1 && saveFlag == 1
    snakeTrace_c.label = label;
    snakeTrace_cRaw.label = label;
    snakeTrace_s.label = label;
    save([pwd '/' analysis_Folder '/' label '_snakeTraceData_cRaw.mat'],'snakeTrace_cRaw');
    save([pwd '/' analysis_Folder '/' label '_snakeTraceData_c.mat'],'snakeTrace_c');
    save([pwd '/' analysis_Folder '/' label '_snakeTraceData_s.mat'],'snakeTrace_s');
elseif loadFlag == 0 && saveFlag ==1
    save([pwd '/' batName '_' dateSesh '_' sessionType '_snakeTraceData_cRaw.mat'],'snakeTrace_cRaw');
    save([pwd '/' batName '_' dateSesh '_' sessionType '_snakeTraceData_c.mat'],'snakeTrace_c');
    save([pwd '/' batName '_' dateSesh '_' sessionType '_snakeTraceData_s.mat'],'snakeTrace_s');
end

%snakeTrace.normMeanTraceAllSmooth = normMeanTraceAllSmooth;
