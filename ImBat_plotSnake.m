function [snakeTrace] = ImBat_plotSnake(cellData,flightPaths,alignment,varargin)
%global batName dateSesh sessionType


batName = [];
dateSesh = [];
sessionType = [];
loadFlag = 0; %do you want to load and save the data individually outside of ImBatAnalyze
saveFlag = 0;
offset = 0.1; % account for slow calcium estimation ~move locations back 100ms in time... This is the knob to turn for 'prospective' coding...
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

prePad = 2*cellData.results.Fs; %number of frames (seconds*freq) to include in the trace extraction
postPad = 6*cellData.results.Fs; %add 2 seconds to the end of the plots to include delay in peak time
%meanTrace = cell(1,length(flightPaths.clusterIndex));
meanTraceAll = []; %initialize the variable to concatenate all traces for uniform zscore

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
        trace = zeros(length(flightPaths.clusterIndex{clust_i}),maxDur+prePad+postPad+1);
        speed{clust_i} = zeros(length(flightPaths.clusterIndex{clust_i}),maxDurSpeed+prePad+postPad+1);
        for trace_i = 1:length(flightPaths.clusterIndex{clust_i})
            try
                trace(trace_i,:) = cellData.results.C(cell_i,closestIndexStart(trace_i) - prePad:closestIndexEnd(trace_i) + (maxDur-dur(trace_i)) + postPad);
                speed{clust_i}(trace_i,:) = flightPaths.batSpeed(flightPaths.flight_starts_idx(flightPaths.clusterIndex{clust_i}(trace_i)) - prePad:flightPaths.flight_ends_idx(flightPaths.clusterIndex{clust_i}(trace_i)) + (maxDurSpeed-durSpeed(trace_i)) + postPad);
                smoothSpeedRaw{clust_i}(trace_i,:) = smooth(speed{clust_i}(trace_i,:),100);
            catch
                sizeToRecordingEnd = size(cellData.results.C(cell_i,closestIndexStart(trace_i) - prePad:end),2);
                sizeToTraceEnd = size(trace(trace_i,:),2);
                trace(trace_i,:) = (cellData.results.C(cell_i,closestIndexStart(trace_i) - prePad:end)+(zeros(1,sizeToTraceEnd - sizeToRecordingEnd)));
                speed{clust_i}(trace_i,:) = (flightPaths.batSpeed(closestIndexStart(trace_i) - prePad:end)+(zeros(1,sizeToTraceEnd - sizeToRecordingEnd)));
                
            end
        end
        
        %calculate the mean neural activity across all flights in a cluster for each cell
        meanTrace{clust_i}(cell_i,:) = mean(trace);
        meanTraceOdd{clust_i}(cell_i,:) = mean(trace(1:2:end,:));
        meanTraceEven{clust_i}(cell_i,:) = mean(trace(2:2:end,:));
        meanSpeed{clust_i} = mean(speed{clust_i});
        
        %smooth and zscore the neural data. subtract the min of the zscore so the
        %min is 0 rather than mean 0
        normTrace{clust_i}(cell_i,:) = zscore(smooth(meanTrace{clust_i}(cell_i,:),3));
        normTrace{clust_i}(cell_i,:) = normTrace{clust_i}(cell_i,:) - min(normTrace{clust_i}(cell_i,:));
        normTraceOdd{clust_i}(cell_i,:) = zscore(smooth(meanTraceOdd{clust_i}(cell_i,:),3));
        normTraceOdd{clust_i}(cell_i,:) = normTraceOdd{clust_i}(cell_i,:) - min(normTraceOdd{clust_i}(cell_i,:));
        normTraceEven{clust_i}(cell_i,:) = zscore(smooth(meanTraceEven{clust_i}(cell_i,:),3));
        normTraceEven{clust_i}(cell_i,:) = normTraceEven{clust_i}(cell_i,:) - min(normTraceEven{clust_i}(cell_i,:));
        smoothSpeed{clust_i} = smooth(meanSpeed{clust_i},40);
        %find time index of max peaks
        [~,maxNormTrace{clust_i}(cell_i,1)] = max(normTrace{clust_i}(cell_i,:));
        [~,maxNormTraceOdd{clust_i}(cell_i,1)] = max(normTraceOdd{clust_i}(cell_i,:));
        
    end
    %sort each cell by the timing of its peak firing
    [B{clust_i},I{clust_i}] = sort(maxNormTrace{clust_i});
    [B1{clust_i},I1{clust_i}] = sort(maxNormTrace{1}); %sort by cluster 1
    [Bodd{clust_i},Iodd{clust_i}] = sort(maxNormTraceOdd{clust_i});
    %split dataset into even and odd for comparing 2 halves
    %     normTraceOdd{clust_i} = normTrace{clust_i}(1:2:end,:);
    %      normTraceEven{clust_i} = normTrace{clust_i}(2:2:end,:);
    %     %sort odd clusters in ascending order of the peak of the odds
    %     if length(normTraceOdd{clust_i}(:,1))>length(normTraceEven{clust_i}(:,1)) %if number of odd elements is greater than even
    %     [Bodd{clust_i},Iodd{clust_i}] = sort(maxNormTrace{clust_i}(1:2:end-1));
    %     else
    %     [Bodd{clust_i},Iodd{clust_i}] = sort(maxNormTrace{clust_i}(1:2:end));
    %     end
    
    %concatenate all clusters together and then zscore to equalize the
    %peaks across all trials
    meanTraceAll = [meanTraceAll meanTrace{clust_i}]; %concatenate all traces
    
end
%% this is to smooth, zscore, and sort the entire cell data by their preferred flight according to a homemade k-means (max dff across flights)
%zscore the full data set and subtract min to start at 0
for cell_ii = 1:length(cellData.results.C(:,1))
    normMeanTraceAll(cell_ii,:) = zscore(smooth(meanTraceAll(cell_ii,:),10));
    normMeanTraceAll(cell_ii,:) = normMeanTraceAll(cell_ii,:) - min(normMeanTraceAll(cell_ii,:));
end
%split the data back into its clusters
traceIndConcat = [1 length(normTrace{1}(1,:))]; %initialize index variable to hold the indices for start/stop of each cluster
%find the start and stop indices for 2nd through n cluster
for clust_ii = 2:clust_i
    traceIndConcat = vertcat(traceIndConcat,[traceIndConcat(clust_ii-1,2)+1 traceIndConcat(clust_ii-1,2)+length(normTrace{clust_ii}(1,:))]);
end
%regroup by each flight cluster
for clust_iii = 1:clust_i
    normTraceAll{clust_iii} = normMeanTraceAll(:,traceIndConcat(clust_iii,1):traceIndConcat(clust_iii,2));
    for cell_iii = 1:length(cellData.results.C(:,1))
        [~,maxNormAll{clust_iii}(cell_iii,1)] = max(normTraceAll{clust_iii}(cell_iii,:));
    end
    %sort each cell by the timing of its peak firing
    [BnormAll{clust_iii},InormAll{clust_iii}] = sort(maxNormAll{clust_iii});
    [B1normAll{clust_iii},I1normAll{clust_iii}] = sort(maxNormAll{1}); %sort by cluster 1
end


%find max and sort based on the peaks of the cell across the whole
%timeseries
[~,maxAll] = max(normMeanTraceAll,[],2);
%kPeaks = kmeans(maxAll,nClusters); %take k-means cluster of the peaks of each cell whole time
%sort based on each clustered peaks
rng(2);
for n = 1:size(maxAll,1);
    for ii = 1:nClusters
        if maxAll(n) >= traceIndConcat(ii,1) & maxAll(n) < traceIndConcat(ii,2); %if the peak is within each cluster, sort it into that particular cluster
            kPeaks(:,n) = 1;
        end
    end
end

%sort the clusters according to their peak times
[BkPeaks,IkPeaks] = sort(kPeaks);
normMeanTraceSort = normMeanTraceAll(IkPeaks,:);
%for each cluster, find the peaks of the cluster (aMat aka aTemp), tke the max
%(maxInd), and sort that max (aSort) from the cluster, add this to the
%normMeanTraceSort variable
for c = 1:nClusters
    aMat = find(BkPeaks == c);
    aTemp = normMeanTraceSort(aMat,:);
    [~,maxInd] = max(aTemp,[],2);
    [Bmax,Imax] = sort(maxInd);
    aSort = aTemp(Imax,:);
    normMeanTraceSort(aMat,:) = aSort;
    clear aMat aTemp maxInd Imax Bax aSort
end
%divide the normMeanTraceSort overall variable into the respective clusters
for c = 1:nClusters
    normMeanTraceEach{c} = normMeanTraceSort(:,traceIndConcat(c,1):traceIndConcat(c,2));
end

%% Plot the clusters normalized and clustered by their preferred flight (max dff across flights)
snakePlot_prefEach = figure();
for p = 1:nClusters
    
    p1 = subplot(4,clust_i,p);
    plot(1:length(smoothSpeed{p}),smoothSpeed{p},'k');
    hold on
    for cell_ii = 1:size(smoothSpeedRaw{p},1)
        plot(1:length(smoothSpeedRaw{p}(cell_ii,:)),smoothSpeedRaw{p}(cell_ii,:));
    end
    title(['Cluster ' num2str(p)]);
    if p == 1
        ylabel('velocity (cm/s)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
    else
        set(gca,'xticklabel',{[]},'yticklabel',{[]});
    end
    p2 = subplot(4,clust_i,[clust_i+p,2*clust_i+p,3*clust_i+p]);
    imagesc(normMeanTraceEach{p},[0.5 5.5]);
    colormap(hot);
    %make labels for first left plot only
    if p == 1
        ylabel('cell number');
    else
        set(gca,'yticklabel',{[]});
    end
    xt = get(gca, 'XTick');
    set(gca,'XTick',xt,'XTickLabel',round(xt/cellData.results.Fs,1));
    xlabel('time (s)');
    %hold off
    sgtitle(['Spatial selectivity sort by flight preference: ' batName ' ' dateSesh ' ' sessionType]);
end

%plot all the flight trajectories concatenated sorted by their preference
snakePlot_prefAll = figure();
imagesc(normMeanTraceSort);
colormap(hot);
ylabel('cell number');
set(gca,'yticklabel',{[]});
title(['Spatial selectivity sort by flight preference: ' batName ' ' dateSesh ' ' sessionType]);
xt = get(gca, 'XTick');
set(gca,'XTick',xt,'XTickLabel',round(xt/cellData.results.Fs,1));
xlabel('time (s)');


%%
%plot the cells according to their peak for each cluster with velocity on top
snakePlot_clust = figure();
for p = 1:clust_i
    p1 = subplot(4,clust_i,p);
    plot(1:length(smoothSpeed{p}),smoothSpeed{p},'k');
    hold on
    for cell_ii = 1:size(smoothSpeedRaw{p},1)
        plot(1:length(smoothSpeedRaw{p}(cell_ii,:)),smoothSpeedRaw{p}(cell_ii,:));
    end
    title(['Cluster ' num2str(p)]);
    if p == 1
        ylabel('velocity (cm/s)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
    else
        set(gca,'xticklabel',{[]},'yticklabel',{[]});
    end
    p2 = subplot(4,clust_i,[clust_i+p,2*clust_i+p,3*clust_i+p]);
    imagesc(normTrace{p}(I{p},:));
    colormap(hot);
    %make labels for first left plot only
    if p == 1
        ylabel('cell number');
    else
        set(gca,'yticklabel',{[]});
    end
    xt = get(gca, 'XTick');
    set(gca,'XTick',xt,'XTickLabel',round(xt/cellData.results.Fs,1));
    xlabel('time (s)');
    %hold off
    sgtitle(['Spatial selectivity sort by peak: ' batName ' ' dateSesh ' ' sessionType]);
    
end

%plot the cells by cluster sorted by cluster1 order
snakePlot_clustBy1 = figure();
for p = 1:clust_i
    p1 = subplot(4,clust_i,p);
    plot(1:length(smoothSpeed{p}),smoothSpeed{p},'k');
    hold on
    for cell_ii = 1:size(smoothSpeedRaw{p},1)
        plot(1:length(smoothSpeedRaw{p}(cell_ii,:)),smoothSpeedRaw{p}(cell_ii,:));
    end
    title(['Cluster ' num2str(p)]);
    if p == 1
        ylabel('velocity (cm/s)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
    else
        set(gca,'xticklabel',{[]},'yticklabel',{[]});
    end
    p2 = subplot(4,clust_i,[clust_i+p,2*clust_i+p,3*clust_i+p]);
    imagesc(normTrace{p}(I1{p},:));
    colormap(hot);
    %make labels for first left plot only
    if p == 1
        ylabel('cell number');
    else
        set(gca,'yticklabel',{[]});
    end
    xt = get(gca, 'XTick');
    set(gca,'XTick',xt,'XTickLabel',round(xt/cellData.results.Fs,1));
    xlabel('time (s)');
    %hold off
    sgtitle(['Spatial selectivity sort by cluster 1: ' batName ' ' dateSesh ' ' sessionType]);
    
end

%plot half the odd cells sorted and then plot the even cells according to
%the odd sorting
snakePlot_clustOddEven = figure();
for p = 1:clust_i
    p1 = subplot(4,2*clust_i,2*p-1);
    plot(1:length(smoothSpeed{p}),smoothSpeed{p},'k');
    hold on
    for cell_ii = 1:size(smoothSpeedRaw{p},1)
        plot(1:length(smoothSpeedRaw{p}(cell_ii,:)),smoothSpeedRaw{p}(cell_ii,:));
    end
    title(['Cluster ' num2str(p)]);
    if p == 1
        ylabel('velocity (cm/s)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
    else
        set(gca,'xticklabel',{[]},'yticklabel',{[]});
    end
    p2 = subplot(4,2*clust_i,[2*clust_i+(2*p-1),4*clust_i+(2*p-1),6*clust_i+(2*p-1)]);
    imagesc(normTraceOdd{p}(Iodd{p},:));
    colormap(hot);
    %make labels for first left plot only
    if p == 1
        ylabel('cell number');
    else
        set(gca,'yticklabel',{[]});
    end
    title(['Cluster ' num2str(p) ' Odd']);
    xt = get(gca, 'XTick');
    set(gca,'XTick',xt,'XTickLabel',round(xt/cellData.results.Fs,1));
    xlabel('time (s)'); %ylabel('cell number');
    p3 = subplot(4,2*clust_i,[2*clust_i+(2*p),4*clust_i+(2*p),6*clust_i+(2*p)]);
    imagesc(normTraceEven{p}(Iodd{p},:));
    colormap(hot);
    title(['Cluster ' num2str(p) ' Even']);
    xt = get(gca, 'XTick');
    set(gca,'XTick',xt,'XTickLabel',round(xt/cellData.results.Fs,1),'yticklabel',{[]});
    xlabel('time (s)'); %ylabel('cell number');
    sgtitle(['Spatial selectivity Odd (sorted) vs Even Trials: ' batName ' ' dateSesh ' ' sessionType]);
end

snakeTrace.meanTrace = meanTrace;
snakeTrace.normTrace = normTrace;
snakeTrace.maxNormTrace = maxNormTrace;
snakeTrace.sortedTrace = B;
snakeTrace.sortedIndex = I;
snakeTrace.sortedTraceOdd = Bodd;
snakeTrace.sortedIndexOdd = Iodd;
snakeTrace.sortedTraceClustBy1 = B1;
snakeTrace.sortedIndexClustBy1 = I1;
snakeTrace.normTraceEven = normTraceEven;
snakeTrace.normTraceOdd = normTraceOdd;
snakeTrace.meanSpeed = meanSpeed;
snakeTrace.smoothSpeed = smoothSpeed;
snakeTrace.snakePlot_clustAll = snakePlot_clust;
snakeTrace.snakePlot_clustOddEven = snakePlot_clustOddEven;
snakeTrace.snakePlot_clustBy1 = snakePlot_clustBy1;
snakeTrace.smoothSpeedRaw = smoothSpeedRaw;
snakeTrace.normMeanTraceEach = normMeanTraceEach;
snakeTrace.normMeanTraceSort = normMeanTraceSort;
snakeTrace.normMeanTraceAllSmooth = normMeanTraceAllSmooth;
snakeTrace.normMeanTraceAll = normMeanTraceAll;
snakeTrace.snakePlot_prefEach = snakePlot_prefEach;
snakeTrace.snakePlot_prefAll = snakePlot_prefAll;

%%
%plot the fully normalized cells according to their peak for each cluster with velocity on top
snakePlot_clust_allNorm = figure();
for p = 1:clust_i
    p1 = subplot(4,clust_i,p);
    plot(1:length(smoothSpeed{p}),smoothSpeed{p},'k');
    hold on
    for cell_ii = 1:size(smoothSpeedRaw{p},1)
        plot(1:length(smoothSpeedRaw{p}(cell_ii,:)),smoothSpeedRaw{p}(cell_ii,:));
    end
    title(['Cluster ' num2str(p)]);
    if p == 1
        ylabel('velocity (cm/s)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
    else
        set(gca,'xticklabel',{[]},'yticklabel',{[]});
    end
    p2 = subplot(4,clust_i,[clust_i+p,2*clust_i+p,3*clust_i+p]);
    imagesc(normTraceAll{p}(InormAll{p},:));
    colormap(hot);
    %make labels for first left plot only
    if p == 1
        ylabel('cell number');
    else
        set(gca,'yticklabel',{[]});
    end
    xt = get(gca, 'XTick');
    set(gca,'XTick',xt,'XTickLabel',round(xt/cellData.results.Fs,1));
    xlabel('time (s)');
    %hold off
    sgtitle(['Spatial selectivity (norm) sort by peak: ' batName ' ' dateSesh ' ' sessionType]);
    
end

%plot the cells by cluster sorted by cluster1 order
snakePlot_clustBy1_normAll = figure();
for p = 1:clust_i
    p1 = subplot(4,clust_i,p);
    plot(1:length(smoothSpeed{p}),smoothSpeed{p},'k');
    hold on
    for cell_ii = 1:size(smoothSpeedRaw{p},1)
        plot(1:length(smoothSpeedRaw{p}(cell_ii,:)),smoothSpeedRaw{p}(cell_ii,:));
    end
    title(['Cluster ' num2str(p)]);
    if p == 1
        ylabel('velocity (cm/s)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
    else
        set(gca,'xticklabel',{[]},'yticklabel',{[]});
    end
    p2 = subplot(4,clust_i,[clust_i+p,2*clust_i+p,3*clust_i+p]);
    imagesc(normTraceAll{p}(I1normAll{p},:));
    colormap(hot);
    %make labels for first left plot only
    if p == 1
        ylabel('cell number');
    else
        set(gca,'yticklabel',{[]});
    end
    xt = get(gca, 'XTick');
    set(gca,'XTick',xt,'XTickLabel',round(xt/cellData.results.Fs,1));
    xlabel('time (s)');
    %hold off
    sgtitle(['Spatial selectivity (norm) sort by cluster 1: ' batName ' ' dateSesh ' ' sessionType]);
    
end
%%
if saveFlag == 1
    saveas(snakeTrace.snakePlot_clustAll, [pwd '\analysis\snakePlots\' label '_snakePlots_clustAll.svg']);
    saveas(snakeTrace.snakePlot_clustOddEven, [pwd '\analysis\snakePlots\' label '_snakePlots_clustOddEven.svg']);
    saveas(snakeTrace.snakePlot_clustBy1, [pwd '\analysis\snakePlots\' label '_snakePlots_clustBy1.svg']);
    saveas(snakeTrace.snakePlot_prefEach, [pwd '\analysis\snakePlots\' label '_snakePlots_prefEach.svg']);
    saveas(snakeTrace.snakePlot_prefAll, [pwd '\analysis\snakePlots\' label '_snakePlots_prefAll.svg']);
    saveas(snakeTrace.snakePlot_clustAll, [pwd '\analysis\snakePlots\' label '_snakePlots_clustAll.tif']);
    saveas(snakeTrace.snakePlot_clustOddEven, [pwd '\analysis\snakePlots\' label '_snakePlots_clustOddEven.tif']);
    saveas(snakeTrace.snakePlot_clustBy1, [pwd '\analysis\snakePlots\' label '_snakePlots_clustBy1.tif']);
    saveas(snakeTrace.snakePlot_prefEach, [pwd '\analysis\snakePlots\' label '_snakePlots_prefEach.tif']);
    saveas(snakeTrace.snakePlot_prefAll, [pwd '\analysis\snakePlots\' label '_snakePlots_prefAll.tif']);
    save([pwd '/analysis/' label '_snakePlotData.mat'],'snakeTrace');
    savefig(snakeTrace.snakePlot_clustAll, [pwd '/analysis/snakePlots/' label '_snakePlots_clustAll.fig']);
    savefig(snakeTrace.snakePlot_clustOddEven, [pwd '/analysis/snakePlots/' label '_snakePlots_clustOddEven.fig']);
    savefig(snakeTrace.snakePlot_clustBy1, [pwd '/analysis/snakePlots/' label '_snakePlots_clustBy1.fig']);
    savefig(snakeTrace.snakePlot_prefEach, [pwd '/analysis/snakePlots/' label '_snakePlots_prefEach.fig']);
    savefig(snakeTrace.snakePlot_prefAll, [pwd '/analysis/snakePlots/' label '_snakePlots_prefAll.fig']);
end