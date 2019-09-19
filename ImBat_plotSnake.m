function [snakeTrace] = ImBat_plotSnake(cellData,flightPaths,alignment,varargin)
%global batName dateSesh sessionType


batName = [];
dateSesh = [];
sessionType = [];
saveFlag = 0; %do you want to load and save the data individually outside of ImBatAnalyze


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
        case 'saveflag'
            saveFlag = varargin{i+1};
    end
end

%labels for loading and saving data if running independent fromImBat_analyze
if saveFlag == 1
date = strcat(lower(batName(1:2)),dateSesh);
label = [batName '_' dateSesh '_' sessionType];

cellData = load([pwd '/processed/Motion_corrected_Data_DS_results.mat']);
alignment = load([pwd '/processed/Alignment.mat']);
load([pwd '/analysis/' label '_flightPaths.mat']);
end

prePad = 2*cellData.results.Fs; %number of frames (seconds*freq) to include in the trace extraction
postPad = 6*cellData.results.Fs; %add 2 seconds to the end of the plots to include delay in peak time
%meanTrace = cell(1,length(flightPaths.clusterIndex));

%for each cluster type
for clust_i = 1:4 %length(flightPaths.clusterIndex)
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
        end
        %calculate max duration for each cluster of trajectories
        maxDur = max(dur);
        %meanTrace{clust_i}=zeros(length(cellData.results.C(:,1)),maxDur+preWindow+1);
        %initialize the vector to store the neural activity of each flight
        trace = zeros(length(flightPaths.clusterIndex{clust_i}),maxDur+prePad+postPad+1);
        speed = zeros(length(flightPaths.clusterIndex{clust_i}),maxDur+prePad+postPad+1);
        for trace_i = 1:length(flightPaths.clusterIndex{clust_i})
            try
            trace(trace_i,:) = cellData.results.C(cell_i,closestIndexStart(trace_i) - prePad:closestIndexEnd(trace_i) + (maxDur-dur(trace_i)) + postPad);
            speed(trace_i,:) = flightPaths.batSpeed(closestIndexStart(trace_i) - prePad:closestIndexEnd(trace_i) + (maxDur-dur(trace_i)) + postPad); 
            catch
                sizeToRecordingEnd = size(cellData.results.C(cell_i,closestIndexStart(trace_i) - prePad:end),2);
                sizeToTraceEnd = size(trace(trace_i,:),2);
                trace(trace_i,:) = (cellData.results.C(cell_i,closestIndexStart(trace_i) - prePad:end)+(zeros(1,sizeToTraceEnd - sizeToRecordingEnd)));
                speed(trace_i,:) = (flightPaths.batSpeed(closestIndexStart(trace_i) - prePad:end)+(zeros(1,sizeToTraceEnd - sizeToRecordingEnd))); 
            
            end    
        end
        
        %calculate the mean neural activity across all flights in a cluster for each cell
        meanTrace{clust_i}(cell_i,:) = mean(trace);
        meanTraceOdd{clust_i}(cell_i,:) = mean(trace(1:2:end,:));
        meanTraceEven{clust_i}(cell_i,:) = mean(trace(2:2:end,:));
        meanSpeed{clust_i} = mean(speed);
        %smooth and zscore the neural data. subtract the min of the zscore so the
        %min is 0 rather than mean 0
        normTrace{clust_i}(cell_i,:) = zscore(smooth(meanTrace{clust_i}(cell_i,:),3));
        normTrace{clust_i}(cell_i,:) = normTrace{clust_i}(cell_i,:) - min(normTrace{clust_i}(cell_i,:));
        normTraceOdd{clust_i}(cell_i,:) = zscore(smooth(meanTraceOdd{clust_i}(cell_i,:),3));
        normTraceOdd{clust_i}(cell_i,:) = normTraceOdd{clust_i}(cell_i,:) - min(normTraceOdd{clust_i}(cell_i,:));
        normTraceEven{clust_i}(cell_i,:) = zscore(smooth(meanTraceEven{clust_i}(cell_i,:),3));
        normTraceEven{clust_i}(cell_i,:) = normTraceEven{clust_i}(cell_i,:) - min(normTraceEven{clust_i}(cell_i,:));
        smoothSpeed{clust_i} = smooth(meanTrace{clust_i}(cell_i,:),20);
        %find time index of max peaks
        [~,maxNormTrace{clust_i}(cell_i,1)] = max(normTrace{clust_i}(cell_i,:));
        [~,maxNormTraceOdd{clust_i}(cell_i,1)] = max(normTraceOdd{clust_i}(cell_i,:));

    end
    
    %sort each cell by the timing of its peak firing
    [B{clust_i},I{clust_i}] = sort(maxNormTrace{clust_i});
    [B1{clust_i},I1{clust_i}] = sort(maxNormTrace{1});
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
end



%%
%plot the cells according to their peak for each cluster with velocity on top
snakePlot_clust = figure();
for p = 1:clust_i
    p1 = subplot(4,clust_i,p);
    plot(1:length(smoothSpeed{p}),smoothSpeed{p});
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
    plot(1:length(smoothSpeed{p}),smoothSpeed{p});
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
    plot(1:length(smoothSpeed{p}),smoothSpeed{p});
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

%%
if saveFlag == 1
saveas(snakeTrace.snakePlot_clustAll, [pwd '\analysis\snakePlots\' label '_snakePlots_clustAll.svg']);
saveas(snakeTrace.snakePlot_clustOddEven, [pwd '\analysis\snakePlots\' label '_snakePlots_clustOddEven.svg']);
saveas(snakeTrace.snakePlot_clustBy1, [pwd '\analysis\snakePlots\' label '_snakePlots_clustBy1.svg']);
saveas(snakeTrace.snakePlot_clustAll, [pwd '\analysis\snakePlots\' label '_snakePlots_clustAll.tif']);
saveas(snakeTrace.snakePlot_clustOddEven, [pwd '\analysis\snakePlots\' label '_snakePlots_clustOddEven.tif']);
saveas(snakeTrace.snakePlot_clustBy1, [pwd '\analysis\snakePlots\' label '_snakePlots_clustBy1.tif']);
save([pwd '/analysis/' label '_snakePlotData.mat'],'snakeTrace');
savefig(snakeTrace.snakePlot_clustAll, [pwd '/analysis/snakePlots/' label '_snakePlots_clustAll.fig']);
savefig(snakeTrace.snakePlot_clustOddEven, [pwd '/analysis/snakePlots/' label '_snakePlots_clustOddEven.fig']);
savefig(snakeTrace.snakePlot_clustBy1, [pwd '/analysis/snakePlots/' label '_snakePlots_clustBy1.fig']);
end