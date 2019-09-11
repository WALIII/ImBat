function [snakeTrace] = ImBat_plotSnake(cellData,flightPaths,alignment,varargin)

%global batName dateSesh sessionType
prePad = 0.5*cellData.results.Fs; %number of frames to include in the trace extraction
postPad = 2*cellData.results.Fs; %add 2 seconds to the end of the plots to include delay in peak time
%meanTrace = cell(1,length(flightPaths.clusterIndex));

batName = [];
dateSesh = [];
sessionType = [];

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
    end
end

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
        for trace_i = 1:length(flightPaths.clusterIndex{clust_i})
            trace(trace_i,:) = cellData.results.C(cell_i,closestIndexStart(trace_i) - prePad:closestIndexEnd(trace_i) + (maxDur-dur(trace_i)) + postPad);
        end
        
        %calculate the mean neural activity across all flights in a cluster for each cell
        meanTrace{clust_i}(cell_i,:) = mean(trace);
        %smooth and zscore the neural data. subtract the min of the zscore so the
        %min is 0 rather than mean 0
        normTrace{clust_i}(cell_i,:) = zscore(smooth(meanTrace{clust_i}(cell_i,:),3));
        normTrace{clust_i}(cell_i,:) = normTrace{clust_i}(cell_i,:) - min(normTrace{clust_i}(cell_i,:));
        %find time index of 
        [~,maxNormTrace{clust_i}(cell_i,1)] = max(normTrace{clust_i}(cell_i,:));
    end
    
    %sort each cell by the timing of its peak firing
    [B{clust_i},I{clust_i}] = sort(maxNormTrace{clust_i});
    %split dataset into even and odd for comparing 2 halves 
    normTraceOdd{clust_i} = normTrace{clust_i}(1:2:end,:);
    normTraceEven{clust_i} = normTrace{clust_i}(2:2:end,:);
    %sort odd clusters in ascending order of the peak of the odds
    if length(normTraceOdd{clust_i}(:,1))>length(normTraceEven{clust_i}(:,1)) %if number of odd elements is greater than even
    [Bodd{clust_i},Iodd{clust_i}] = sort(maxNormTrace{clust_i}(1:2:end-1));
    else
    [Bodd{clust_i},Iodd{clust_i}] = sort(maxNormTrace{clust_i}(1:2:end));
    end

    %     %plot the cells according to their peak
    %     strcat('snakePlot_fig' num2str(clust_i)) = figure();
    %     imagesc(normTrace{clust_i}(I{clust_i},:));
    %     colormap(hot);
    %     hold on
    %     %set(gcf,'Name',figName, 'NumberTitle','Off');
    %     title(['Spatial selectivity cluster ' num2str(clust_i) ' : '  batName ' ' dateSesh ' ' sessionType]);
    %     xt = get(gca, 'XTick');
    %     set(gca,'XTick',xt,'XTickLabel',round(xt/cellData.results.Fs,1));
%     xlabel('time (s)'); ylabel('cell number');
%     hold off
end
snakeTrace.meanTrace = meanTrace;
snakeTrace.normTrace = normTrace;
snakeTrace.maxNormTrace = maxNormTrace;
snakeTrace.sortedTrace = B;
snakeTrace.sortedIndex = I;
snakeTrace.sortedTraceOdd = Bodd;
snakeTrace.sortedIndexOdd = Iodd;
snakeTrace.normTraceEven = normTraceEven;
snakeTrace.normTraceOdd = normTraceOdd;

%%
%plot the cells according to their peak for each cluster
figure();
snakeTrace.snakePlot_fig1 = imagesc(snakeTrace.normTrace{1}(snakeTrace.sortedIndex{1},:));
colormap(hot);
title(['Spatial selectivity cluster 1: ' batName ' ' dateSesh ' ' sessionType]);
xt = get(gca, 'XTick');
set(gca,'XTick',xt,'XTickLabel',round(xt/cellData.results.Fs,1));
xlabel('time (s)'); ylabel('cell number');
%hold off

figure();
snakeTrace.snakePlot_fig2 = imagesc(snakeTrace.normTrace{2}(snakeTrace.sortedIndex{2},:));
colormap(hot);
title(['Spatial selectivity cluster 2: ' batName ' ' dateSesh ' ' sessionType]);
xt = get(gca, 'XTick');
set(gca,'XTick',xt,'XTickLabel',round(xt/cellData.results.Fs,1));
xlabel('time (s)'); ylabel('cell number');
%hold off

figure();
snakeTrace.snakePlot_fig3 = imagesc(snakeTrace.normTrace{3}(snakeTrace.sortedIndex{3},:));
colormap(hot);
title(['Spatial selectivity cluster 3: ' batName ' ' dateSesh ' ' sessionType]);
xt = get(gca, 'XTick');
set(gca,'XTick',xt,'XTickLabel',round(xt/cellData.results.Fs,1));
xlabel('time (s)'); ylabel('cell number');
%hold off

figure();
snakeTrace.snakePlot_fig4 = imagesc(snakeTrace.normTrace{4}(snakeTrace.sortedIndex{4},:));
colormap(hot);
title(['Spatial selectivity cluster 4: ' batName ' ' dateSesh ' ' sessionType]);
xt = get(gca, 'XTick');
set(gca,'XTick',xt,'XTickLabel',round(xt/cellData.results.Fs,1));
xlabel('time (s)'); ylabel('cell number');
%hold off

%%
figure();
snakeTrace.snakePlot_fig1_odd = imagesc(snakeTrace.normTraceOdd{1}(snakeTrace.sortedIndexOdd{1},:));
colormap(hot);
title(['Spatial selectivity cluster 1 Odd: ' batName ' ' dateSesh ' ' sessionType]);
xt = get(gca, 'XTick');
set(gca,'XTick',xt,'XTickLabel',round(xt/cellData.results.Fs,1));
xlabel('time (s)'); ylabel('cell number');
%hold off

figure();
snakeTrace.snakePlot_fig1_even = imagesc(snakeTrace.normTraceEven{1}(snakeTrace.sortedIndexOdd{1},:));
colormap(hot);
title(['Spatial selectivity cluster 1 Even: ' batName ' ' dateSesh ' ' sessionType]);
xt = get(gca, 'XTick');
set(gca,'XTick',xt,'XTickLabel',round(xt/cellData.results.Fs,1));
xlabel('time (s)'); ylabel('cell number');

figure();
snakeTrace.snakePlot_fig2_odd = imagesc(snakeTrace.normTraceOdd{2}(snakeTrace.sortedIndexOdd{2},:));
colormap(hot);
title(['Spatial selectivity cluster 2 Odd: ' batName ' ' dateSesh ' ' sessionType]);
xt = get(gca, 'XTick');
set(gca,'XTick',xt,'XTickLabel',round(xt/cellData.results.Fs,1));
xlabel('time (s)'); ylabel('cell number');
%hold off

figure();
snakeTrace.snakePlot_fig2_even = imagesc(snakeTrace.normTraceEven{2}(snakeTrace.sortedIndexOdd{2},:));
colormap(hot);
title(['Spatial selectivity cluster 2 Even: ' batName ' ' dateSesh ' ' sessionType]);
xt = get(gca, 'XTick');
set(gca,'XTick',xt,'XTickLabel',round(xt/cellData.results.Fs,1));
xlabel('time (s)'); ylabel('cell number');

figure();
snakeTrace.snakePlot_fig3_odd = imagesc(snakeTrace.normTraceOdd{3}(snakeTrace.sortedIndexOdd{3},:));
colormap(hot);
title(['Spatial selectivity cluster 3 Odd: ' batName ' ' dateSesh ' ' sessionType]);
xt = get(gca, 'XTick');
set(gca,'XTick',xt,'XTickLabel',round(xt/cellData.results.Fs,1));
xlabel('time (s)'); ylabel('cell number');
%hold off

figure();
snakeTrace.snakePlot_fig3_even = imagesc(snakeTrace.normTraceEven{3}(snakeTrace.sortedIndexOdd{3},:));
colormap(hot);
title(['Spatial selectivity cluster 3 Even: ' batName ' ' dateSesh ' ' sessionType]);
xt = get(gca, 'XTick');
set(gca,'XTick',xt,'XTickLabel',round(xt/cellData.results.Fs,1));
xlabel('time (s)'); ylabel('cell number');

figure();
snakeTrace.snakePlot_fig4_odd = imagesc(snakeTrace.normTraceOdd{4}(snakeTrace.sortedIndexOdd{4},:));
colormap(hot);
title(['Spatial selectivity cluster 4 Odd: ' batName ' ' dateSesh ' ' sessionType]);
xt = get(gca, 'XTick');
set(gca,'XTick',xt,'XTickLabel',round(xt/cellData.results.Fs,1));
xlabel('time (s)'); ylabel('cell number');
%hold off

figure();
snakeTrace.snakePlot_fig4_even = imagesc(snakeTrace.normTraceEven{4}(snakeTrace.sortedIndexOdd{4},:));
colormap(hot);
title(['Spatial selectivity cluster 4 Even: ' batName ' ' dateSesh ' ' sessionType]);
xt = get(gca, 'XTick');
set(gca,'XTick',xt,'XTickLabel',round(xt/cellData.results.Fs,1));
xlabel('time (s)'); ylabel('cell number');
