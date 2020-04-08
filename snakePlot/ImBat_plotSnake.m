function [snakeTrace] = ImBat_plotSnake(snakeTrace,varargin)
% User inputs overrides
saveFlag = 0;

nparams=length(varargin);
for i=1:2:nparams
    switch lower(varargin{i})
        case 'saveflag'
            saveFlag = varargin{i+1};
        case 'analysisfolder'
            analysis_Folder = varargin{i+1};
    end
end

%%
% plot clusters normalized and clustered by preferred flight/pre/post (max dff across flights)
snakePlot_prefEachPrePostFlight = figure();
for p = 1:snakeTrace.nClusters
    
    p1 = subplot(12,snakeTrace.nClusters,p);
    plot(1:length(snakeTrace.smoothSpeedPre{p}),snakeTrace.smoothSpeedPre{p},'k');
    hold on
    for cell_ii = 1:size(snakeTrace.smoothSpeedRawPre{p},1)
        plot(1:length(snakeTrace.smoothSpeedRawPre{p}(cell_ii,:)),snakeTrace.smoothSpeedRawPre{p}(cell_ii,:));
    end
    title(['Cluster ' num2str(p)]);
    if p == 1
        ylabel('velocity (cm/s)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
    else
        %set(gca,'xticklabel',{[]},'yticklabel',{[]});
    end
    
    p2 = subplot(12,snakeTrace.nClusters,[snakeTrace.nClusters+p,2*snakeTrace.nClusters+p,3*snakeTrace.nClusters+p]);
    imagesc(snakeTrace.normMeanTraceEachPre{p},[0.5 5.5]);
    colormap(hot);
    %make labels for first left plot only
    if p == 1
        %ylabel('cell number');
        ylabel('pre-flight');
    else
        set(gca,'yticklabel',{[]});
    end
    xt = get(gca, 'XTick');
    %set(gca,'XTick',xt,'XTickLabel',round(xt/cellData.results.Fs,1));
    xlabel('time (s)');
    %hold off
    
    p3 = subplot(12,snakeTrace.nClusters,p+20);
    plot(1:length(snakeTrace.smoothSpeedFlight{p}),snakeTrace.smoothSpeedFlight{p},'k');
    hold on
    for cell_ii = 1:size(snakeTrace.smoothSpeedRawFlight{p},1)
        plot(1:length(snakeTrace.smoothSpeedRawFlight{p}(cell_ii,:)),snakeTrace.smoothSpeedRawFlight{p}(cell_ii,:));
    end
    
    p4 = subplot(12,snakeTrace.nClusters,[snakeTrace.nClusters+p+20,2*snakeTrace.nClusters+p+20,3*snakeTrace.nClusters+p+20]);
    imagesc(snakeTrace.normMeanTraceEachFlight{p},[0.5 5.5]);
    colormap(hot);
    if p == 1
        %ylabel('cell number');
        ylabel('during flight');
    else
        set(gca,'yticklabel',{[]});
    end
    xt = get(gca, 'XTick');
   
    p5 = subplot(12,snakeTrace.nClusters,p+40);
    plot(1:length(snakeTrace.smoothSpeedPost{p}),snakeTrace.smoothSpeedPost{p},'k');
    hold on
    for cell_ii = 1:size(snakeTrace.smoothSpeedRawPost{p},1)
        plot(1:length(snakeTrace.smoothSpeedRawPost{p}(cell_ii,:)),snakeTrace.smoothSpeedRawPost{p}(cell_ii,:));
    end
    
    p6 = subplot(12,snakeTrace.nClusters,[snakeTrace.nClusters+p+40,2*snakeTrace.nClusters+p+40,3*snakeTrace.nClusters+p+40]);
    imagesc(snakeTrace.normMeanTraceEachPost{p},[0.5 5.5]);
    colormap(hot);
    if p == 1
        %ylabel('cell number');
        ylabel('post-flight');
    else
        set(gca,'yticklabel',{[]});
    end
    xt = get(gca, 'XTick');
    
    
    sgtitle(['Selectivity sort by pre/post/flight preference: ' snakeTrace.batName ' ' snakeTrace.dateSesh ' ' snakeTrace.sessionType]);
end

%plot all the flight trajectories concatenated sorted by their preference
snakePlot_prefAll = figure();
subplot(3,1,1)
imagesc(snakeTrace.normMeanTraceSortPre);
colormap(hot);
ylabel('pre-flight');
set(gca,'yticklabel',{[]});
set(gca,'xticklabel',{[]});
hold on

subplot(3,1,2)
imagesc(snakeTrace.normMeanTraceSortFlight);
colormap(hot);
ylabel('during flight');
set(gca,'yticklabel',{[]});
set(gca,'xticklabel',{[]});

subplot(3,1,3)
imagesc(snakeTrace.normMeanTraceSortPost);
colormap(hot);
ylabel('post-flight');
set(gca,'yticklabel',{[]});
xt = get(gca, 'XTick');
set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
xlabel('time (s)');
sgtitle(['Spatial selectivity sort by flight preference (pre/post/flight): ' snakeTrace.batName ' ' snakeTrace.dateSesh ' ' snakeTrace.sessionType]);



%% Plot the FLIGHT clusters normalized and clustered by their preferred flight (max dff across flights)
snakePlot_prefEach = figure();
for p = 1:snakeTrace.nClusters
    
    p1 = subplot(4,snakeTrace.nClusters,p);
    plot(1:length(snakeTrace.smoothSpeedFlight{p}),snakeTrace.smoothSpeedFlight{p},'k');
    hold on
    for cell_ii = 1:size(snakeTrace.smoothSpeedRawFlight{p},1)
        plot(1:length(snakeTrace.smoothSpeedRawFlight{p}(cell_ii,:)),snakeTrace.smoothSpeedRawFlight{p}(cell_ii,:));
    end
    title(['Cluster ' num2str(p)]);
    if p == 1
        ylabel('velocity (cm/s)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
    else
        %set(gca,'xticklabel',{[]},'yticklabel',{[]});
    end
    p2 = subplot(4,snakeTrace.nClusters,[snakeTrace.nClusters+p,2*snakeTrace.nClusters+p,3*snakeTrace.nClusters+p]);
    imagesc(snakeTrace.normMeanTraceEachFlight{p},[0.5 5.5]);
    colormap(hot);
    %make labels for first left plot only
    if p == 1
        ylabel('cell number');
    else
        set(gca,'yticklabel',{[]});
    end
    xt = get(gca, 'XTick');
    %set(gca,'XTick',xt,'XTickLabel',round(xt/cellData.results.Fs,1));
    xlabel('time (s)');
    %hold off
    sgtitle(['Spatial selectivity sort by flight preference: ' snakeTrace.batName ' ' snakeTrace.dateSesh ' ' snakeTrace.sessionType]);
end

%plot all the flight trajectories concatenated sorted by their preference
snakePlot_prefAll = figure();
imagesc(snakeTrace.normMeanTraceSortFlight);
colormap(hot);
ylabel('cell number');
set(gca,'yticklabel',{[]});
title(['Spatial selectivity sort by flight preference: ' snakeTrace.batName ' ' snakeTrace.dateSesh ' ' snakeTrace.sessionType]);
xt = get(gca, 'XTick');
set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
xlabel('time (s)');
%% plot the cells according to their peak for each cluster with velocity on top (pre/post/flight)
snakePlot_clustPrePostFlight = figure();
for p = 1:snakeTrace.nClusters
    p1 = subplot(12,snakeTrace.nClusters,p);
    plot(1:length(snakeTrace.smoothSpeedPre{p}),snakeTrace.smoothSpeedPre{p},'k');
    hold on
    for cell_ii = 1:size(snakeTrace.smoothSpeedRawFlight{p},1)
        plot(1:length(snakeTrace.smoothSpeedRawPre{p}(cell_ii,:)),snakeTrace.smoothSpeedRawPre{p}(cell_ii,:));
    end
    title(['Cluster ' num2str(p)]);
    if p == 1
        ylabel('velocity (cm/s)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
    else
        set(gca,'xticklabel',{[]},'yticklabel',{[]});
    end
    p2 = subplot(12,snakeTrace.nClusters,[snakeTrace.nClusters+p,2*snakeTrace.nClusters+p,3*snakeTrace.nClusters+p]);
    imagesc(snakeTrace.normTracePre{p}(snakeTrace.IPre{p},:));
    colormap(hot);
    %make labels for first left plot only
    if p == 1
        ylabel('pre-flight');
    else
        set(gca,'yticklabel',{[]});
    end
    set(gca,'xticklabel',[]);
    hold on
    
    p3 = subplot(12,snakeTrace.nClusters,p+20);
    plot(1:length(snakeTrace.smoothSpeedFlight{p}),snakeTrace.smoothSpeedFlight{p},'k');
    hold on
    for cell_ii = 1:size(snakeTrace.smoothSpeedRawFlight{p},1)
        plot(1:length(snakeTrace.smoothSpeedRawFlight{p}(cell_ii,:)),snakeTrace.smoothSpeedRawFlight{p}(cell_ii,:));
    end
    title(['Cluster ' num2str(p)]);
    if p == 1
        ylabel('velocity (cm/s)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
    else
        set(gca,'xticklabel',{[]},'yticklabel',{[]});
    end
    p4 = subplot(12,snakeTrace.nClusters,[snakeTrace.nClusters+p+20,2*snakeTrace.nClusters+p+20,3*snakeTrace.nClusters+p+20]);
    imagesc(snakeTrace.normTraceFlight{p}(snakeTrace.IFlight{p},:));
    colormap(hot);
    %make labels for first left plot only
    if p == 1
        ylabel('during flight');
    else
        set(gca,'yticklabel',{[]});
    end
    set(gca,'xticklabel',[]);

    
    p5 = subplot(12,snakeTrace.nClusters,p+40);
    plot(1:length(snakeTrace.smoothSpeedPost{p}),snakeTrace.smoothSpeedPost{p},'k');
    hold on
    for cell_ii = 1:size(snakeTrace.smoothSpeedRawPost{p},1)
        plot(1:length(snakeTrace.smoothSpeedRawPost{p}(cell_ii,:)),snakeTrace.smoothSpeedRawPost{p}(cell_ii,:));
    end
    title(['Cluster ' num2str(p)]);
    if p == 1
        ylabel('velocity (cm/s)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
    else
        set(gca,'xticklabel',{[]},'yticklabel',{[]});
    end
    p6 = subplot(12,snakeTrace.nClusters,[snakeTrace.nClusters+p+40,2*snakeTrace.nClusters+p+40,3*snakeTrace.nClusters+p+40]);
    imagesc(snakeTrace.normTracePost{p}(snakeTrace.IPost{p},:));
    colormap(hot);
    %make labels for first left plot only
    if p == 1
        ylabel('post-flight');
    else
        set(gca,'yticklabel',{[]});
    end
    xt = get(gca, 'XTick');
    set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
    xlabel('time (s)');
    
    sgtitle(['Spatial selectivity sort by peak (pre/post/flight): ' snakeTrace.batName ' ' snakeTrace.dateSesh ' ' snakeTrace.sessionType]);
    
end

%plot the cells by cluster sorted by cluster1 order
snakePlot_clustBy1PrePostFlight = figure();
for p = 1:snakeTrace.nClusters
    p1 = subplot(12,snakeTrace.nClusters,p);
    plot(1:length(snakeTrace.smoothSpeedPre{p}),snakeTrace.smoothSpeedPre{p},'k');
    hold on
    for cell_ii = 1:size(snakeTrace.smoothSpeedRawPre{p},1)
        plot(1:length(snakeTrace.smoothSpeedRawPre{p}(cell_ii,:)),snakeTrace.smoothSpeedRawPre{p}(cell_ii,:));
    end
    title(['Cluster ' num2str(p)]);
    if p == 1
        ylabel('velocity (cm/s)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
    else
        set(gca,'xticklabel',{[]},'yticklabel',{[]});
    end
    p2 = subplot(12,snakeTrace.nClusters,[snakeTrace.nClusters+p,2*snakeTrace.nClusters+p,3*snakeTrace.nClusters+p]);
    imagesc(snakeTrace.normTracePre{p}(snakeTrace.I1Pre{p},:));
    colormap(hot);
    %make labels for first left plot only
    if p == 1
        ylabel('pre-flight');
    else
        set(gca,'yticklabel',{[]});
    end
    set(gca,'xticklabel',{[]});
    
    p3 = subplot(12,snakeTrace.nClusters,p+20);
    plot(1:length(snakeTrace.smoothSpeedFlight{p}),snakeTrace.smoothSpeedFlight{p},'k');
    hold on
    for cell_ii = 1:size(snakeTrace.smoothSpeedRawFlight{p},1)
        plot(1:length(snakeTrace.smoothSpeedRawFlight{p}(cell_ii,:)),snakeTrace.smoothSpeedRawFlight{p}(cell_ii,:));
    end
    title(['Cluster ' num2str(p)]);
    if p == 1
        ylabel('velocity (cm/s)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
    else
        set(gca,'xticklabel',{[]},'yticklabel',{[]});
    end
    p4 = subplot(12,snakeTrace.nClusters,[snakeTrace.nClusters+p+20,2*snakeTrace.nClusters+p+20,3*snakeTrace.nClusters+p+20]);
    imagesc(snakeTrace.normTraceFlight{p}(snakeTrace.I1Flight{p},:));
    colormap(hot);
    %make labels for first left plot only
    if p == 1
        ylabel('during flight');
    else
        set(gca,'yticklabel',{[]});
    end
    set(gca,'xticklabel',{[]});

    p5 = subplot(12,snakeTrace.nClusters,p+40);
    plot(1:length(snakeTrace.smoothSpeedPost{p}),snakeTrace.smoothSpeedPost{p},'k');
    hold on
    for cell_ii = 1:size(snakeTrace.smoothSpeedRawPost{p},1)
        plot(1:length(snakeTrace.smoothSpeedRawPost{p}(cell_ii,:)),snakeTrace.smoothSpeedRawPost{p}(cell_ii,:));
    end
    title(['Cluster ' num2str(p)]);
    if p == 1
        ylabel('velocity (cm/s)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
    else
        set(gca,'xticklabel',{[]},'yticklabel',{[]});
    end
    p6 = subplot(12,snakeTrace.nClusters,[snakeTrace.nClusters+p+40,2*snakeTrace.nClusters+p+40,3*snakeTrace.nClusters+p+40]);
    imagesc(snakeTrace.normTracePost{p}(snakeTrace.I1Post{p},:));
    colormap(hot);
    %make labels for first left plot only
    if p == 1
        ylabel('post-flight');
    else
        set(gca,'yticklabel',{[]});
    end
    xt = get(gca, 'XTick');
    set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
    xlabel('time (s)');
    sgtitle(['Spatial selectivity sort by cluster 1 Pre/post/flight: ' snakeTrace.batName ' ' snakeTrace.dateSesh ' ' snakeTrace.sessionType]);
    
end

%%
%plot the cells according to their peak for each cluster with velocity on top
snakePlot_clust = figure();
for p = 1:snakeTrace.nClusters
    p1 = subplot(4,snakeTrace.nClusters,p);
    plot(1:length(snakeTrace.smoothSpeedFlight{p}),snakeTrace.smoothSpeedFlight{p},'k');
    hold on
    for cell_ii = 1:size(snakeTrace.smoothSpeedRawFlight{p},1)
        plot(1:length(snakeTrace.smoothSpeedRawFlight{p}(cell_ii,:)),snakeTrace.smoothSpeedRawFlight{p}(cell_ii,:));
    end
    title(['Cluster ' num2str(p)]);
    if p == 1
        ylabel('velocity (cm/s)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
    else
        set(gca,'xticklabel',{[]},'yticklabel',{[]});
    end
    p2 = subplot(4,snakeTrace.nClusters,[snakeTrace.nClusters+p,2*snakeTrace.nClusters+p,3*snakeTrace.nClusters+p]);
    imagesc(snakeTrace.normTraceFlight{p}(snakeTrace.IFlight{p},:));
    colormap(hot);
    %make labels for first left plot only
    if p == 1
        ylabel('cell number');
    else
        set(gca,'yticklabel',{[]});
    end
    xt = get(gca, 'XTick');
    set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
    xlabel('time (s)');
    %hold off
    sgtitle(['Spatial selectivity sort by peak: ' snakeTrace.batName ' ' snakeTrace.dateSesh ' ' snakeTrace.sessionType]);
    
end

%plot the cells by cluster sorted by cluster1 order
snakePlot_clustBy1 = figure();
for p = 1:snakeTrace.nClusters
    p1 = subplot(4,snakeTrace.nClusters,p);
    plot(1:length(snakeTrace.smoothSpeedFlight{p}),snakeTrace.smoothSpeedFlight{p},'k');
    hold on
    for cell_ii = 1:size(snakeTrace.smoothSpeedRawFlight{p},1)
        plot(1:length(snakeTrace.smoothSpeedRawFlight{p}(cell_ii,:)),snakeTrace.smoothSpeedRawFlight{p}(cell_ii,:));
    end
    title(['Cluster ' num2str(p)]);
    if p == 1
        ylabel('velocity (cm/s)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
    else
        set(gca,'xticklabel',{[]},'yticklabel',{[]});
    end
    p2 = subplot(4,snakeTrace.nClusters,[snakeTrace.nClusters+p,2*snakeTrace.nClusters+p,3*snakeTrace.nClusters+p]);
    imagesc(snakeTrace.normTraceFlight{p}(snakeTrace.I1Flight{p},:));
    colormap(hot);
    %make labels for first left plot only
    if p == 1
        ylabel('cell number');
    else
        set(gca,'yticklabel',{[]});
    end
    xt = get(gca, 'XTick');
    set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
    xlabel('time (s)');
    %hold off
    sgtitle(['Spatial selectivity sort by cluster 1: ' snakeTrace.batName ' ' snakeTrace.dateSesh ' ' snakeTrace.sessionType]);
    
end

%plot half the odd cells sorted and then plot the even cells according to
%the odd sorting
snakePlot_clustOddEven = figure();
for p = 1:snakeTrace.nClusters
    p1 = subplot(4,2*snakeTrace.nClusters,2*p-1);
    plot(1:length(snakeTrace.smoothSpeedFlight{p}),snakeTrace.smoothSpeedFlight{p},'k');
    hold on
    for cell_ii = 1:size(snakeTrace.smoothSpeedRawFlight{p},1)
        plot(1:length(snakeTrace.smoothSpeedRawFlight{p}(cell_ii,:)),snakeTrace.smoothSpeedRawFlight{p}(cell_ii,:));
    end
    title(['Cluster ' num2str(p)]);
    if p == 1
        ylabel('velocity (cm/s)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
    else
        set(gca,'xticklabel',{[]},'yticklabel',{[]});
    end
    p2 = subplot(4,2*snakeTrace.nClusters,[2*snakeTrace.nClusters+(2*p-1),4*snakeTrace.nClusters+(2*p-1),6*snakeTrace.nClusters+(2*p-1)]);
    imagesc(snakeTrace.normTraceOdd{p}(snakeTrace.Iodd{p},:));
    colormap(hot);
    %make labels for first left plot only
    if p == 1
        ylabel('cell number');
    else
        set(gca,'yticklabel',{[]});
    end
    title(['Cluster ' num2str(p) ' Odd']);
    xt = get(gca, 'XTick');
    set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
    xlabel('time (s)'); %ylabel('cell number');
    p3 = subplot(4,2*snakeTrace.nClusters,[2*snakeTrace.nClusters+(2*p),4*snakeTrace.nClusters+(2*p),6*snakeTrace.nClusters+(2*p)]);
    imagesc(snakeTrace.normTraceEven{p}(snakeTrace.Iodd{p},:));
    colormap(hot);
    title(['Cluster ' num2str(p) ' Even']);
    xt = get(gca, 'XTick');
    set(gca,'XTick',xt,'XTickLabel',round(xt/30,1),'yticklabel',{[]});
    xlabel('time (s)'); %ylabel('cell number');
    sgtitle(['Spatial selectivity Odd (sorted) vs Even Trials: ' snakeTrace.batName ' ' snakeTrace.dateSesh ' ' snakeTrace.sessionType]);
end


%%
%plot the fully normalized cells according to their peak for each cluster with velocity on top
snakePlot_clust_allNorm = figure();
for p = 1:snakeTrace.nClusters
    p1 = subplot(4,snakeTrace.nClusters,p);
    plot(1:length(snakeTrace.smoothSpeedFlight{p}),snakeTrace.smoothSpeedFlight{p},'k');
    hold on
    for cell_ii = 1:size(snakeTrace.smoothSpeedRawFlight{p},1)
        plot(1:length(snakeTrace.smoothSpeedRawFlight{p}(cell_ii,:)),snakeTrace.smoothSpeedRawFlight{p}(cell_ii,:));
    end
    title(['Cluster ' num2str(p)]);
    if p == 1
        ylabel('velocity (cm/s)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
    else
        set(gca,'xticklabel',{[]},'yticklabel',{[]});
    end
    p2 = subplot(4,snakeTrace.nClusters,[snakeTrace.nClusters+p,2*snakeTrace.nClusters+p,3*snakeTrace.nClusters+p]);
    imagesc(snakeTrace.normTraceFlightAll{p}(snakeTrace.InormFlightAll{p},:));
    colormap(hot);
    %make labels for first left plot only
    if p == 1
        ylabel('cell number');
    else
        set(gca,'yticklabel',{[]});
    end
    xt = get(gca, 'XTick');
    set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
    xlabel('time (s)');
    %hold off
    sgtitle(['Spatial selectivity (norm) sort by peak: ' snakeTrace.batName ' ' snakeTrace.dateSesh ' ' snakeTrace.sessionType]);
    
end

%plot the cells by cluster sorted by cluster1 order
snakePlot_clustBy1_normAll = figure();
for p = 1:snakeTrace.nClusters
    p1 = subplot(4,snakeTrace.nClusters,p);
    plot(1:length(snakeTrace.smoothSpeedFlight{p}),snakeTrace.smoothSpeedFlight{p},'k');
    hold on
    for cell_ii = 1:size(snakeTrace.smoothSpeedRawFlight{p},1)
        plot(1:length(snakeTrace.smoothSpeedRawFlight{p}(cell_ii,:)),snakeTrace.smoothSpeedRawFlight{p}(cell_ii,:));
    end
    title(['Cluster ' num2str(p)]);
    if p == 1
        ylabel('velocity (cm/s)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
    else
        set(gca,'xticklabel',{[]},'yticklabel',{[]});
    end
    p2 = subplot(4,snakeTrace.nClusters,[snakeTrace.nClusters+p,2*snakeTrace.nClusters+p,3*snakeTrace.nClusters+p]);
    imagesc(snakeTrace.normTraceFlightAll{p}(snakeTrace.I1normFlightAll{p},:));
    colormap(hot);
    %make labels for first left plot only
    if p == 1
        ylabel('cell number');
    else
        set(gca,'yticklabel',{[]});
    end
    xt = get(gca, 'XTick');
    set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
    xlabel('time (s)');
    %hold off
    sgtitle(['Spatial selectivity (norm) sort by cluster 1: ' snakeTrace.batName ' ' snakeTrace.dateSesh ' ' snakeTrace.sessionType]);
    
end

%% save plotVariables to snakeTrace
snakeTrace.snakePlot_clust = snakePlot_clust;
snakeTrace.snakePlot_clustOddEven = snakePlot_clustOddEven; 
snakeTrace.snakePlot_clustBy1 = snakePlot_clustBy1;
snakeTrace.snakePlot_prefEachPrePostFlight = snakePlot_prefEachPrePostFlight;
snakeTrace.snakePlot_clustPrePostFlight = snakePlot_clustPrePostFlight;
snakeTrace.snakePlot_clust = snakePlot_clust;
snakeTrace.snakePlot_prefAll = snakePlot_prefAll;
snakeTrace.snakePlot_prefEach = snakePlot_prefEach;
snakeTrace.snakePlot_clustBy1PrePostFlight = snakePlot_clustBy1PrePostFlight;
%% save figures and matlab variable
if saveFlag == 1
    saveas(snakePlot_prefEachPrePostFlight, [pwd '\' analysis_Folder '\snakePlots\' label '_snakePlot_prefEachPrePostFlight.svg']);
    saveas(snakePlot_prefEachPrePostFlight, [pwd '\' analysis_Folder '\snakePlots\' label '_snakePlot_prefEachPrePostFlight.tif']);
    savefig(snakePlot_prefEachPrePostFlight, [pwd '/' analysis_Folder '/snakePlots/' label '_snakePlot_prefEachPrePostFlight.fig']);
    saveas(snakePlot_clust, [pwd '\' analysis_Folder '\snakePlots\' label '_snakePlots_clustAll.svg']);
    saveas(snakePlot_clustOddEven, [pwd '\' analysis_Folder '\snakePlots\' label '_snakePlots_clustOddEven.svg']);
    saveas(snakePlot_clustBy1, [pwd '\' analysis_Folder '\snakePlots\' label '_snakePlots_clustBy1.svg']);
    saveas(snakePlot_prefEach, [pwd '\' analysis_Folder '\snakePlots\' label '_snakePlots_prefEach.svg']);
    saveas(snakePlot_prefAll, [pwd '\' analysis_Folder '\snakePlots\' label '_snakePlots_prefAll.svg']);
    saveas(snakePlot_prefEachPrePostFlight, [pwd '/' analysis_Folder '/snakePlots/' fileName '_snakePlot_prefEachPrePostFlight.svg']);
    saveas(snakePlot_clust, [pwd '\' analysis_Folder '\snakePlots\' label '_snakePlots_clustAll.tif']);
    saveas(snakePlot_clustOddEven, [pwd '\' analysis_Folder '\snakePlots\' label '_snakePlots_clustOddEven.tif']);
    saveas(snakePlot_clustBy1, [pwd '\' analysis_Folder '\snakePlots\' label '_snakePlots_clustBy1.tif']);
    saveas(snakePlot_prefEach, [pwd '\' analysis_Folder '\snakePlots\' label '_snakePlots_prefEach.tif']);
    saveas(snakePlot_prefAll, [pwd '\' analysis_Folder '\snakePlots\' label '_snakePlots_prefAll.tif']);
    saveas(snakePlot_prefEachPrePostFlight, [pwd '\' analysis_Folder '\snakePlots\' fileName '_snakePlot_prefEachPrePostFlight.tif']);
    saveas(snakePlot_clustPrePostFlight, [pwd '/' analysis_Folder '/snakePlots/' fileName '_snakePlot_clustPrePostFlight.tif']);
    saveas(snakePlot_clustPrePostFlight, [pwd  '/' analysis_Folder '/snakePlots/' fileName '_snakePlot_clustPrePostFlight.svg']);
    saveas(snakePlot_prefAll, [pwd '/' analysis_Folder '/snakePlots/' fileName '_snakePlot_prefAll.tif']);
    saveas(snakePlot_prefAll, [pwd '/' analysis_Folder '/snakePlots/' fileName '_snakePlot_prefAll.svg']);
    saveas(snakePlot_clustBy1PrePostFlight, [pwd '/' analysis_Folder '/snakePlots/' fileName '_snakePlot_clustBy1PrePostFlight.tif']);
    saveas(snakePlot_clustBy1PrePostFlight, [pwd '/' analysis_Folder '/snakePlots/' fileName '_snakePlot_clustBy1PrePostFlight.svg']);

    %save([pwd '/' analysis_Folder '/' label '_snakePlotData.mat'],'snakeTrace');
    savefig(snakePlot_prefEachPrePostFlight, [pwd '\' analysis_Folder '\snakePlots\' label '_snakePlot_prefEachPrePostFlight.fig']);
    savefig(snakePlot_clustOddEven, [pwd '\' analysis_Folder '\snakePlots\' label '_snakePlots_clustOddEven.fig']);
    savefig(snakePlot_clustBy1, [pwd '\' analysis_Folder '\snakePlots\' label '_snakePlots_clustBy1.fig']);
    savefig(snakePlot_prefEach, [pwd '\' analysis_Folder '\snakePlots\' label '_snakePlots_prefEach.fig']);
    savefig(snakePlot_prefAll, [pwd '\' analysis_Folder '\snakePlots\' label '_snakePlots_prefAll.fig']);
    savefig(snakePlot_clustPrePostFlight, [pwd '/' analysis_Folder '/snakePlots/' fileName '_snakePlot_clustPrePostFlight.fig']);
    savefig(snakePlot_clustBy1PrePostFlight, [pwd '/' analysis_Folder '/snakePlots/' label '_snakePlot_clustBy1PrePostFlight.fig']);

end