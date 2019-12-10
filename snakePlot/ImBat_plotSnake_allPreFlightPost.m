function [snakeTrace] = ImBat_plotSnake_allPreFlightPost(snakeTrace)

% plot the cells according to their peak for each cluster with velocity on top (pre/post/flight)
snakePlot_allPrePostFlight = figure();
for p = 1:snakeTrace.nClusters
    p1 = subplot(snakeTrace.nClusters*3,3,((p-1)*9)+1);
    plot(1:length(snakeTrace.smoothSpeedPre{p}),snakeTrace.smoothSpeedPre{p},'k');
    hold on
    for cell_ii = 1:size(snakeTrace.smoothSpeedRawPre{p},1)
        plot(1:length(snakeTrace.smoothSpeedRawPre{p}(cell_ii,:)),snakeTrace.smoothSpeedRawPre{p}(cell_ii,:));
    end
    if p ==1
    title('Pre-Flight');
    end
    ylabel('velocity (cm/s)');
    yt = get(gca,'YTick');
    set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
    set(gca,'xticklabel',{[]});
    xlim([1 length(snakeTrace.smoothSpeedRawPre{p}(cell_ii,:))]);
    
    p2 = subplot(snakeTrace.nClusters*3,3,((p-1)*9)+2)
    plot(1:length(snakeTrace.smoothSpeedFlight{p}),snakeTrace.smoothSpeedFlight{p},'k');
    hold on
    for cell_ii = 1:size(snakeTrace.smoothSpeedRawFlight{p},1)
        plot(1:length(snakeTrace.smoothSpeedRawFlight{p}(cell_ii,:)),snakeTrace.smoothSpeedRawFlight{p}(cell_ii,:));
    end
    if p ==1
    title('Flight');
    end
    set(gca,'xticklabel',{[]});
    xlim([1 length(snakeTrace.smoothSpeedRawPre{p}(cell_ii,:))]);
    ylabel('velocity (cm/s)');
    yt = get(gca,'YTick');
    set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
    
    p3 = subplot(snakeTrace.nClusters*3,3,((p-1)*9)+3)
    plot(1:length(snakeTrace.smoothSpeedPost{p}),snakeTrace.smoothSpeedPost{p},'k');
    hold on
    for cell_ii = 1:size(snakeTrace.smoothSpeedRawPost{p},1)
        plot(1:length(snakeTrace.smoothSpeedRawPost{p}(cell_ii,:)),snakeTrace.smoothSpeedRawPost{p}(cell_ii,:));
    end
    if p ==1
    title('Post-Flight');
    end
    set(gca,'xticklabel',{[]});
    xlim([1 length(snakeTrace.smoothSpeedRawPre{p}(cell_ii,:))]);
    ylabel('velocity (cm/s)');
    yt = get(gca,'YTick');
    set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});

    
    p4 = subplot(snakeTrace.nClusters*3,3,[((p-1)*9+4) ((p-1)*9+7)])
    imagesc(snakeTrace.normMeanTraceEachPreFP{p});
    colormap(hot);
    %make labels for first left plot only
    ylabel('ROI #');
    xt = get(gca, 'XTick');
    set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
    xlabel('time (s)');
    
    p5 = subplot(snakeTrace.nClusters*3,3,[((p-1)*9+5) ((p-1)*9+8)])
    imagesc(snakeTrace.normMeanTraceEachPFlightP{p});
    colormap(hot);
    %make labels for first left plot only
    set(gca,'yticklabel',{[]});
    xt = get(gca, 'XTick');
    set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
    xlabel('time (s)');
    
    p6 = subplot(snakeTrace.nClusters*3,3,[((p-1)*9+6) ((p-1)*9+9)])
    imagesc(snakeTrace.normMeanTraceEachPFPost{p});
    colormap(hot);
    %make labels for first left plot only
    set(gca,'yticklabel',{[]});
    xt = get(gca, 'XTick');
    set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
    xlabel('time (s)');
    
    sgtitle('ROI Activity for all cells resorted for each cluster')

end

snakeTrace.snakePlot_allPrePostFlight = snakePlot_allPrePostFlight;