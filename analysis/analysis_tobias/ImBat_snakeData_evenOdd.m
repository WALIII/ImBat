function [snakeTrace] = ImBat_snakeData_evenOdd(snakeTrace)

for c = 1:snakeTrace.nClusters

    snakeTrace.normMeanTraceSortPreFlightPost_odd{c} = snakeTrace.normMeanTraceSortPreFlightPost{c}(1:2:end,:);
    snakeTrace.normMeanTraceSortPreFlightPost_even{c} = snakeTrace.normMeanTraceSortPreFlightPost{c}(2:2:end,:);
    
    snakeTrace.smoothSpeedPreFlightPost{c} = vertcat(snakeTrace.smoothSpeedPre{c},snakeTrace.smoothSpeedFlight{c},snakeTrace.smoothSpeedPost{c});
    snakeTrace.smoothSpeedPreFlightPost_odd{c} = snakeTrace.smoothSpeedPreFlightPost{c}(1:2:end,:);
    snakeTrace.smoothSpeedPreFlightPost_even{c} = snakeTrace.smoothSpeedPreFlightPost{c}(2:2:end,:);
    
    snakeTrace.smoothSpeedRawPreFlightPost{c} = [snakeTrace.smoothSpeedRawPre{c} snakeTrace.smoothSpeedRawFlight{c} snakeTrace.smoothSpeedRawPost{c}];
    snakeTrace.smoothSpeedRawPreFlightPost_odd{c} = snakeTrace.smoothSpeedRawPreFlightPost{c}(1:2:end,:);
    snakeTrace.smoothSpeedRawPreFlightPost_even{c} = snakeTrace.smoothSpeedRawPreFlightPost{c}(2:2:end,:);
    
    
for cell_i = 1:length(snakeTrace.normMeanTraceSortPreFlightPost_odd{c}(:,1))
    [~,snakeTrace.maxNormTracePreFlightPostOdd{c}(cell_i,1)] = max(snakeTrace.normMeanTraceSortPreFlightPost_odd{c}(cell_i,:)); 

end

[snakeTrace.BoddPreFlightPost{c},snakeTrace.IoddPreFlightPost{c}] = sort(snakeTrace.maxNormTracePreFlightPostOdd{c});

end

plotSnake_oddPreFlightPost = figure();
 
for p = 1:snakeTrace.nClusters
    p1 = subplot(snakeTrace.nClusters*3,2,((p-1)*6)+1);
    plot(1:length(snakeTrace.smoothSpeedPreFlightPost_odd{p}),snakeTrace.smoothSpeedPreFlightPost_odd{p},'k');
    hold on
    for cell_ii = 1:size(snakeTrace.smoothSpeedRawPreFlightPost_odd{p},1)
        plot(1:length(snakeTrace.smoothSpeedRawPreFlightPost_odd{p}(cell_ii,:)),snakeTrace.smoothSpeedRawPreFlightPost_odd{p}(cell_ii,:));
    end
    if p ==1
    title('Odd flights');
    end
    ylabel('velocity (cm/s)');
    yt = get(gca,'YTick');
    set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
    set(gca,'xticklabel',{[]});
    xlim([1 length(snakeTrace.smoothSpeedRawPreFlightPost_odd{p}(cell_ii,:))]);
    
    p2 = subplot(snakeTrace.nClusters*3,2,((p-1)*6)+2);
    plot(1:length(snakeTrace.smoothSpeedPreFlightPost_even{p}),snakeTrace.smoothSpeedPreFlightPost_even{p},'k');
    hold on
    for cell_ii = 1:size(snakeTrace.smoothSpeedRawPreFlightPost_even{p},1)
        plot(1:length(snakeTrace.smoothSpeedRawPreFlightPost_even{p}(cell_ii,:)),snakeTrace.smoothSpeedRawPreFlightPost_even{p}(cell_ii,:));
    end
    if p ==1
    title('Even flights');
    end
    ylabel('velocity (cm/s)');
    yt = get(gca,'YTick');
    set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
    set(gca,'xticklabel',{[]});
    xlim([1 length(snakeTrace.smoothSpeedRawPreFlightPost_even{p}(cell_ii,:))]);
    
    p3 = subplot(snakeTrace.nClusters*3,2,[((p-1)*6+3) ((p-1)*6+5)]);
    imagesc(snakeTrace.normMeanTraceSortPreFlightPost_odd{p}(snakeTrace.IoddPreFlightPost{c},:));
    colormap(hot);
    %make labels for first left plot only
    ylabel('ROI #');
    xt = get(gca, 'XTick');
    set(gca,'XTick',xt,'XTickLabel',round(xt/30,1)); 
    
    p3 = subplot(snakeTrace.nClusters*3,2,[((p-1)*6+4) ((p-1)*6+6)]);
    imagesc(snakeTrace.normMeanTraceSortPreFlightPost_even{p}(snakeTrace.IoddPreFlightPost{c},:));
    colormap(hot);
    %make labels for first left plot only
    ylabel('ROI #');
    xt = get(gca, 'XTick');
    set(gca,'XTick',xt,'XTickLabel',round(xt/30,1)); 
end
sgtitle(['Sort Odd vs Even Pre/post/flight: ' snakeTrace.batName ' ' snakeTrace.dateSesh ' ' snakeTrace.sessionType]);



