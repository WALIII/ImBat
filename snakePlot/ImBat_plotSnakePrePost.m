function [snakeTracePrePost] = ImBat_plotSnakePrePost(snakeTracePrePost,varargin)
% User inputs overrides
saveFlag = 0;

nparams=length(varargin);
for i=1:2:nparams
    switch lower(varargin{i})
        case 'saveflag'
            saveFlag = varargin{i+1};
    end
end


%plot all the flight trajectories concatenated sorted by their preference
snakePlot_prefAll = figure();
imagesc(snakeTracePrePost.normMeanTraceSortPreFlightPost);
colormap(hot);
ylabel('cell number');
set(gca,'yticklabel',{[]});
xt = get(gca, 'XTick');
set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
xlabel('time (s)');
sgtitle(['Selectivity sort by group preference (pre/post/flight): ' snakeTracePrePost.batName ' ' snakeTracePrePost.dateSesh ' ' snakeTracePrePost.sessionType]);


%% plot the cells according to their peak within each grouping with velocity on top (pre/post/flight)
snakePlot_noClust_withinPrePostFlight = figure();

    p1 = subplot(4,3,1);
    plot(1:length(snakeTracePrePost.smoothSpeedPre),snakeTracePrePost.smoothSpeedPre,'k');
    hold on
    for cell_ii = 1:size(snakeTracePrePost.smoothSpeedRawFlight,1)
        plot(1:length(snakeTracePrePost.smoothSpeedRawPre(cell_ii,:)),snakeTracePrePost.smoothSpeedRawPre(cell_ii,:));
    end
    title('Pre-flight');
    ylabel('velocity (cm/s)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
    
        set(gca,'xticklabel',{[]},'yticklabel',{[]});
    
    p2 = subplot(4,3,[4 7 10]);
    imagesc(snakeTracePrePost.normTracePre(snakeTracePrePost.IPre,:));
    colormap(hot);
    %make labels for first left plot only
    ylabel('ROI number');
    %set(gca,'yticklabel',{[]});
    xt = get(gca, 'XTick');
    set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
    xlabel('time (s)');
    
    p3 = subplot(4,3,2);
    plot(1:length(snakeTracePrePost.smoothSpeedFlight),snakeTracePrePost.smoothSpeedFlight,'k');
    hold on
    for cell_ii = 1:size(snakeTracePrePost.smoothSpeedRawFlight,1)
        plot(1:length(snakeTracePrePost.smoothSpeedRawFlight(cell_ii,:)),snakeTracePrePost.smoothSpeedRawFlight(cell_ii,:));
    end
    title('During flight');
    set(gca,'xticklabel',{[]},'yticklabel',{[]});

    p4 = subplot(4,3,[5 8 11]);
    imagesc(snakeTracePrePost.normTraceFlight(snakeTracePrePost.IFlight,:));
    colormap(hot);
    %make labels for first left plot only
    xt = get(gca, 'XTick');
    set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
    xlabel('time (s)');
    set(gca,'yticklabel',{[]});

    
    p5 = subplot(4,3,3);
    plot(1:length(snakeTracePrePost.smoothSpeedPost),snakeTracePrePost.smoothSpeedPost,'k');
    hold on
    for cell_ii = 1:size(snakeTracePrePost.smoothSpeedRawPost,1)
        plot(1:length(snakeTracePrePost.smoothSpeedRawPost(cell_ii,:)),snakeTracePrePost.smoothSpeedRawPost(cell_ii,:));
    end
    title('Post-flight');
    set(gca,'xticklabel',{[]},'yticklabel',{[]});

    p6 = subplot(4,3,[6 9 12]);
    imagesc(snakeTracePrePost.normTracePost(snakeTracePrePost.IPost,:));
    colormap(hot);
    %make labels for first left plot only
    set(gca,'yticklabel',{[]});
    xt = get(gca, 'XTick');
    set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
    xlabel('time (s)');
    
    sgtitle(['Selectivity sort within groupings (pre/post/flight): ' snakeTracePrePost.batName ' ' snakeTracePrePost.dateSesh ' ' snakeTracePrePost.sessionType]);
  
%% plot the cells according to their peak across the groupings with velocity on top (pre/post/flight)
snakePlot_noClust_acrossPrePostFlight = figure();

    p1 = subplot(4,3,1);
    plot(1:length(snakeTracePrePost.smoothSpeedPre),snakeTracePrePost.smoothSpeedPre,'k');
    hold on
    for cell_ii = 1:size(snakeTracePrePost.smoothSpeedRawFlight,1)
        plot(1:length(snakeTracePrePost.smoothSpeedRawPre(cell_ii,:)),snakeTracePrePost.smoothSpeedRawPre(cell_ii,:));
    end
    title('Pre-flight');
    ylabel('velocity (cm/s)');
    %yt = get(gca,'YTick');
    %set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
    set(gca,'xticklabel',{[]},'yticklabel',{[]});
    
    p2 = subplot(4,3,[4 7 10]);
    imagesc(snakeTracePrePost.normMeanTraceEachPre);
    colormap(hot);
    %make labels for first left plot only
    ylabel('ROI number');
    %set(gca,'yticklabel',{[]});
    xt = get(gca, 'XTick');
    set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
    xlabel('time (s)');
    
    p3 = subplot(4,3,2);
    plot(1:length(snakeTracePrePost.smoothSpeedFlight),snakeTracePrePost.smoothSpeedFlight,'k');
    hold on
    for cell_ii = 1:size(snakeTracePrePost.smoothSpeedRawFlight,1)
        plot(1:length(snakeTracePrePost.smoothSpeedRawFlight(cell_ii,:)),snakeTracePrePost.smoothSpeedRawFlight(cell_ii,:));
    end
    title('During flight');
    set(gca,'xticklabel',{[]},'yticklabel',{[]});

    p4 = subplot(4,3,[5 8 11]);
    imagesc(snakeTracePrePost.normMeanTraceEachFlight);
    colormap(hot);
    %make labels for first left plot only
    set(gca,'yticklabel',{[]});
    xt = get(gca, 'XTick');
    set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
    xlabel('time (s)');

    
    p5 = subplot(4,3,3);
    plot(1:length(snakeTracePrePost.smoothSpeedPost),snakeTracePrePost.smoothSpeedPost,'k');
    hold on
    for cell_ii = 1:size(snakeTracePrePost.smoothSpeedRawPost,1)
        plot(1:length(snakeTracePrePost.smoothSpeedRawPost(cell_ii,:)),snakeTracePrePost.smoothSpeedRawPost(cell_ii,:));
    end
    title('Post-flight');
    set(gca,'xticklabel',{[]},'yticklabel',{[]});

    p6 = subplot(4,3,[6 9 12]);
    imagesc(snakeTracePrePost.normMeanTraceEachPost);
    colormap(hot);
    %make labels for first left plot only
    set(gca,'yticklabel',{[]});
    xt = get(gca, 'XTick');
    set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
    xlabel('time (s)');
    
    sgtitle(['Selectivity sort across groups (pre/post/flight): ' snakeTracePrePost.batName ' ' snakeTracePrePost.dateSesh ' ' snakeTracePrePost.sessionType]);
      