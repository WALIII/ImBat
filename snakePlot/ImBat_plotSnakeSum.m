function [plotSnakeSum] = ImBat_plotSnakeSum(snakeTrace,flightPaths,goodCellIdx,varargin)
% User inputs overrides
saveFlag = 1;
plotSumFlag = 1;
plotCumSumFlag = 1;
plotSumClusterFlag = 1;

nparams=length(varargin);
for i=1:2:nparams
    switch lower(varargin{i})
        case 'saveflag'
            saveFlag = varargin{i+1};
    end
end

%this is to merge specific clusters if 2 or more of them are the same but
%were separated by the flight k-means
flightPaths.clusterIndex{2}=cat(2,flightPaths.clusterIndex{2},flightPaths.clusterIndex{3});
flightPaths.clusterIndex{3} =[];
flightPaths.clusterIndex= flightPaths.clusterIndex(~cellfun('isempty',flightPaths.clusterIndex));

%calculate the sums and cumsums for each flight trajectory
%initialize variables for sum and cumsum
sumSnakeFlight = cell(size(snakeTrace.meanTraceFlight,2),1);
sumSnakePost = cell(size(snakeTrace.meanTraceFlight,2),1);
cumSumSnakePre = cell(size(snakeTrace.meanTraceFlight,2),1);
cumSumSnakeFlight = cell(size(snakeTrace.meanTraceFlight,2),1);
cumSumSnakePost = cell(size(snakeTrace.meanTraceFlight,2),1);
%calculate the sums and cumsums for each trajectory pre/flight/post
for traj_i = 1:size(snakeTrace.meanTraceFlight,2)
    sumSnakePre{traj_i} = sum(snakeTrace.meanTracePre{traj_i}(goodCellIdx.goodCellIndex,:),1);
    sumSnakeFlight{traj_i} = sum(snakeTrace.meanTraceFlight{traj_i}(goodCellIdx.goodCellIndex,:),1);
    sumSnakePost{traj_i} = sum(snakeTrace.meanTracePost{traj_i}(goodCellIdx.goodCellIndex,:),1);
    
    cumSumSnakePre{traj_i} = cumsum(snakeTrace.meanTracePre{traj_i}(goodCellIdx.goodCellIndex,:),1);
    cumSumSnakeFlight{traj_i} = cumsum(snakeTrace.meanTraceFlight{traj_i}(goodCellIdx.goodCellIndex,:),1);
    cumSumSnakePost{traj_i} = cumsum(snakeTrace.meanTracePost{traj_i}(goodCellIdx.goodCellIndex,:),1);
end

%% plot all the sums of calcium activity for pre/flight/post separated for
if plotSumFlag ==1
    %each flight cluster
    plotSnakeSum = figure();
    for traj_i = 1:size(snakeTrace.meanTraceFlight,2)
        %plot sums for preflight times
        p1 = subplot(size(snakeTrace.meanTraceFlight,2),3,(2*traj_i)+(traj_i-2));
        plot(sumSnakePre{traj_i});
        hold on
        set(gca,'xticklabel',{[]});
        ylabel({['Traj ' num2str(traj_i)];'sum dff'});
        if traj_i == 1
            title('Pre-flight DF Sum')
        elseif traj_i == size(snakeTrace.meanTraceFlight,2)
            xt = get(gca, 'XTick');
            set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
            xlabel('time (s)');
        end
        hold off
        
        %plot sums for flight time
        p2 = subplot(size(snakeTrace.meanTraceFlight,2),3,(2*traj_i)+(traj_i-1));
        plot(sumSnakeFlight{traj_i});
        hold on
        set(gca,'xticklabel',{[]});
        if traj_i == 1
            title('Flight DF Sum')
        elseif traj_i == size(snakeTrace.meanTraceFlight,2)
            xt = get(gca, 'XTick');
            set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
            xlabel('time (s)');
        end
        hold off
        
        %plot sums for post flight
        p3 = subplot(size(snakeTrace.meanTraceFlight,2),3,(2*traj_i)+traj_i)
        plot(sumSnakePost{traj_i})
        hold on
        set(gca,'xticklabel',{[]});
        if traj_i == 1
            title('Post-flight DF Sum')
        elseif traj_i == size(snakeTrace.meanTraceFlight,2)
            xt = get(gca, 'XTick');
            set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
            xlabel('time (s)');
        end
        hold off
    end
    sgtitle(['Sum of calcium activity: ' snakeTrace.batName ' ' snakeTrace.dateSesh ' ' snakeTrace.sessionType]);
    
    if saveFlag ==1
        saveas(plotSnakeSum, [pwd '/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_plotSumPreFlightPost.svg']);
        saveas(plotSnakeSum, [pwd '/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_plotSumPreFlightPost.tif']);
        savefig(plotSnakeSum, [pwd '/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_plotSumPreFlightPost.fig']);
    end
end
%% plot all the cumsums of calcium activity for pre/flight/post separated for each flight cluster
if plotCumSumFlag == 1
    plotSnakeCumSum = figure();
    for traj_i = 1:size(snakeTrace.meanTraceFlight,2)
        %plot cumsum of preflight activity
        p1 = subplot(size(snakeTrace.meanTraceFlight,2),3,(2*traj_i)+(traj_i-2))
        plot(cumSumSnakePre{traj_i})
        hold on
        set(gca,'xticklabel',{[]});
        ylabel({['Traj ' num2str(traj_i)];'sum dff'});
        if traj_i == 1
            title('Pre-flight DF Sum')
        elseif traj_i == size(snakeTrace.meanTraceFlight,2)
            xt = get(gca, 'XTick');
            set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
            xlabel('time (s)');
        end
        hold off
        
        %plot cumsum of flight activity
        p2 = subplot(size(snakeTrace.meanTraceFlight,2),3,(2*traj_i)+(traj_i-1))
        plot(cumSumSnakeFlight{traj_i})
        hold on
        set(gca,'xticklabel',{[]});
        if traj_i == 1
            title('Flight DF Sum')
        elseif traj_i == size(snakeTrace.meanTraceFlight,2)
            xt = get(gca, 'XTick');
            set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
            xlabel('time (s)');
        end
        hold off
        
        %plot cumsum of postflight activity
        p3 = subplot(size(snakeTrace.meanTraceFlight,2),3,(2*traj_i)+traj_i)
        plot(cumSumSnakePost{traj_i})
        hold on
        set(gca,'xticklabel',{[]});
        if traj_i == 1
            title('Post-flight DF Sum')
        elseif traj_i == size(snakeTrace.meanTraceFlight,2)
            xt = get(gca, 'XTick');
            set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
            xlabel('time (s)');
        end
        hold off
    end
    sgtitle(['Cummulative sum of calcium activity: ' snakeTrace.batName ' ' snakeTrace.dateSesh ' ' snakeTrace.sessionType]);
    if saveFlag ==1
        saveas(plotSnakeCumSum, [pwd '/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_plotCumSumPreFlightPost.svg']);
        saveas(plotSnakeCumSum, [pwd '/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_plotCumSumPreFlightPost.tif']);
        savefig(plotSnakeCumSum, [pwd '/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_plotCumSumPreFlightPost.fig']);
    end
end

%% plot the sums with the flight trajectories, speed, and snakeplot included
% plot the cells according to their peak within each grouping with velocity on top (pre/post/flight)

if plotSumClusterFlag ==1
    jj = jet(size(snakeTrace.meanTraceFlight,2));
    %for each flight cluster, plot a new figure
    for traj_i = 1:size(snakeTrace.meanTraceFlight,2)
        plotSnakeSum_noClust = figure();
        
        for flight_i = 1:length(snakeTrace.smoothSpeedRawPre{traj_i}(:,1))
            %plot trajectories preflight
            p1 = subplot(7,3,1);
            plot(flightPaths.trajectories_continuous(1,flightPaths.flight_starts_idx(flightPaths.clusterIndex{traj_i}(flight_i))-snakeTrace.preFlightPadSpeed:flightPaths.flight_starts_idx(flightPaths.clusterIndex{traj_i}(flight_i))),flightPaths.trajectories_continuous(2,flightPaths.flight_starts_idx(flightPaths.clusterIndex{traj_i}(flight_i))-snakeTrace.preFlightPadSpeed:flightPaths.flight_starts_idx(flightPaths.clusterIndex{traj_i}(flight_i))),'LineWidth',1,'Color',jj(traj_i,:));
            hold on
            %scatter(flightPaths.flight_starts_xyz(flightPaths.clusterIndex{traj_i}(flight_i),1),flightPaths.flight_starts_xyz(flightPaths.clusterIndex{traj_i}(flight_i),2),50,'r','filled')
            title('Pre-Flight: A to B');
            ylabel('y (cm)');
            yt = get(gca,'YTick');
            set(gca,'YTick',yt,'YTickLabel',yt/10);
            xlabel('x (cm)');
            xt = get(gca,'XTick');
            set(gca,'XTick',xt,'XTickLabel',xt/10);
            %plot trajectories of flight cluster (traj_i)
            p2 = subplot(7,3,2);
            plot(flightPaths.trajectories_continuous(1,flightPaths.flight_starts_idx(flightPaths.clusterIndex{traj_i}(flight_i)):flightPaths.flight_ends_idx(flightPaths.clusterIndex{traj_i}(flight_i))),flightPaths.trajectories_continuous(2,flightPaths.flight_starts_idx(flightPaths.clusterIndex{traj_i}(flight_i)):flightPaths.flight_ends_idx(flightPaths.clusterIndex{traj_i}(flight_i))),'LineWidth',1,'Color',jj(traj_i,:));
            hold on
            scatter(flightPaths.flight_starts_xyz(flightPaths.clusterIndex{traj_i}(flight_i),1),flightPaths.flight_starts_xyz(flightPaths.clusterIndex{traj_i}(flight_i),2),50,'r','filled')
            title('Pre-Flight: A to B');
            %ylabel('y (cm)');
            yt = get(gca,'YTick');
            set(gca,'YTick',yt,'YTickLabel',yt/10);
            %xlabel('x (cm)');
            xt = get(gca,'XTick');
            set(gca,'XTick',xt,'XTickLabel',xt/10);
            %plot trajectories post flight
            p3 = subplot(7,3,3);
            plot(flightPaths.trajectories_continuous(1,flightPaths.flight_ends_idx(flightPaths.clusterIndex{traj_i}(flight_i)):flightPaths.flight_ends_idx(flightPaths.clusterIndex{traj_i}(flight_i))+snakeTrace.postFlightPadSpeed),flightPaths.trajectories_continuous(2,flightPaths.flight_ends_idx(flightPaths.clusterIndex{traj_i}(flight_i)):flightPaths.flight_ends_idx(flightPaths.clusterIndex{traj_i}(flight_i))+snakeTrace.postFlightPadSpeed),'LineWidth',1,'Color',jj(traj_i,:));
            hold on
            %scatter(flightPaths.flight_starts_xyz(flightPaths.clusterIndex{traj_i}(flight_i),1),flightPaths.flight_starts_xyz(flightPaths.clusterIndex{traj_i}(flight_i),2),50,'r','filled')
            title('Pre-Flight: A to B');
            %ylabel('y (cm)');
            yt = get(gca,'YTick');
            set(gca,'YTick',yt,'YTickLabel',yt/10);
            %xlabel('x (cm)');
            xt = get(gca,'XTick');
            set(gca,'XTick',xt,'XTickLabel',xt/10);
            
            %plot speed pre flight
            p4 = subplot(7,3,4);
            plot(1:length(snakeTrace.smoothSpeedRawPre{traj_i}(flight_i,:)),snakeTrace.smoothSpeedRawPre{traj_i}(flight_i,:));
            hold on
            ylabel('cm/s');
            yt = get(gca,'YTick');
            set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
            set(gca,'xticklabel',{[]});
            ylim([0 6]);
            %plot speed during flight
            p5 = subplot(7,3,5); %plot speed pre flight A to B
            plot(1:length(snakeTrace.smoothSpeedRawFlight{traj_i}(flight_i,:)),snakeTrace.smoothSpeedRawFlight{traj_i}(flight_i,:));
            hold on
            yt = get(gca,'YTick');
            set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
            set(gca,'xticklabel',{[]});
            ylim([0 6]);
            %plot speed post flight
            p6 = subplot(7,3,6); %plot speed pre flight A to B
            plot(1:length(snakeTrace.smoothSpeedRawPost{traj_i}(flight_i,:)),snakeTrace.smoothSpeedRawPost{traj_i}(flight_i,:));
            hold on
            yt = get(gca,'YTick');
            set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
            ylim([0 6]);
        end
        
        %get rid of bad cells
        for bad_i = 1:length(goodCellIdx.badCellIndex)
        snakeTrace.IPre{traj_i}(find(snakeTrace.IPre{traj_i}==goodCellIdx.badCellIndex(bad_i)))=[];
        snakeTrace.IFlight{traj_i}(find(snakeTrace.IFlight{traj_i}==goodCellIdx.badCellIndex(bad_i)))=[];
        snakeTrace.IPost{traj_i}(find(snakeTrace.IPost{traj_i}==goodCellIdx.badCellIndex(bad_i)))=[];
        end
        
        %plot sums for preflight times
        p7 = subplot(7,3,7);
        plot(sumSnakePre{traj_i});
        hold on
        ylabel('sum dff');
        xt = get(gca, 'XTick');
        set(gca,'XTick',xt,'XTickLabel',[]);
        hold off
        %plot sums for flight time
        p8 = subplot(7,3,8);
        plot(sumSnakeFlight{traj_i});
        hold on
        set(gca,'xticklabel',{[]});
        hold off
        %plot sums for post flight
        p9 = subplot(7,3,9);
        plot(sumSnakePost{traj_i});
        hold on
        set(gca,'xticklabel',{[]});
        hold off
        
        %plot sums for preflight times
        p7 = subplot(7,3,10);
        plot(cumSumSnakePre{traj_i});
        hold on
        ylabel('cumsum dff');
        xt = get(gca, 'XTick');
        set(gca,'XTick',xt,'XTickLabel',[]);
        hold off
        %plot sums for flight time
        p8 = subplot(7,3,11);
        plot(cumSumSnakeFlight{traj_i});
        hold on
        set(gca,'xticklabel',{[]});
        hold off
        %plot sums for post flight
        p9 = subplot(7,3,12);
        plot(cumSumSnakePost{traj_i});
        hold on
        set(gca,'xticklabel',{[]});
        hold off
        
        %plot heat maps preflight
        p10 = subplot(7,3,[13 16 19]);
        imagesc(snakeTrace.normTracePre{traj_i}(snakeTrace.IPre{traj_i},:));
        colormap(hot);
        ylabel('ROI number');
        xt = get(gca, 'XTick');
        set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
        xlabel('time (s)');
        %plot head maps during flight
        p11 = subplot(7,3,[14 17 20]);
        imagesc(snakeTrace.normTraceFlight{traj_i}(snakeTrace.IFlight{traj_i},:));
        colormap(hot);
        xt = get(gca, 'XTick');
        set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
        xlabel('time (s)');
        set(gca,'yticklabel',{[]});
        %plot head maps post flight
        p12 = subplot(7,3,[15 18 21]);
        imagesc(snakeTrace.normTracePost{traj_i}(snakeTrace.IPost{traj_i},:));
        colormap(hot);
        xt = get(gca, 'XTick');
        set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
        xlabel('time (s)');
        set(gca,'yticklabel',{[]});
        
        sgtitle(['Trajectory ' num2str(traj_i) ': ' snakeTrace.batName ' ' snakeTrace.dateSesh ' ' snakeTrace.sessionType]);
        if saveFlag ==1
            saveas(plotSnakeSum_noClust, [pwd '/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_plotSnakeSum_noClust_traj' num2str(traj_i) '.svg']);
            saveas(plotSnakeSum_noClust, [pwd '/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_plotSnakeSum_noClust_traj' num2str(traj_i) '.tif']);
            savefig(plotSnakeSum_noClust, [pwd '/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_plotSnakeSum_noClust_traj' num2str(traj_i) '.fig']);
        end
        
    end
end
