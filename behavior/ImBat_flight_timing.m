function ImBat_flight_timing(flightPaths,snakeTrace,flightTransitions,varargin)

%% Plot the time differences between flights A to B & BtoA
histAtoBFlag = 1;
boxplotAtoBFlag = 1;
plotFlightTimesFlag = 1;
saveFlag = 1;

nparams=length(varargin);
for i=1:2:nparams
    switch lower(varargin{i})
        case 'saveflag'
            saveFlag = varargin{i+1};
    end
end


%calculate timing difference from start/end of flight A to start of flight B
for flight_i = 1:length(flightTransitions.AtoB(1,:))
    %timing from flight B to A
    flightStartA(flight_i) = flightPaths.flight_starts_idx(flightTransitions.AtoB(flight_i));
    flightStartB(flight_i) = flightPaths.flight_starts_idx(flightTransitions.AtoB(flight_i)+1);
    flightEndA(flight_i) = flightPaths.flight_ends_idx(flightTransitions.AtoB(flight_i));
    flightEndB(flight_i) = flightPaths.flight_ends_idx(flightTransitions.AtoB(flight_i)+1);
    
    
    diffStartAtoStartB(flight_i) = flightStartB(flight_i) - flightStartA(flight_i);
    diffStartAtoStartB(flight_i) = diffStartAtoStartB(flight_i)/120; %convert to seconds
    diffEndAtoStartB(flight_i) = flightStartB(flight_i) - flightEndA(flight_i);
    diffEndAtoStartB(flight_i) = diffEndAtoStartB(flight_i)/120; %convert to seconds
end

%calculate timing difference from start/end of flight B to start of flight A
for flight_i = 1:length(flightTransitions.BtoA(1,:))
    %timing from flight B to A
    flightSB(flight_i) = flightPaths.flight_starts_idx(flightTransitions.BtoA(flight_i));
    flightSA(flight_i) = flightPaths.flight_starts_idx(flightTransitions.BtoA(flight_i)+1);
    flightEB(flight_i) = flightPaths.flight_ends_idx(flightTransitions.BtoA(flight_i));
    flightEA(flight_i) = flightPaths.flight_ends_idx(flightTransitions.BtoA(flight_i)+1);
    
    diffSBtoSA(flight_i) = flightSA(flight_i) - flightSB(flight_i);
    diffSBtoSA(flight_i) = diffSBtoSA(flight_i)/120; %convert to seconds
    diffEBtoSA(flight_i) = flightSA(flight_i) - flightEB(flight_i);
    diffEBtoSA(flight_i) = diffEBtoSA(flight_i)/120; %convert to seconds
    
    %     flightConcat = [flightConcat flightPaths.trajectories_continuous(1,flightPaths.flight_starts_idx(flightTransitions.AtoB(flight_i))-snakeTrace.preFlightPadSpeed...
    %         :flightPaths.flight_ends_idx(flightTransitions.AtoB(flight_i))),flightPaths.trajectories_continuous(2,flightPaths.flight_starts_idx(flightTransitions.AtoB(flight_i))-snakeTrace.preFlightPadSpeed:flightPaths.flight_ends_idx(flightTransitions.AtoB(flight_i))+snakeTrace.postFlightPadSpeed),'LineWidth',1,'Color','b')%,jj(traj*100,:))
    
end

%plot the 6 histograms of timing from A to B and from B to A
if histAtoBFlag == 1
    nBins = 9; %number of bins to use for the histogram
    hist_flightTiming_AtoB = figure('units','normalized','outerposition',[0.5 0 0.5 1]);
    subplot(3,2,1);
    histogram(diffStartAtoStartB,nBins);
    title('Timing Start Flight A to Start Flight B');
    ylabel('Num flights');
    xlabel('Time between flights (sec)');
    
    subplot(3,2,2);
    histogram(diffEndAtoStartB,nBins);
    title('Timing End Flight A to Start Flight B');
    ylabel('Num flights');
    xlabel('Time between flights (sec)');
    
    subplot(3,2,3);
    histogram(diffSBtoSA,nBins,'FaceColor','r');
    title('Timing Start Flight B to Start Flight A');
    ylabel('Num flights');
    xlabel('Time between flights (sec)');
    
    subplot(3,2,4);
    histogram(diffEBtoSA,nBins,'FaceColor','r');
    title('Timing End Flight B to Start Flight A');
    ylabel('Num flights');
    xlabel('Time between flights (sec)');
    
    subplot(3,2,5);
    histogram(diffSBtoSA(find(diffSBtoSA<40)),6,'FaceColor','r');
    title('Timing Start Flight B to Start Flight A');
    ylabel('Num flights');
    xlabel('Time between flights (sec)');
    
    subplot(3,2,6);
    histogram(diffEBtoSA(find(diffEBtoSA<40)),6,'FaceColor','r');
    title('Timing End Flight B to Start Flight A');
    ylabel('Num flights');
    xlabel('Time between flights (sec)');
    
    sgtitle(['Flight timing/jitter for A/B flights: ' snakeTrace.batName ' ' snakeTrace.dateSesh ' ' snakeTrace.sessionType]);
    if saveFlag ==1
        saveas(hist_flightTiming_AtoB, [pwd '/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_hist_flightTiming_AtoB.svg']);
        saveas(hist_flightTiming_AtoB, [pwd '/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_hist_flightTiming_AtoB.tif']);
        savefig(hist_flightTiming_AtoB, [pwd '/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_hist_flightTiming_AtoB.fig']);
    end
    
end

%plot the 6 boxplots of timing from A to B and from B to A
if boxplotAtoBFlag == 1    
    boxplot_flightTiming_AtoB = figure('units','normalized','outerposition',[0 0 0.5 1]);
    subplot(3,2,1);
    boxplot(diffStartAtoStartB);
    title('Timing Start Flight A to Start Flight B');
    ylabel('Num flights');
    xlabel('Time between flights (sec)');
    
    subplot(3,2,2);
    boxplot(diffEndAtoStartB);
    title('Timing End Flight A to Start Flight B');
    ylabel('Num flights');
    xlabel('Time between flights (sec)');
    
    subplot(3,2,3);
    boxplot(diffSBtoSA);
    title('Timing Start Flight B to Start Flight A');
    ylabel('Num flights');
    xlabel('Time between flights (sec)');
    
    subplot(3,2,4);
    boxplot(diffEBtoSA);
    title('Timing End Flight B to Start Flight A');
    ylabel('Num flights');
    xlabel('Time between flights (sec)');
    
    subplot(3,2,5);
    boxplot(diffSBtoSA(find(diffSBtoSA<40)));
    title('Timing Start Flight B to Start Flight A');
    ylabel('Num flights');
    xlabel('Time between flights (sec)');
    
    subplot(3,2,6);
    boxplot(diffEBtoSA(find(diffEBtoSA<40)));
    title('Timing End Flight B to Start Flight A');
    ylabel('Num flights');
    xlabel('Time between flights (sec)');
    
    sgtitle(['Flight timing/jitter for A/B flights: ' snakeTrace.batName ' ' snakeTrace.dateSesh ' ' snakeTrace.sessionType]);
    if saveFlag ==1
        saveas(boxplot_flightTiming_AtoB, [pwd '/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_boxplot_flightTiming_AtoB.svg']);
        saveas(boxplot_flightTiming_AtoB, [pwd '/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_boxplot_flightTiming_AtoB.tif']);
        savefig(boxplot_flightTiming_AtoB, [pwd '/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_boxplot_flightTiming_AtoB.fig']);
    end
    
end

%% Plot each of the flights as a star in a timeline with a different color for each traj cluster

if plotFlightTimesFlag ==1
    jj = jet(length(flightPaths.clusterIndex)*1000); %make color vector
    plotFlightTimes = figure('units','normalized','outerposition',[0 0 1 0.5]);
    for clust_i = 1:length(flightPaths.clusterIndex)
        %make height vector for the plot
        ht = ones(1,length(flightPaths.clusterIndex{clust_i}));
        %plot each flight in the color of the trajectory
        for flight_i = 1:length(flightPaths.clusterIndex{clust_i})
            legendData{clust_i} = plot(flightPaths.flight_starts_idx(flightPaths.clusterIndex{clust_i}(flight_i)),ht,'*','Color',jj(clust_i*1000,:));
            hold on;
        end
        %make dummy points for the legend so the colors match the trajectories
        l(clust_i) = plot([NaN,NaN],'*', 'color', jj(clust_i*1000,:));
        
    end
    title(['Flight times of each trajectory: ' snakeTrace.batName ' ' snakeTrace.dateSesh ' ' snakeTrace.sessionType]);
    xt = get(gca, 'XTick');
    set(gca,'XTick',xt,'XTickLabel',round(xt/120,1));
    xlabel('time (s)');
    set(gca,'yticklabel',{[]});
    legend([l(1) l(2) l(3) l(4) l(5) l(6) l(7) l(8) l(9) l(10) l(11) l(12)],{'traj 1','traj 2','traj 3','traj 4','traj 5','traj 6','traj 7','traj 8','traj 9','traj 10','traj 11','traj 12'});
    
    if saveFlag ==1
        saveas(plotFlightTimes, [pwd '/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_plotFlightTimes.svg']);
        saveas(plotFlightTimes, [pwd '/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_plotFlightTimes.tif']);
        savefig(plotFlightTimes, [pwd '/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_plotFlightTimes.fig']);
    end
end