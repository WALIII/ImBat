function [traceAnalysis] = ImBat_plotFlightvsCells_transitions(snakeTrace,flightTransitions,flightPaths,goodCellIdx,flightTiming,varargin)
plotEachFlag = 0;
plotPreFlightPostFlag = 0;
plotMeanSpikeFlag = 0;
plotMeanSortedFlag = 0;
histoTraceStatsFlag = 0;
plotSortedFWHMFlag = 0;
plotTimeJitterSortFlag = 1;

batName = snakeTrace.batName;
dateSesh = snakeTrace.dateSesh;
sessionType = snakeTrace.sessionType;
loadFlag = 0; %do you want to load and save the data individually outside of ImBatAnalyze
saveFlag = 0; %do you want to load and save the data individually outside of ImBatAnalyze
meanSmooth = 30; %number of video frames to smooth the calcium mean spiking activity
out = 0; %save output variables?

traceDataAll = snakeTrace.normTraceRawPreFlightPost;
meanTraceDataAll = snakeTrace.meanTracePreFlightPost;
traceDataPre = snakeTrace.meanTracePre;
traceDataFlight = snakeTrace.meanTraceFlight;
traceDataPost = snakeTrace.meanTracePost;

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

%%plot the pre/post/and flight traces for each cell by their transitions
if plotEachFlag == 1
    plotFlightvsCells_transitions = figure('units','normalized','outerposition',[0 0 1 1]);
    %for each flight in the subgroup AtoB
    for flight_i = 1:length(flightTransitions.AtoB)
        
        p1 = subplot(3,6,1); %plot trajectories pre-flight
        plot(flightPaths.trajectories_continuous(1,flightPaths.flight_starts_idx(flightTransitions.AtoB(flight_i))-snakeTrace.preFlightPadSpeed:flightPaths.flight_starts_idx(flightTransitions.AtoB(flight_i))),flightPaths.trajectories_continuous(2,flightPaths.flight_starts_idx(flightTransitions.AtoB(flight_i))-snakeTrace.preFlightPadSpeed:flightPaths.flight_starts_idx(flightTransitions.AtoB(flight_i))),'LineWidth',1,'Color','b')%,jj(traj*100,:))
        hold on
        scatter(flightPaths.flight_starts_xyz(flightTransitions.AtoB(flight_i),1),flightPaths.flight_starts_xyz(flightTransitions.AtoB(flight_i),2),50,'r','filled')
        title('Pre-Flight: A to B');
        ylabel('y (cm)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt/10);
        xlabel('x (cm)');
        xt = get(gca,'XTick');
        set(gca,'XTick',xt,'XTickLabel',xt/10);
        
        p2 = subplot(3,6,2); %plot trajectories during flight
        plot(flightPaths.trajectories_continuous(1,flightPaths.flight_starts_idx(flightTransitions.AtoB(flight_i)):flightPaths.flight_ends_idx(flightTransitions.AtoB(flight_i))),flightPaths.trajectories_continuous(2,flightPaths.flight_starts_idx(flightTransitions.AtoB(flight_i)):flightPaths.flight_ends_idx(flightTransitions.AtoB(flight_i))),'LineWidth',1,'Color','k')%,jj(traj*100,:))
        %scatter3(fstartxyz(nf,1),fstartxyz(nf,2),fstartxyz(nf,3),25,'r','filled')
        hold on
        scatter(flightPaths.flight_ends_xyz(flightTransitions.AtoB(flight_i),1),flightPaths.flight_ends_xyz(flightTransitions.AtoB(flight_i),2),50,'k','filled')
        title('During Flight: A to B');
        ylabel('y (cm)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt/10);
        xlabel('x (cm)');
        xt = get(gca,'XTick');
        set(gca,'XTick',xt,'XTickLabel',xt/10);
        
        p3 = subplot(3,6,3); %plot trajectories post-flight A to B
        plot(flightPaths.trajectories_continuous(1,flightPaths.flight_ends_idx(flightTransitions.AtoB(flight_i)):flightPaths.flight_ends_idx(flightTransitions.AtoB(flight_i))+snakeTrace.postFlightPadSpeed),flightPaths.trajectories_continuous(2,flightPaths.flight_ends_idx(flightTransitions.AtoB(flight_i)):flightPaths.flight_ends_idx(flightTransitions.AtoB(flight_i))+snakeTrace.postFlightPadSpeed),'LineWidth',1,'Color','m')%,jj(traj*100,:))
        %scatter3(fstartxyz(nf,1),fstartxyz(nf,2),fstartxyz(nf,3),25,'r','filled')
        hold on
        scatter(flightPaths.flight_ends_xyz(flightTransitions.AtoB(flight_i),1),flightPaths.flight_ends_xyz(flightTransitions.AtoB(flight_i),2),50,'k','filled')
        title('Post-Flight: A to B');
        ylabel('y (cm)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt/10);
        xlabel('x (cm)');
        xt = get(gca,'XTick');
        set(gca,'XTick',xt,'XTickLabel',xt/10);
        
        p7 = subplot(3,6,7); %plot speed pre flight A to B
        %plot(1:length(snakeTrace.smoothSpeedPre{p}),snakeTrace.smoothSpeedPre{p},'k');
        %hold on
        plot(1:length(snakeTrace.smoothSpeedRawPre{1}(find([flightPaths.clusterIndex{:}]==flightTransitions.AtoB(flight_i)),:)),snakeTrace.smoothSpeedRawPre{1}(find([flightPaths.clusterIndex{:}]==flightTransitions.AtoB(flight_i)),:))
        hold on
        ylabel('velocity (cm/s)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
        set(gca,'xticklabel',{[]});
        ylim([0 6]);
        %xlim([1 length(snakeTrace.smoothSpeedRawPre{1}(flightPaths.clusterIndex{1}(flight_i,:)))]);
        
        p8 = subplot(3,6,8); %plot speed pre flight A to B
        %plot(1:length(snakeTrace.smoothSpeedPre{p}),snakeTrace.smoothSpeedPre{p},'k');
        %hold on
        plot(1:length(snakeTrace.smoothSpeedRawFlight{1}(find([flightPaths.clusterIndex{:}]==flightTransitions.AtoB(flight_i)),:)),snakeTrace.smoothSpeedRawFlight{1}(find([flightPaths.clusterIndex{:}]==flightTransitions.AtoB(flight_i)),:))
        hold on
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
        set(gca,'xticklabel',{[]});
        ylim([0 6]);
        
        %xlim([1 length(snakeTrace.smoothSpeedRawFlight{1}(flightPaths.clusterIndex{1}(flight_i,:)))]);
        
        p9 = subplot(3,6,9); %plot speed pre flight A to B
        %plot(1:length(snakeTrace.smoothSpeedPre{p}),snakeTrace.smoothSpeedPre{p},'k');
        %hold on
        plot(1:length(snakeTrace.smoothSpeedRawPost{1}(find([flightPaths.clusterIndex{:}]==flightTransitions.AtoB(flight_i)),:)),snakeTrace.smoothSpeedRawPost{1}(find([flightPaths.clusterIndex{:}]==flightTransitions.AtoB(flight_i)),:))
        hold on
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
        set(gca,'xticklabel',{[]});
        ylim([0 6]);
        %xlim([1 length(snakeTrace.smoothSpeedRawPost{1}(flightPaths.clusterIndex{1}(flight_i,:)))]);
    end
    
    %for each flight in the subgroup B to A
    for flight_i = 1:length(flightTransitions.BtoA)
        p4 = subplot(3,6,4); %plot trajectories pre-flight B to A
        plot(flightPaths.trajectories_continuous(1,flightPaths.flight_starts_idx(flightTransitions.BtoA(flight_i))-snakeTrace.preFlightPadSpeed:flightPaths.flight_starts_idx(flightTransitions.BtoA(flight_i))),flightPaths.trajectories_continuous(2,flightPaths.flight_starts_idx(flightTransitions.BtoA(flight_i))-snakeTrace.preFlightPadSpeed:flightPaths.flight_starts_idx(flightTransitions.BtoA(flight_i))),'LineWidth',1,'Color','b')%,jj(traj*100,:))
        hold on
        scatter(flightPaths.flight_starts_xyz(flightTransitions.BtoA(flight_i),1),flightPaths.flight_starts_xyz(flightTransitions.BtoA(flight_i),2),50,'r','filled')
        title('Pre-Flight: B to A');
        ylabel('y (cm)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt/10);
        xlabel('x (cm)');
        xt = get(gca,'XTick');
        set(gca,'XTick',xt,'XTickLabel',xt/10);
        
        p5 = subplot(3,6,5); %plot trajectories during flight B to A
        plot(flightPaths.trajectories_continuous(1,flightPaths.flight_starts_idx(flightTransitions.BtoA(flight_i)):flightPaths.flight_ends_idx(flightTransitions.BtoA(flight_i))),flightPaths.trajectories_continuous(2,flightPaths.flight_starts_idx(flightTransitions.BtoA(flight_i)):flightPaths.flight_ends_idx(flightTransitions.BtoA(flight_i))),'LineWidth',1,'Color','k')%,jj(traj*100,:))
        %scatter3(fstartxyz(nf,1),fstartxyz(nf,2),fstartxyz(nf,3),25,'r','filled')
        hold on
        scatter(flightPaths.flight_ends_xyz(flightTransitions.BtoA(flight_i),1),flightPaths.flight_ends_xyz(flightTransitions.BtoA(flight_i),2),50,'k','filled')
        title('During Flight: B to A');
        ylabel('y (cm)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt/10);
        xlabel('x (cm)');
        xt = get(gca,'XTick');
        set(gca,'XTick',xt,'XTickLabel',xt/10);
        
        p6 = subplot(3,6,6); %plot trajectories post flight B to A
        plot(flightPaths.trajectories_continuous(1,flightPaths.flight_ends_idx(flightTransitions.BtoA(flight_i)):flightPaths.flight_ends_idx(flightTransitions.BtoA(flight_i))+snakeTrace.postFlightPadSpeed),flightPaths.trajectories_continuous(2,flightPaths.flight_ends_idx(flightTransitions.BtoA(flight_i)):flightPaths.flight_ends_idx(flightTransitions.BtoA(flight_i))+snakeTrace.postFlightPadSpeed),'LineWidth',1,'Color','m')%,jj(traj*100,:))
        %scatter3(fstartxyz(nf,1),fstartxyz(nf,2),fstartxyz(nf,3),25,'r','filled')
        hold on
        scatter(flightPaths.flight_ends_xyz(flightTransitions.BtoA(flight_i),1),flightPaths.flight_ends_xyz(flightTransitions.BtoA(flight_i),2),50,'k','filled')
        title('Post-Flight: B to A');
        ylabel('y (cm)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt/10);
        xlabel('x (cm)');
        xt = get(gca,'XTick');
        set(gca,'XTick',xt,'XTickLabel',xt/10);
        
        p10 = subplot(3,6,10); %plot speed pre flight A to B
        %plot(1:length(snakeTrace.smoothSpeedPre{p}),snakeTrace.smoothSpeedPre{p},'k');
        %hold on
        plot(1:length(snakeTrace.smoothSpeedRawPre{2}(find([flightPaths.clusterIndex{:}]==flightTransitions.BtoA(flight_i))-length(snakeTrace.smoothSpeedRawPre{1}(:,1)),:)),snakeTrace.smoothSpeedRawPre{2}(find([flightPaths.clusterIndex{:}]==flightTransitions.BtoA(flight_i))-length(snakeTrace.smoothSpeedRawPre{1}(:,1)),:))
        hold on
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
        set(gca,'xticklabel',{[]});
        ylim([0 6]);
        %xlim([1 length(snakeTrace.smoothSpeedRawPre{1}(flightPaths.clusterIndex{1}(flight_i,:)))]);
        
        p11 = subplot(3,6,11); %plot speed pre flight A to B
        %plot(1:length(snakeTrace.smoothSpeedPre{p}),snakeTrace.smoothSpeedPre{p},'k');
        %hold on
        plot(1:length(snakeTrace.smoothSpeedRawFlight{2}(find([flightPaths.clusterIndex{:}]==flightTransitions.BtoA(flight_i))-length(snakeTrace.smoothSpeedRawFlight{1}(:,1)),:)),snakeTrace.smoothSpeedRawFlight{2}(find([flightPaths.clusterIndex{:}]==flightTransitions.BtoA(flight_i))-length(snakeTrace.smoothSpeedRawFlight{1}(:,1)),:))
        hold on
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
        set(gca,'xticklabel',{[]});
        ylim([0 6]);
        
        %xlim([1 length(snakeTrace.smoothSpeedRawFlight{1}(flightPaths.clusterIndex{1}(flight_i,:)))]);
        
        p12 = subplot(3,6,12); %plot speed pre flight A to B
        %plot(1:length(snakeTrace.smoothSpeedPre{p}),snakeTrace.smoothSpeedPre{p},'k');
        %hold on
        plot(1:length(snakeTrace.smoothSpeedRawPost{2}(find([flightPaths.clusterIndex{:}]==flightTransitions.BtoA(flight_i))-length(snakeTrace.smoothSpeedRawPost{1}(:,1)),:)),snakeTrace.smoothSpeedRawPost{2}(find([flightPaths.clusterIndex{:}]==flightTransitions.BtoA(flight_i))-length(snakeTrace.smoothSpeedRawPost{1}(:,1)),:))
        hold on
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
        set(gca,'xticklabel',{[]});
        ylim([0 6]);
        %xlim([1 length(snakeTrace.smoothSpeedRawPost{1}(flightPaths.clusterIndex{1}(flight_i,:)))]);
        
    end
    
    for cell_i = 1:10%length(snakeTrace.smoothTraceRawFlight{1}(:,1,:))
        for flight_i = 1:length(flightTransitions.AtoB)
            p13 = subplot(3,6,13); %plot neuron traces for each pre-flight A to B
            plot(1:length(traceDataPre{1}(1,:,1)),traceDataPre{1}(find([flightPaths.clusterIndex{:}]==flightTransitions.AtoB(flight_i)),:,cell_i))
            hold on
            title(['Cell # ' num2str(cell_i)]);
            ylabel('Norm df/f');
            yt = get(gca,'YTick');
            %set(gca,'YTick',yt,'YTickLabel',yt/10);
            xt = get(gca, 'XTick');
            set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
            xlabel('time (s)');
            
            p14 = subplot(3,6,14); %plot neuron traces for each pre-flight A to B
            plot(1:length(traceDataFlight{1}(1,:,1)),traceDataFlight{1}(find([flightPaths.clusterIndex{:}]==flightTransitions.AtoB(flight_i)),:,cell_i))
            hold on
            title(['Cell # ' num2str(cell_i)]);
            ylabel('Norm df/f');
            yt = get(gca,'YTick');
            %set(gca,'YTick',yt,'YTickLabel',yt/10);
            xt = get(gca, 'XTick');
            set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
            xlabel('time (s)');
            
            p15 = subplot(3,6,15); %plot neuron traces for each pre-flight A to B
            plot(1:length(traceDataPost{1}(1,:,1)),traceDataPost{1}(find([flightPaths.clusterIndex{:}]==flightTransitions.AtoB(flight_i)),:,cell_i))
            hold on
            title(['Cell # ' num2str(cell_i)]);
            ylabel('Norm df/f');
            yt = get(gca,'YTick');
            %set(gca,'YTick',yt,'YTickLabel',yt/10);
            xt = get(gca, 'XTick');
            set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
            xlabel('time (s)');
        end
        
        for flight_i = 1:length(flightTransitions.BtoA)
            p16 = subplot(3,6,16); %plot neuron traces for each pre-flight B to A
            plot(1:length(traceDataPre{2}(1,:,1)),traceDataPre{2}(find([flightPaths.clusterIndex{:}]==flightTransitions.BtoA(flight_i))-length(traceDataPre{1}(:,1)),:,cell_i))
            hold on
            title(['Cell # ' num2str(cell_i)]);
            ylabel('Norm df/f');
            yt = get(gca,'YTick');
            %set(gca,'YTick',yt,'YTickLabel',yt/10);
            xt = get(gca, 'XTick');
            set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
            xlabel('time (s)');
            
            p17 = subplot(3,6,17); %plot neuron traces for each during flight B to A
            plot(1:length(traceDataFlight{2}(1,:,1)),traceDataFlight{2}(find([flightPaths.clusterIndex{:}]==flightTransitions.BtoA(flight_i))-length(traceDataFlight{1}(:,1)),:,cell_i))
            hold on
            title(['Cell # ' num2str(cell_i)]);
            ylabel('Norm df/f');
            yt = get(gca,'YTick');
            %set(gca,'YTick',yt,'YTickLabel',yt/10);
            xt = get(gca, 'XTick');
            set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
            xlabel('time (s)');
            
            p18 = subplot(3,6,18); %plot neuron traces for each post-flight B to A
            plot(1:length(traceDataPost{2}(1,:,1)),traceDataPost{2}(find([flightPaths.clusterIndex{:}]==flightTransitions.BtoA(flight_i))-length(traceDataPost{1}(:,1)),:,cell_i))
            hold on
            title(['Cell # ' num2str(cell_i)]);
            ylabel('Norm df/f');
            yt = get(gca,'YTick');
            %set(gca,'YTick',yt,'YTickLabel',yt/10);
            xt = get(gca, 'XTick');
            set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
            xlabel('time (s)');
        end
        
        pause
        cla(p13)
        cla(p14)
        cla(p15)
        cla(p16)
        cla(p17)
        cla(p18)
    end
end

%% plot flight vs cell traces for each cluster for each cell with the pre/flight/post concatenated together in 1 plot
if plotPreFlightPostFlag == 1
    
    plotPreFlightPostvsCells_transitions = figure('units','normalized','outerposition',[0 0 1 1]);
    %for each flight in the subgroup AtoB
    for flight_i = 1:length(flightTransitions.AtoB)
        
        p1 = subplot(5,2,1); %plot trajectories pre/flight/post for A to B
        plot(flightPaths.trajectories_continuous(1,flightPaths.flight_starts_idx(flightTransitions.AtoB(flight_i))-snakeTrace.preFlightPadSpeed:flightPaths.flight_ends_idx(flightTransitions.AtoB(flight_i))+snakeTrace.postFlightPadSpeed),flightPaths.trajectories_continuous(2,flightPaths.flight_starts_idx(flightTransitions.AtoB(flight_i))-snakeTrace.preFlightPadSpeed:flightPaths.flight_ends_idx(flightTransitions.AtoB(flight_i))+snakeTrace.postFlightPadSpeed),'LineWidth',1,'Color','b')%,jj(traj*100,:))
        hold on
        scatter(flightPaths.flight_starts_xyz(flightTransitions.AtoB(flight_i),1),flightPaths.flight_starts_xyz(flightTransitions.AtoB(flight_i),2),50,'r','filled')
        scatter(flightPaths.flight_ends_xyz(flightTransitions.AtoB(flight_i),1),flightPaths.flight_ends_xyz(flightTransitions.AtoB(flight_i),2),50,'k','filled')
        title('A to B (r=start,b=ends)');
        ylabel('y (cm)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt/10);
        xlabel('x (cm)');
        xt = get(gca,'XTick');
        set(gca,'XTick',xt,'XTickLabel',xt/10);
        
        p2 = subplot(5,2,3); %plot speed pre/flight/post A to B
        plot(1:length(snakeTrace.smoothSpeedRawPreFlightPost{1}(find([flightPaths.clusterIndex{:}]==flightTransitions.AtoB(flight_i)),:)),snakeTrace.smoothSpeedRawPreFlightPost{1}(find([flightPaths.clusterIndex{:}]==flightTransitions.AtoB(flight_i)),:))
        hold on
        ylabel('velocity (cm/s)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
        ylim([0 6]);
        xlim([1 length(snakeTrace.smoothSpeedRawPreFlightPost{1}(find([flightPaths.clusterIndex{:}]==flightTransitions.AtoB(flight_i)),:))]);
    end
    
    %for each flight in the subgroup BtoA
    for flight_i = 1:length(flightTransitions.BtoA)
        
        p3 = subplot(5,2,2); %plot trajectories pre/flight/post for A to A
        plot(flightPaths.trajectories_continuous(1,flightPaths.flight_starts_idx(flightTransitions.BtoA(flight_i))-snakeTrace.preFlightPadSpeed:flightPaths.flight_ends_idx(flightTransitions.BtoA(flight_i))+snakeTrace.postFlightPadSpeed),flightPaths.trajectories_continuous(2,flightPaths.flight_starts_idx(flightTransitions.BtoA(flight_i))-snakeTrace.preFlightPadSpeed:flightPaths.flight_ends_idx(flightTransitions.BtoA(flight_i))+snakeTrace.postFlightPadSpeed),'LineWidth',1,'Color','m')%,jj(traj*100,:))
        hold on
        scatter(flightPaths.flight_starts_xyz(flightTransitions.BtoA(flight_i),1),flightPaths.flight_starts_xyz(flightTransitions.BtoA(flight_i),2),50,'r','filled')
        scatter(flightPaths.flight_ends_xyz(flightTransitions.BtoA(flight_i),1),flightPaths.flight_ends_xyz(flightTransitions.BtoA(flight_i),2),50,'k','filled')
        title('B to A (r=start,b=ends)');
        ylabel('y (cm)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt/10);
        xlabel('x (cm)');
        xt = get(gca,'XTick');
        set(gca,'XTick',xt,'XTickLabel',xt/10);
        
        p4 = subplot(5,2,4); %plot speed pre/flight/post B to A
        plot(1:length(snakeTrace.smoothSpeedRawPreFlightPost{2}(find([flightPaths.clusterIndex{:}]==flightTransitions.BtoA(flight_i))-length(snakeTrace.smoothSpeedRawPreFlightPost{1}(:,1)),:)),snakeTrace.smoothSpeedRawPreFlightPost{2}(find([flightPaths.clusterIndex{:}]==flightTransitions.BtoA(flight_i))-length(snakeTrace.smoothSpeedRawPreFlightPost{1}(:,1)),:))
        hold on
        ylabel('velocity (cm/s)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
        ylim([0 6]);
        xlim([1 length(snakeTrace.smoothSpeedRawPreFlightPost{2}(find([flightPaths.clusterIndex{:}]==flightTransitions.BtoA(flight_i))-length(snakeTrace.smoothSpeedRawPreFlightPost{1}(:,1)),:))]);
    end
    
    for cell_i = 1:length(snakeTrace.smoothTraceRawPreFlightPost{1}(:,1,:))
        for flight_i = 1:length(flightTransitions.AtoB)
            p5 = subplot(5,2,[5 7]); %plot neuron traces for each pre-flight A to B
            normTraceRawAllAtoB(flight_i,:) = traceDataAll{1}(find([flightPaths.clusterIndex{:}]==flightTransitions.AtoB(flight_i)),:,cell_i);
            plot(1:length(traceDataAll{1}(1,:,1)),traceDataAll{1}(find([flightPaths.clusterIndex{:}]==flightTransitions.AtoB(flight_i)),:,cell_i)+flight_i*2)
            hold on
            title(['Cell # ' num2str(cell_i)]);
            ylabel('Norm df/f');
            yt = get(gca,'YTick');
            %set(gca,'YTick',yt,'YTickLabel',yt/10);
            set(gca,'xticklabel',{[]});
            xlim([1 length(traceDataAll{1}(1,:,1))]);
        end
        p7 = subplot(5,2,9);
        imagesc(normTraceRawAllAtoB)
        colormap('hot');
        xt = get(gca, 'XTick');
        set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
        xlabel('time (s)');
        ylabel('flight trials');
        
        for flight_i = 1:length(flightTransitions.BtoA)
            p6 = subplot(5,2,[6 8]); %plot neuron traces for each pre-flight B to A
            normTraceRawAllBtoA(flight_i,:) = traceDataAll{2}(find([flightPaths.clusterIndex{:}]==flightTransitions.BtoA(flight_i))-length(traceDataAll{1}(:,1,1)),:,cell_i);
            plot(1:length(traceDataAll{2}(1,:,1)),traceDataAll{2}(find([flightPaths.clusterIndex{:}]==flightTransitions.BtoA(flight_i))-length(traceDataAll{1}(:,1,1)),:,cell_i)+flight_i*2)
            hold on
            title(['Cell # ' num2str(cell_i)]);
            ylabel('Norm df/f');
            yt = get(gca,'YTick');
            %set(gca,'YTick',yt,'YTickLabel',yt/10);
            set(gca,'xticklabel',{[]});
            xlim([1 length(traceDataAll{2}(1,:,1))]);
        end
        p8 = subplot(5,2,10);
        imagesc(normTraceRawAllBtoA)
        colormap('hot');
        xt = get(gca, 'XTick');
        set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
        xlabel('time (s)');
        ylabel('flight trials');
        linkaxes([p5, p7], 'x');
        linkaxes([p6, p8], 'x');
        sgtitle(['Flights aligned to A or B: ' snakeTrace.batName ' ' snakeTrace.dateSesh ' ' snakeTrace.sessionType]);
        if saveFlag ==1
            % Save 'each cell' as jpg and fig files..
            set(findall(gcf,'-property','FontSize'),'FontSize',20);
            saveas(gcf,[pwd '\' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_transAtoB_cell' num2str(cell_i) '.tif']);
            savefig(gcf,[pwd '\' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_transAtoB_cell' num2str(cell_i) '.fig']);
            saveas(gcf,[pwd '\' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_transAtoB_cell' num2str(cell_i) '.svg']);
        else
            pause
        end
        
        cla(p5)
        cla(p6)
        cla(p7)
        cla(p8)
    end
end

%% plot flight vs spike Means for each cluster for each cell with the pre/flight/post concatenated together in 1 plot
concAtoB = [];
concBtoA = [];
concSemAtoB = [];
concSemBtoA = [];
%normalize all the A to B data
for cell_i = 1:length(meanTraceDataAll{1}(:,1))
    concAtoB = [concAtoB meanTraceDataAll{1}(cell_i,:)]; %concatenate all the data, zscore, then split it up
    concSemAtoB = [ concSemAtoB snakeTrace.semTracePreFlightPost{1}(cell_i,:)];
end
normConcAtoB = zscore(concAtoB);
normConcAtoB = normConcAtoB - min(normConcAtoB);
normConcSemAtoB = zscore(concSemAtoB);
normConcSemAtoB = normConcSemAtoB - min(normConcSemAtoB);
for cell_i = goodCellIdx.lowCorrIndex%1:length(meanTraceDataAll{1}(:,1))
    normAtoB(cell_i,:) = concAtoB(1+(length(meanTraceDataAll{1}(1,:))*(cell_i-1)):length(meanTraceDataAll{1}(1,:))*cell_i);
    normSemAtoB(cell_i,:) = concSemAtoB(1+(length(meanTraceDataAll{1}(1,:))*(cell_i-1)):length(meanTraceDataAll{1}(1,:))*cell_i);
end
%normalize all the B to A data
for cell_i = 1:length(meanTraceDataAll{2}(:,1))
    concBtoA = [concBtoA meanTraceDataAll{2}(cell_i,:)]; %concatenate all the data, zscore, then split it up
    concSemBtoA = [ concSemBtoA snakeTrace.semTracePreFlightPost{2}(cell_i,:)];
end
normConcBtoA = zscore(concBtoA);
normConcBtoA = normConcBtoA - min(normConcBtoA);
normConcSemBtoA = zscore(concSemBtoA);
normConcSemBtoA = normConcSemBtoA - min(normConcSemBtoA);
for cell_i = goodCellIdx.lowCorrIndex%1:length(meanTraceDataAll{2}(:,1))
    normBtoA(cell_i,:) = concBtoA(1+(length(meanTraceDataAll{2}(1,:))*(cell_i-1)):length(meanTraceDataAll{2}(1,:))*cell_i);
    normSemBtoA(cell_i,:) = concSemBtoA(1+(length(meanTraceDataAll{2}(1,:))*(cell_i-1)):length(meanTraceDataAll{2}(1,:))*cell_i);
end
%smooth the data
for cell_i = goodCellIdx.lowCorrIndex%1:length(snakeTrace.smoothTraceRawPreFlightPost{1}(:,1,:))
    %smooth data and calculate full-width half max + offset from flight starts
    %smoothAtoB(cell_i,:) = smooth(meanTraceDataAll{1}(cell_i,:),meanSmooth); %smooth the data
    smoothAtoB(cell_i,:) = smooth(normAtoB(cell_i,:),meanSmooth); %smooth the normalized data
    smoothSDAtoB(cell_i,:) = smooth(snakeTrace.sdTracePreFlightPost{1}(cell_i,:),meanSmooth); %smooth the data
    smoothSEMAtoB(cell_i,:) = smooth(snakeTrace.semTracePreFlightPost{1}(cell_i,:),meanSmooth); %smooth the data
    %smoothSEMAtoB(cell_i,:) = smooth(normSemAtoB,meanSmooth);
    [maxAtoB(cell_i),indMaxAtoB(cell_i)] = max(smoothAtoB(cell_i,5:end-5));
    tStartAtoMax(cell_i) = indMaxAtoB(cell_i) - length(traceDataPre{1}(cell_i,:)); %find offset from start of flight A to peak
    halfMaxAtoB(cell_i) = maxAtoB(cell_i)/2; %find halfmax
    try
        indRiseAtoB(cell_i) = find(smoothAtoB(cell_i,5:indMaxAtoB(cell_i))<=halfMaxAtoB(cell_i),1,'last') + 5; %find where data goes above half max (-5 because we start with the 6th element)
    catch
        indRiseAtoB(cell_i) = 1; %in case the activity never drops below the fwhm by beginning of sample
    end
    try
        indFallAtoB(cell_i) = find(smoothAtoB(cell_i,indMaxAtoB(cell_i):end-5)<=halfMaxAtoB(cell_i),1,'first')+indMaxAtoB(cell_i)-1; %find where data drops below half max (-1 because you need to look 1 back from the first below half max)
    catch
        indFallAtoB(cell_i) = length(smoothAtoB(cell_i,:)); %in case the activity never drops below the fwhm by end of sample
    end
    if indFallAtoB(cell_i) < indRiseAtoB(cell_i)
        indFallAtoB(cell_i) = length(smoothAtoB(cell_i,:));
    end
    fwhmAtoB(cell_i) = indRiseAtoB(cell_i) - indFallAtoB(cell_i); %fwhm by subtracting the indices
    
    %plot neuron traces for each pre-flight B to A
    %smooth data and calculate full-width half max + offset from flight starts
    %smoothBtoA(cell_i,:) = smooth(meanTraceDataAll{2}(cell_i,:),meanSmooth);
    smoothBtoA(cell_i,:) = smooth(normBtoA(cell_i,:),meanSmooth); %smooth the normalized data
    smoothSDBtoA(cell_i,:) = smooth(snakeTrace.sdTracePreFlightPost{2}(cell_i,:),meanSmooth);
    smoothSEMBtoA(cell_i,:) = smooth(snakeTrace.semTracePreFlightPost{2}(cell_i,:),meanSmooth);
    %smoothSEMBtoA(cell_i,:) = smooth(normSemBtoA,meanSmooth); %smooth the normalized SEM
    [maxBtoA(cell_i),indMaxBtoA(cell_i)] = max(smoothBtoA(cell_i,5:end-5));
    tStartBtoMax(cell_i) = indMaxBtoA(cell_i) - length(traceDataPre{2}(cell_i,:)); %find offset from start of flight A to peak
    halfMaxBtoA(cell_i) = maxBtoA(cell_i)/2; %find halfmax
    try
        indRiseBtoA(cell_i) = find(smoothBtoA(cell_i,5:indMaxBtoA(cell_i))<=halfMaxBtoA(cell_i),1,'last')+5; %find where data goes above half max
    catch
        indRiseBtoA(cell_i) = 1; %in case the activity never drops below the fwhm by beginning of sample
    end
    try
        indFallBtoA(cell_i) = find(smoothBtoA(cell_i,indMaxBtoA(cell_i):end-5)<=halfMaxBtoA(cell_i),1,'first')+indMaxBtoA(cell_i)-1; %find where data drops below half max
    catch
        indFallBtoA(cell_i) = length(smoothBtoA(cell_i,:)); %in case the activity never drops below the fwhm by end of sample
    end
    if indFallBtoA(cell_i) < indRiseBtoA(cell_i)
        indFallBtoA(cell_i) = length(smoothBtoA(cell_i,:));
    end
    fwhmBtoA(cell_i) = indRiseBtoA(cell_i) - indFallBtoA(cell_i); %fwhm by subtracting the indices
end
%     normAtoB = zscore(meanTraceDataAll{1},0,2);
%     normAtoB = normAtoB - min(normAtoB);
%     normBtoA = zscore(meanTraceDataAll{2},0,2);
%     normBtoA = normBtoA - min(normBtoA);

if plotMeanSpikeFlag == 1
    
    plotMeanSpike_transitions = figure('units','normalized','outerposition',[0 0 1 1]);
    %for each flight in the subgroup AtoB
    for flight_i = 1:length(flightTransitions.AtoB)
        
        p1 = subplot(5,2,1); %plot trajectories pre/flight/post for A to B
        plot(flightPaths.trajectories_continuous(1,flightPaths.flight_starts_idx(flightTransitions.AtoB(flight_i))-snakeTrace.preFlightPadSpeed:flightPaths.flight_ends_idx(flightTransitions.AtoB(flight_i))+snakeTrace.postFlightPadSpeed),flightPaths.trajectories_continuous(2,flightPaths.flight_starts_idx(flightTransitions.AtoB(flight_i))-snakeTrace.preFlightPadSpeed:flightPaths.flight_ends_idx(flightTransitions.AtoB(flight_i))+snakeTrace.postFlightPadSpeed),'LineWidth',1,'Color','b')%,jj(traj*100,:))
        hold on
        scatter(flightPaths.flight_starts_xyz(flightTransitions.AtoB(flight_i),1),flightPaths.flight_starts_xyz(flightTransitions.AtoB(flight_i),2),50,'r','filled')
        scatter(flightPaths.flight_ends_xyz(flightTransitions.AtoB(flight_i),1),flightPaths.flight_ends_xyz(flightTransitions.AtoB(flight_i),2),50,'k','filled')
        title('A to B (r=start,b=ends)');
        ylabel('y (cm)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt/10);
        xlabel('x (cm)');
        xt = get(gca,'XTick');
        set(gca,'XTick',xt,'XTickLabel',xt/10);
        
        p2 = subplot(5,2,3); %plot speed pre/flight/post A to B
        plot(1:length(snakeTrace.smoothSpeedRawPreFlightPost{1}(find([flightPaths.clusterIndex{:}]==flightTransitions.AtoB(flight_i)),:)),snakeTrace.smoothSpeedRawPreFlightPost{1}(find([flightPaths.clusterIndex{:}]==flightTransitions.AtoB(flight_i)),:))
        hold on
        ylabel('velocity (cm/s)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
        ylim([0 6]);
        xlim([1 length(snakeTrace.smoothSpeedRawPreFlightPost{1}(find([flightPaths.clusterIndex{:}]==flightTransitions.AtoB(flight_i)),:))]);
    end
    
    %for each flight in the subgroup BtoA
    for flight_i = 1:length(flightTransitions.BtoA)
        
        p3 = subplot(5,2,2); %plot trajectories pre/flight/post for B to A
        plot(flightPaths.trajectories_continuous(1,flightPaths.flight_starts_idx(flightTransitions.BtoA(flight_i))-snakeTrace.preFlightPadSpeed:flightPaths.flight_ends_idx(flightTransitions.BtoA(flight_i))+snakeTrace.postFlightPadSpeed),flightPaths.trajectories_continuous(2,flightPaths.flight_starts_idx(flightTransitions.BtoA(flight_i))-snakeTrace.preFlightPadSpeed:flightPaths.flight_ends_idx(flightTransitions.BtoA(flight_i))+snakeTrace.postFlightPadSpeed),'LineWidth',1,'Color','m')%,jj(traj*100,:))
        hold on
        scatter(flightPaths.flight_starts_xyz(flightTransitions.BtoA(flight_i),1),flightPaths.flight_starts_xyz(flightTransitions.BtoA(flight_i),2),50,'r','filled')
        scatter(flightPaths.flight_ends_xyz(flightTransitions.BtoA(flight_i),1),flightPaths.flight_ends_xyz(flightTransitions.BtoA(flight_i),2),50,'k','filled')
        title('B to A (r=start,b=ends)');
        ylabel('y (cm)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt/10);
        xlabel('x (cm)');
        xt = get(gca,'XTick');
        set(gca,'XTick',xt,'XTickLabel',xt/10);
        
        p4 = subplot(5,2,4); %plot speed pre/flight/post B to A
        plot(1:length(snakeTrace.smoothSpeedRawPreFlightPost{2}(find([flightPaths.clusterIndex{:}]==flightTransitions.BtoA(flight_i))-length(snakeTrace.smoothSpeedRawPreFlightPost{1}(:,1)),:)),snakeTrace.smoothSpeedRawPreFlightPost{2}(find([flightPaths.clusterIndex{:}]==flightTransitions.BtoA(flight_i))-length(snakeTrace.smoothSpeedRawPreFlightPost{1}(:,1)),:))
        hold on
        ylabel('velocity (cm/s)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
        ylim([0 6]);
        xlim([1 length(snakeTrace.smoothSpeedRawPreFlightPost{2}(find([flightPaths.clusterIndex{:}]==flightTransitions.BtoA(flight_i))-length(snakeTrace.smoothSpeedRawPreFlightPost{1}(:,1)),:))]);
    end
    
    %plot neural traces
    for cell_i = goodCellIdx.lowCorrIndex%1:length(snakeTrace.smoothTraceRawPreFlightPost{1}(:,1,:))
        
        %plot average spike activity for flights A to B
        p5 = subplot(5,2,[5 7]); %plot neuron traces for each pre-flight A to B
        boundedline(1:length(meanTraceDataAll{1}(1,:)),smoothAtoB(cell_i,:),smoothSEMAtoB(cell_i,:),'r');
        hold on
        plot(1:length(meanTraceDataAll{1}(1,:)),meanTraceDataAll{1}(cell_i,:));
        %plot(1:length(meanTraceDataAll{1}(1,:)),smoothAtoB(cell_i,:)); %plot the smoothed trace on top
        plot(indMaxAtoB(cell_i),maxAtoB(cell_i),'o','MarkerFaceColor','r','MarkerSize',10) %plot peak on
        plot(indRiseAtoB(cell_i):indFallAtoB(cell_i),halfMaxAtoB(cell_i),'o') %plot the half width max height
        %titles
        title(['Cell # ' num2str(cell_i)]);
        ylabel('z-scored # spikes (std)');%ylabel('# Spikes');
        yt = get(gca,'YTick');
        %set(gca,'YTick',yt,'YTickLabel',yt/10);
        set(gca,'xticklabel',{[]});
        xlim([1 length(meanTraceDataAll{1}(1,:))]);
        %plot heatmap for each pre-flight A to B
        p7 = subplot(5,2,9);
        imagesc(meanTraceDataAll{1}(cell_i,:))
        colormap('hot');
        xt = get(gca, 'XTick');
        set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
        xlabel('time (s)');
        ylabel('flight trials');
        
        %plots BtoA
        p6 = subplot(5,2,[6 8]);
        boundedline(1:length(meanTraceDataAll{2}(1,:)),smoothBtoA(cell_i,:),smoothSEMBtoA(cell_i,:),'r');
        hold on
        plot(1:length(meanTraceDataAll{2}(1,:)),meanTraceDataAll{2}(cell_i,:)); %plot the mean spiking activity aligned to start of flight B
        %plot(1:length(meanTraceDataAll{2}(1,:)),smoothBtoA(cell_i,:)); %plot the smoothed trace on top
        plot(indMaxBtoA(cell_i),maxBtoA(cell_i),'o','MarkerFaceColor','r','MarkerSize',10) %plot the peak of the activity
        plot(indRiseBtoA(cell_i):indFallBtoA(cell_i),halfMaxBtoA(cell_i),'o') %plot the half width max height
        %titles
        title(['Cell # ' num2str(cell_i)]);
        ylabel('z-scored # spikes (std)');%ylabel('# Spikes');
        yt = get(gca,'YTick');
        %set(gca,'YTick',yt,'YTickLabel',yt/10);
        set(gca,'xticklabel',{[]});
        xlim([1 length(meanTraceDataAll{2}(1,:))]);
        %plot heatmap for each pre-flight B to A
        p8 = subplot(5,2,10);
        imagesc(meanTraceDataAll{2}(cell_i,:))
        colormap('hot');
        xt = get(gca, 'XTick');
        set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
        xlabel('time (s)');
        ylabel('flight trials');
        linkaxes([p5, p7], 'x');
        linkaxes([p6, p8], 'x');
        sgtitle(['Flights aligned to A or B: ' snakeTrace.batName ' ' snakeTrace.dateSesh ' ' snakeTrace.sessionType]);
        
        if saveFlag ==1
            % Save 'each cell' as jpg and fig files..
            set(findall(gcf,'-property','FontSize'),'FontSize',20);
            saveas(gcf,[pwd '/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_transAtoB_sPeaks_cell' num2str(cell_i) '.tif']);
            savefig(gcf,[pwd '/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_transAtoB_sPeaks_cell' num2str(cell_i) '.fig']);
            saveas(gcf, [pwd '/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_transAtoB_sPeaks_cell' num2str(cell_i) '.svg']);
        else
            pause
        end
        
        cla(p5)
        cla(p6)
        cla(p7)
        cla(p8)
    end
end

traceAnalysis.smoothAtoB = smoothAtoB;
traceAnalysis.smoothSDAtoB = smoothSDAtoB;
traceAnalysis.smoothSEMAtoB = smoothSEMAtoB;
traceAnalysis.maxAtoB =maxAtoB;
traceAnalysis.indMaxAtoB = indMaxAtoB;
traceAnalysis.tStartAtoMax = tStartAtoMax;
traceAnalysis.halfMaxAtoB = halfMaxAtoB;
traceAnalysis.indRiseAtoB = indRiseAtoB;
traceAnalysis.indFallAtoB = indFallAtoB;
traceAnalysis.fwhmAtoB = fwhmAtoB;
traceAnalysis.smoothBtoA = smoothBtoA;
traceAnalysis.smoothSDAtoB = smoothSDBtoA;
traceAnalysis.smoothSEMAtoB = smoothSEMBtoA;
traceAnalysis.maxBtoA = maxBtoA;
traceAnalysis.indMaxBtoA = indMaxBtoA;
traceAnalysis.tStartBtoMax = tStartBtoMax;
traceAnalysis.halfMaxBtoA = halfMaxBtoA;
traceAnalysis.indRiseBtoA = indRiseBtoA;
traceAnalysis.indFallBtoA = indFallBtoA;
traceAnalysis.fwhmBtoA = fwhmBtoA;
if saveFlag ==1
    save([pwd '/' batName '_' dateSesh '_' sessionType '_traceAnalysis.mat'],'traceAnalysis');
end
%% sort and plot the cells according to their peak from flight time
for cell_i = goodCellIdx.lowCorrIndex
    traceNormAtoB(cell_i,:) = zscore(traceAnalysis.smoothAtoB(cell_i,:));
    traceNormAtoB(cell_i,:) = traceNormAtoB(cell_i,:) - min(traceNormAtoB(cell_i,:));
    traceNormBtoA(cell_i,:) = zscore(traceAnalysis.smoothBtoA(cell_i,:));
    traceNormBtoA(cell_i,:) = traceNormBtoA(cell_i,:) - min(traceNormBtoA(cell_i,:));
    
    [~,maxSmoothBtoA(cell_i,1)] = max(traceAnalysis.smoothBtoA(cell_i,:));
    [~,maxSmoothAtoB(cell_i,1)] = max(traceAnalysis.smoothAtoB(cell_i,:));
    %[~,maxSmoothAtoB(cell_i,1)] = max(traceNormAtoB(cell_i,:));
    %[~,maxSmoothBtoA(cell_i,1)] = max(traceNormBtoA(cell_i,:));
    
end
[BTraceAtoB,ITraceAtoB] = sort(maxSmoothAtoB);
[BTraceBtoA,ITraceBtoA] = sort(maxSmoothBtoA);

if plotMeanSortedFlag == 1
    
    plotSortedPeaks = figure('units','normalized','outerposition',[0 0 1 1]);
    %for each flight in the subgroup AtoB
    for flight_i = 1:length(flightTransitions.AtoB)
        
        p1 = subplot(6,2,1); %plot trajectories pre/flight/post for A to B
        plot(flightPaths.trajectories_continuous(1,flightPaths.flight_starts_idx(flightTransitions.AtoB(flight_i))-snakeTrace.preFlightPadSpeed:flightPaths.flight_ends_idx(flightTransitions.AtoB(flight_i))+snakeTrace.postFlightPadSpeed),flightPaths.trajectories_continuous(2,flightPaths.flight_starts_idx(flightTransitions.AtoB(flight_i))-snakeTrace.preFlightPadSpeed:flightPaths.flight_ends_idx(flightTransitions.AtoB(flight_i))+snakeTrace.postFlightPadSpeed),'LineWidth',1,'Color','b')%,jj(traj*100,:))
        hold on
        scatter(flightPaths.flight_starts_xyz(flightTransitions.AtoB(flight_i),1),flightPaths.flight_starts_xyz(flightTransitions.AtoB(flight_i),2),50,'r','filled')
        scatter(flightPaths.flight_ends_xyz(flightTransitions.AtoB(flight_i),1),flightPaths.flight_ends_xyz(flightTransitions.AtoB(flight_i),2),50,'k','filled')
        title('A to B (r=start,b=ends)');
        ylabel('y (cm)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt/10);
        xlabel('x (cm)');
        xt = get(gca,'XTick');
        set(gca,'XTick',xt,'XTickLabel',xt/10);
        
        p2 = subplot(6,2,3); %plot speed pre/flight/post A to B
        plot(1:length(snakeTrace.smoothSpeedRawPreFlightPost{1}(find([flightPaths.clusterIndex{:}]==flightTransitions.AtoB(flight_i)),:)),snakeTrace.smoothSpeedRawPreFlightPost{1}(find([flightPaths.clusterIndex{:}]==flightTransitions.AtoB(flight_i)),:))
        hold on
        ylabel('velocity (cm/s)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
        ylim([0 6]);
        xlim([1 length(snakeTrace.smoothSpeedRawPreFlightPost{1}(find([flightPaths.clusterIndex{:}]==flightTransitions.AtoB(flight_i)),:))]);
    end
    
    %for each flight in the subgroup BtoA
    for flight_i = 1:length(flightTransitions.BtoA)
        
        p3 = subplot(6,2,2); %plot trajectories pre/flight/post for B to A
        plot(flightPaths.trajectories_continuous(1,flightPaths.flight_starts_idx(flightTransitions.BtoA(flight_i))-snakeTrace.preFlightPadSpeed:flightPaths.flight_ends_idx(flightTransitions.BtoA(flight_i))+snakeTrace.postFlightPadSpeed),flightPaths.trajectories_continuous(2,flightPaths.flight_starts_idx(flightTransitions.BtoA(flight_i))-snakeTrace.preFlightPadSpeed:flightPaths.flight_ends_idx(flightTransitions.BtoA(flight_i))+snakeTrace.postFlightPadSpeed),'LineWidth',1,'Color','m')%,jj(traj*100,:))
        hold on
        scatter(flightPaths.flight_starts_xyz(flightTransitions.BtoA(flight_i),1),flightPaths.flight_starts_xyz(flightTransitions.BtoA(flight_i),2),50,'r','filled')
        scatter(flightPaths.flight_ends_xyz(flightTransitions.BtoA(flight_i),1),flightPaths.flight_ends_xyz(flightTransitions.BtoA(flight_i),2),50,'k','filled')
        title('B to A (r=start,b=ends)');
        ylabel('y (cm)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt/10);
        xlabel('x (cm)');
        xt = get(gca,'XTick');
        set(gca,'XTick',xt,'XTickLabel',xt/10);
        
        p4 = subplot(6,2,4); %plot speed pre/flight/post B to A
        plot(1:length(snakeTrace.smoothSpeedRawPreFlightPost{2}(find([flightPaths.clusterIndex{:}]==flightTransitions.BtoA(flight_i))-length(snakeTrace.smoothSpeedRawPreFlightPost{1}(:,1)),:)),snakeTrace.smoothSpeedRawPreFlightPost{2}(find([flightPaths.clusterIndex{:}]==flightTransitions.BtoA(flight_i))-length(snakeTrace.smoothSpeedRawPreFlightPost{1}(:,1)),:))
        hold on
        ylabel('velocity (cm/s)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
        ylim([0 6]);
        xlim([1 length(snakeTrace.smoothSpeedRawPreFlightPost{2}(find([flightPaths.clusterIndex{:}]==flightTransitions.BtoA(flight_i))-length(snakeTrace.smoothSpeedRawPreFlightPost{1}(:,1)),:))]);
    end
    %plot heat maps of the high correlated cells aligned to A or B with
    %sorting of A or B by their peaks
    p5 = subplot(6,2,[5 7])
    %imagesc(smoothAtoB(ITraceAtoB,:));
    imagesc(traceNormAtoB(ITraceAtoB,:));
    colormap(hot)
    title(['Aligned A to B: Sorted A to B']);
    xt = get(gca, 'XTick');
    set(gca,'XTick',xt,'XTickLabel',[]);
    xlabel('time (s)');
    ylabel('ROI #');
    
    p6 = subplot(6,2,[6 8])
    %imagesc(smoothBtoA(ITraceBtoA,:));
    imagesc(traceNormBtoA(ITraceBtoA,:));
    colormap(hot)
    title(['Aligned B to A: Sorted B to A']);
    xt = get(gca, 'XTick');
    set(gca,'XTick',xt,'XTickLabel',[]);
    ylabel('ROI #');
    
    p7 = subplot(6,2,[9 11])
    %imagesc(smoothBtoA(ITraceAtoB,:));
    imagesc(traceNormBtoA(ITraceAtoB,:));
    colormap(hot)
    title(['Aligned B to A: Sorted A to B']);
    xt = get(gca, 'XTick');
    set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
    xlabel('time (s)');
    ylabel('ROI #');
    
    p8 = subplot(6,2,[10 12])
    %imagesc(smoothAtoB(ITraceBtoA,:));
    imagesc(traceNormAtoB(ITraceBtoA,:));
    colormap(hot)
    title(['Aligned A to B: Sorted B to A']);
    xt = get(gca, 'XTick');
    set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
    xlabel('time (s)');
    ylabel('ROI #');
    
    if saveFlag ==1
        saveas(plotSortedPeaks, [pwd '/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_heatMap_highCorr_sortedPeaks_NORM_AvB.svg']);
        saveas(plotSortedPeaks, [pwd '/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_heatMap_highCorr_sortedPeaks_NORM_AvB.tif']);
        savefig(plotSortedPeaks, [pwd '/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_heatMap_highCorr_sortedPeaks_NORM_AvB.fig']);
    end
end

traceAnalysis.traceNormAtoB = traceNormAtoB;
traceAnalysis.maxSmoothAtoB = maxSmoothAtoB;
traceAnalysis.ITraceAtoB = ITraceAtoB;
traceAnalysis.BTraceAtoB = BTraceAtoB;
traceAnalysis.traceNormBtoA = traceNormBtoA;
traceAnalysis.maxSmoothBtoA = maxSmoothBtoA;
traceAnalysis.ITraceBtoA = ITraceBtoA;
traceAnalysis.BTraceBtoA = BTraceBtoA;

if saveFlag ==1
    save([pwd '/' batName '_' dateSesh '_' sessionType '_traceAnalysis.mat'],'traceAnalysis');
end

%% plot the stats of the height, time from flight to peak, and full width half max of peaks aligned to A or B
if histoTraceStatsFlag == 1
    
    figure; histogram(traceAnalysis.halfMaxAtoB)
    hold on; histogram(traceAnalysis.halfMaxBtoA)
    title('Height of half-max aligned to A or B')
    legend('Aligned to A','Aligned to B')
    xlabel('Height of half-max (SD)')
    if saveFlag == 1
        saveas(gcf, [pwd '/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_histo_halfMaxHeight_AvB.svg']);
        saveas(gcf, [pwd '/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_histo_halfMaxHeight_AvB.tif']);
        savefig(gcf, [pwd '/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_histo_halfMaxHeight_AvB.fig']);
    end
    
    figure; histogram(abs(traceAnalysis.fwhmAtoB))
    hold on; histogram(abs(traceAnalysis.fwhmBtoA),11)
    title('Width of halfmax aligned to A or B')
    legend('Aligned to A','Aligned to B')
    xt = get(gca, 'XTick');
    set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
    xlabel('time (s)');
    if saveFlag == 1
        saveas(gcf, [pwd '/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_histo_fwhm_AvB.svg']);
        saveas(gcf, [pwd '/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_histo_fwhm_AvB.tif']);
        savefig(gcf, [pwd '/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_histo_fwhm_AvB.fig']);
    end
    
    figure; histogram(traceAnalysis.tStartAtoMax,10)
    hold on; histogram(traceAnalysis.tStartBtoMax)
    title('Time from flight start to peak')
    legend('Aligned to A','Aligned to B')
    xt = get(gca, 'XTick');
    set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
    xlabel('time (s)');
    if saveFlag == 1
        saveas(gcf, [pwd '/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_histo_timeStartPeak_AvB.svg']);
        saveas(gcf, [pwd '/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_histo_timeStartPeak_AvB.tif']);
        savefig(gcf, [pwd '/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_histo_timeStartPeak_AvB.fig']);
    end
end

%% Plot the fwhm across time sorted by the width of the half max
if plotSortedFWHMFlag ==1
    [BfwhmAtoB,IfwhmAtoB] = sort(traceAnalysis.fwhmAtoB);
    [BfwhmBtoA,IfwhmBtoA] = sort(traceAnalysis.fwhmBtoA);
    
    plotSortedFWHM = figure('units','normalized','outerposition',[0 0 1 1]);
    %for each flight in the subgroup AtoB
    for flight_i = 1:length(flightTransitions.AtoB)
        
        p1 = subplot(6,2,1); %plot trajectories pre/flight/post for A to B
        plot(flightPaths.trajectories_continuous(1,flightPaths.flight_starts_idx(flightTransitions.AtoB(flight_i))-snakeTrace.preFlightPadSpeed:flightPaths.flight_ends_idx(flightTransitions.AtoB(flight_i))+snakeTrace.postFlightPadSpeed),flightPaths.trajectories_continuous(2,flightPaths.flight_starts_idx(flightTransitions.AtoB(flight_i))-snakeTrace.preFlightPadSpeed:flightPaths.flight_ends_idx(flightTransitions.AtoB(flight_i))+snakeTrace.postFlightPadSpeed),'LineWidth',1,'Color','b')%,jj(traj*100,:))
        hold on
        scatter(flightPaths.flight_starts_xyz(flightTransitions.AtoB(flight_i),1),flightPaths.flight_starts_xyz(flightTransitions.AtoB(flight_i),2),50,'r','filled')
        scatter(flightPaths.flight_ends_xyz(flightTransitions.AtoB(flight_i),1),flightPaths.flight_ends_xyz(flightTransitions.AtoB(flight_i),2),50,'k','filled')
        title('A to B (r=start,b=ends)');
        ylabel('y (cm)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt/10);
        xlabel('x (cm)');
        xt = get(gca,'XTick');
        set(gca,'XTick',xt,'XTickLabel',xt/10);
        
        p2 = subplot(6,2,3); %plot speed pre/flight/post A to B
        plot(1:length(snakeTrace.smoothSpeedRawPreFlightPost{1}(find([flightPaths.clusterIndex{:}]==flightTransitions.AtoB(flight_i)),:)),snakeTrace.smoothSpeedRawPreFlightPost{1}(find([flightPaths.clusterIndex{:}]==flightTransitions.AtoB(flight_i)),:))
        hold on
        ylabel('velocity (cm/s)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
        ylim([0 6]);
        xlim([1 length(snakeTrace.smoothSpeedRawPreFlightPost{1}(find([flightPaths.clusterIndex{:}]==flightTransitions.AtoB(flight_i)),:))]);
    end
    
    %for each flight in the subgroup BtoA
    for flight_i = 1:length(flightTransitions.BtoA)
        
        p3 = subplot(6,2,2); %plot trajectories pre/flight/post for B to A
        plot(flightPaths.trajectories_continuous(1,flightPaths.flight_starts_idx(flightTransitions.BtoA(flight_i))-snakeTrace.preFlightPadSpeed:flightPaths.flight_ends_idx(flightTransitions.BtoA(flight_i))+snakeTrace.postFlightPadSpeed),flightPaths.trajectories_continuous(2,flightPaths.flight_starts_idx(flightTransitions.BtoA(flight_i))-snakeTrace.preFlightPadSpeed:flightPaths.flight_ends_idx(flightTransitions.BtoA(flight_i))+snakeTrace.postFlightPadSpeed),'LineWidth',1,'Color','m')%,jj(traj*100,:))
        hold on
        scatter(flightPaths.flight_starts_xyz(flightTransitions.BtoA(flight_i),1),flightPaths.flight_starts_xyz(flightTransitions.BtoA(flight_i),2),50,'r','filled')
        scatter(flightPaths.flight_ends_xyz(flightTransitions.BtoA(flight_i),1),flightPaths.flight_ends_xyz(flightTransitions.BtoA(flight_i),2),50,'k','filled')
        title('B to A (r=start,b=ends)');
        ylabel('y (cm)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt/10);
        xlabel('x (cm)');
        xt = get(gca,'XTick');
        set(gca,'XTick',xt,'XTickLabel',xt/10);
        
        p4 = subplot(6,2,4); %plot speed pre/flight/post B to A
        plot(1:length(snakeTrace.smoothSpeedRawPreFlightPost{2}(find([flightPaths.clusterIndex{:}]==flightTransitions.BtoA(flight_i))-length(snakeTrace.smoothSpeedRawPreFlightPost{1}(:,1)),:)),snakeTrace.smoothSpeedRawPreFlightPost{2}(find([flightPaths.clusterIndex{:}]==flightTransitions.BtoA(flight_i))-length(snakeTrace.smoothSpeedRawPreFlightPost{1}(:,1)),:))
        hold on
        ylabel('velocity (cm/s)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
        ylim([0 6]);
        xlim([1 length(snakeTrace.smoothSpeedRawPreFlightPost{2}(find([flightPaths.clusterIndex{:}]==flightTransitions.BtoA(flight_i))-length(snakeTrace.smoothSpeedRawPreFlightPost{1}(:,1)),:))]);
    end
    
    %plot the fwhm aligned by flight A or B
    for cell_i = goodCellIdx.lowCorrIndex
        p5 = subplot(6,2,[5 7]);
        plot(traceAnalysis.indRiseAtoB(cell_i):traceAnalysis.indFallAtoB(cell_i),find(IfwhmAtoB==cell_i),'o')
        hold on;
        title('Aligned A to B: Sorted A to B');
        xt = get(gca, 'XTick');
        set(gca,'XTick',xt,'XTickLabel',[]);
        ylabel('ROI #');
        
        p6 = subplot(6,2,[6 8]);
        plot(traceAnalysis.indRiseBtoA(cell_i):traceAnalysis.indFallBtoA(cell_i),find(IfwhmBtoA==cell_i),'o')
        hold on;
        title('Aligned B to A: Sorted B to A');
        xt = get(gca, 'XTick');
        set(gca,'XTick',xt,'XTickLabel',[]);
        ylabel('ROI #');
        
        p7 = subplot(6,2,[9 11]);
        plot(traceAnalysis.indRiseAtoB(cell_i):traceAnalysis.indFallAtoB(cell_i),find(IfwhmBtoA==cell_i),'o')
        hold on;
        title('Aligned A to B: Sorted B to A');
        xt = get(gca, 'XTick');
        set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
        xlabel('time (s)');
        ylabel('ROI #');
        
        p8 = subplot(6,2,[10 12]);
        plot(traceAnalysis.indRiseBtoA(cell_i):traceAnalysis.indFallBtoA(cell_i),find(IfwhmAtoB==cell_i),'o')
        hold on;
        title('Aligned B to A: Sorted A to B');
        xt = get(gca, 'XTick');
        set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
        xlabel('time (s)');
        ylabel('ROI #');
        %pause
    end
    %     subplot(2,2,1);
    %     imagesc(traceAnalysis.smoothAtoB(IfwhmAtoB,:));
    %     subplot(2,2,2);
    %     imagesc(traceAnalysis.smoothBtoA(IfwhmBtoA,:));
    %     subplot(2,2,3);
    %     imagesc(traceAnalysis.smoothAtoB(IfwhmBtoA,:));
    %     subplot(2,2,4);
    %     imagesc(traceAnalysis.smoothBtoA(IfwhmAtoB,:));
    if saveFlag == 1
        saveas(gcf, [pwd '/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_SortedFWHM_AvB.svg']);
        saveas(gcf, [pwd '/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_SortedFWHM_AvB.tif']);
        savefig(gcf, [pwd '/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_SortedFWHM_AvB.fig']);
    end
    
    
end

%%
if plotTimeJitterSortFlag == 1
    
    [BdiffStartAtoB,IdiffStartAtoB] = sort(flightTiming.diffStartAtoStartB);
    [BdiffSBtoSA,IdiffSBtoSA] = sort(flightTiming.diffSBtoSA);
    
    plotSortedFWHM = figure('units','normalized','outerposition',[0 0 1 1]);
    %for each flight in the subgroup AtoB
    for flight_i = 1:length(flightTransitions.AtoB)
        
        p1 = subplot(6,2,1); %plot trajectories pre/flight/post for A to B
        plot(flightPaths.trajectories_continuous(1,flightPaths.flight_starts_idx(flightTransitions.AtoB(flight_i))-snakeTrace.preFlightPadSpeed:flightPaths.flight_ends_idx(flightTransitions.AtoB(flight_i))+snakeTrace.postFlightPadSpeed),flightPaths.trajectories_continuous(2,flightPaths.flight_starts_idx(flightTransitions.AtoB(flight_i))-snakeTrace.preFlightPadSpeed:flightPaths.flight_ends_idx(flightTransitions.AtoB(flight_i))+snakeTrace.postFlightPadSpeed),'LineWidth',1,'Color','b')%,jj(traj*100,:))
        hold on
        scatter(flightPaths.flight_starts_xyz(flightTransitions.AtoB(flight_i),1),flightPaths.flight_starts_xyz(flightTransitions.AtoB(flight_i),2),50,'r','filled')
        scatter(flightPaths.flight_ends_xyz(flightTransitions.AtoB(flight_i),1),flightPaths.flight_ends_xyz(flightTransitions.AtoB(flight_i),2),50,'k','filled')
        title('A to B (r=start,b=ends)');
        ylabel('y (cm)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt/10);
        xlabel('x (cm)');
        xt = get(gca,'XTick');
        set(gca,'XTick',xt,'XTickLabel',xt/10);
        
        p2 = subplot(6,2,3); %plot speed pre/flight/post A to B
        plot(1:length(snakeTrace.smoothSpeedRawPreFlightPost{1}(find([flightPaths.clusterIndex{:}]==flightTransitions.AtoB(flight_i)),:)),snakeTrace.smoothSpeedRawPreFlightPost{1}(find([flightPaths.clusterIndex{:}]==flightTransitions.AtoB(flight_i)),:))
        hold on
        ylabel('velocity (cm/s)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
        ylim([0 6]);
        xlim([1 length(snakeTrace.smoothSpeedRawPreFlightPost{1}(find([flightPaths.clusterIndex{:}]==flightTransitions.AtoB(flight_i)),:))]);
    end
    
    %for each flight in the subgroup BtoA
    for flight_i = 1:length(flightTransitions.BtoA)
        
        p3 = subplot(6,2,2); %plot trajectories pre/flight/post for B to A
        plot(flightPaths.trajectories_continuous(1,flightPaths.flight_starts_idx(flightTransitions.BtoA(flight_i))-snakeTrace.preFlightPadSpeed:flightPaths.flight_ends_idx(flightTransitions.BtoA(flight_i))+snakeTrace.postFlightPadSpeed),flightPaths.trajectories_continuous(2,flightPaths.flight_starts_idx(flightTransitions.BtoA(flight_i))-snakeTrace.preFlightPadSpeed:flightPaths.flight_ends_idx(flightTransitions.BtoA(flight_i))+snakeTrace.postFlightPadSpeed),'LineWidth',1,'Color','m')%,jj(traj*100,:))
        hold on
        scatter(flightPaths.flight_starts_xyz(flightTransitions.BtoA(flight_i),1),flightPaths.flight_starts_xyz(flightTransitions.BtoA(flight_i),2),50,'r','filled')
        scatter(flightPaths.flight_ends_xyz(flightTransitions.BtoA(flight_i),1),flightPaths.flight_ends_xyz(flightTransitions.BtoA(flight_i),2),50,'k','filled')
        title('B to A (r=start,b=ends)');
        ylabel('y (cm)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt/10);
        xlabel('x (cm)');
        xt = get(gca,'XTick');
        set(gca,'XTick',xt,'XTickLabel',xt/10);
        
        p4 = subplot(6,2,4); %plot speed pre/flight/post B to A
        plot(1:length(snakeTrace.smoothSpeedRawPreFlightPost{2}(find([flightPaths.clusterIndex{:}]==flightTransitions.BtoA(flight_i))-length(snakeTrace.smoothSpeedRawPreFlightPost{1}(:,1)),:)),snakeTrace.smoothSpeedRawPreFlightPost{2}(find([flightPaths.clusterIndex{:}]==flightTransitions.BtoA(flight_i))-length(snakeTrace.smoothSpeedRawPreFlightPost{1}(:,1)),:))
        hold on
        ylabel('velocity (cm/s)');
        yt = get(gca,'YTick');
        set(gca,'YTick',yt,'YTickLabel',yt*100,'xticklabel',{[]});
        ylim([0 6]);
        xlim([1 length(snakeTrace.smoothSpeedRawPreFlightPost{2}(find([flightPaths.clusterIndex{:}]==flightTransitions.BtoA(flight_i))-length(snakeTrace.smoothSpeedRawPreFlightPost{1}(:,1)),:))]);
    end
    for cell_i = goodCellIdx.lowCorrIndex
        p5 = subplot(6,2,[5 7]);
        imagesc(snakeTrace.tracePreFlightPost{1}(IdiffStartAtoB,:,cell_i));
        hold on;
        title('Aligned A to B: Sorted from start of A to start of B');
        xt = get(gca, 'XTick');
        set(gca,'XTick',xt,'XTickLabel',[]);
        ylabel('Trial #');
        colormap('hot');
        
        
        p6 = subplot(6,2,[6 8]);
                imagesc(snakeTrace.tracePreFlightPost{2}(IdiffSBtoSA,:,cell_i));
        hold on;
        title('Aligned B to A: Sorted from start of B to start of A');
        xt = get(gca, 'XTick');
        set(gca,'XTick',xt,'XTickLabel',[]);
        ylabel('Trial #');
        colormap('hot');
        
        
        p7 = subplot(6,2,[9 11]);
        imagesc(snakeTrace.tracePreFlightPost{1}(IdiffSBtoSA,:,cell_i));
        hold on;
        title('Aligned A to B: Sorted from start of B to start of A');
        xt = get(gca, 'XTick');
        set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
        xlabel('time (s)');
        ylabel('Trial #');
        colormap('hot');
        
        
        p8 = subplot(6,2,[10 12]);
        imagesc(snakeTrace.tracePreFlightPost{2}(IdiffStartAtoB,:,cell_i));
        hold on;
        title('Aligned B to A: Sorted from start of A to start of B');
        xt = get(gca, 'XTick');
        set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
        xlabel('time (s)');
        ylabel('Trial #');
        colormap('hot');
        sgtitle(['Sorting each trial from short to long interval: Cell # ' num2str(cell_i) ' : ' snakeTrace.batName ' ' snakeTrace.dateSesh ' ' snakeTrace.sessionType]);
        
        
        
        if saveFlag == 1
            saveas(gcf, [pwd '/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_heatMapTimeJitterSort_AvB-cell' num2str(cell_i) '.svg']);
            saveas(gcf, [pwd '/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_heatMapTimeJitterSort_AvB-cell' num2str(cell_i) '.tif']);
            savefig(gcf, [pwd '/' snakeTrace.batName '_' snakeTrace.dateSesh '_' snakeTrace.sessionType '_heatMapTimeJitterSort_AvB-cell' num2str(cell_i) '.fig']);
        else
            pause
        end
        cla(p5)
        cla(p6)
        cla(p7)
        cla(p8)
    end
end

