function ImBat_plotFlightvsCells_transitions(snakeTrace,flightTransitions,flightPaths,varargin)
plotEachFlag = 0;
plotPreFlightPostFlag = 0;
plotMeanSpikeFlag = 1;

batName = [];
dateSesh = [];
sessionType = [];
loadFlag = 0; %do you want to load and save the data individually outside of ImBatAnalyze
saveFlag = 0; %do you want to load and save the data individually outside of ImBatAnalyze
meanSmooth = 30; %number of video frames to smooth the calcium mean spiking activity

traceDataAll = snakeTrace.meanTracePreFlightPost;
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
    
    for cell_i = 1:length(snakeTrace.smoothTraceRawFlight{1}(:,1,:))
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
            plot(1:length(traceDataPre{2}(1,:,1)),traceDataPre{2}(find([flightPaths.clusterIndex{:}]==flightTransitions.BtoA(flight_i))-length(traceDataPost{1}(:,1)),:,cell_i))
            hold on
            title(['Cell # ' num2str(cell_i)]);
            ylabel('Norm df/f');
            yt = get(gca,'YTick');
            %set(gca,'YTick',yt,'YTickLabel',yt/10);
            xt = get(gca, 'XTick');
            set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
            xlabel('time (s)');
            
            p17 = subplot(3,6,17); %plot neuron traces for each during flight B to A
            plot(1:length(traceDataFlight{2}(1,:,1)),traceDataFlight{2}(find([flightPaths.clusterIndex{:}]==flightTransitions.BtoA(flight_i))-length(traceDataPost{1}(:,1)),:,cell_i))
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
    
elseif plotPreFlightPostFlag == 1
    %% plot flight vs cell traces for each cluster for each cell with the pre/flight/post concatenated together in 1 plot
    
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
            normTraceRawAllAtoB(flight_i,:) = [flight_i,traceDataAll{1}(find([flightPaths.clusterIndex{:}]==flightTransitions.AtoB(flight_i)),:,cell_i)];
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
            normTraceRawAllBtoA(flight_i,:) = [traceDataAll{2}(find([flightPaths.clusterIndex{:}]==flightTransitions.BtoA(flight_i))-length(traceDataAll{1}(:,1)),:,cell_i)];
            plot(1:length(traceDataAll{2}(1,:,1)),traceDataAll{2}(find([flightPaths.clusterIndex{:}]==flightTransitions.BtoA(flight_i))-length(traceDataAll{1}(:,1)),:,cell_i)+flight_i*2)
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
    
    
    elseif plotMeanSpikeFlag == 1
    %% plot flight vs spike Means for each cluster for each cell with the pre/flight/post concatenated together in 1 plot
     
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
    
    %plot neural traces
    for cell_i = 1:length(snakeTrace.smoothTraceRawPreFlightPost{1}(:,1,:))
        %plot average spike activity for flights A to B
        p5 = subplot(5,2,[5 7]); %plot neuron traces for each pre-flight A to B
        plot(1:length(traceDataAll{1}(1,:)),traceDataAll{1}(cell_i,:));
        hold on
        plot(1:length(traceDataAll{1}(1,:)),smooth(traceDataAll{1}(cell_i,:),meanSmooth)); %plot the smoothed trace on top

        %normTraceRawAllAtoB = [flight_i,traceDataAll{1}(find([flightPaths.clusterIndex{:}]==flightTransitions.AtoB(flight_i)),:,cell_i)];
        %plot(1:length(traceDataAll{1}(1,:,1)),traceDataAll{1}(find([flightPaths.clusterIndex{:}]==flightTransitions.AtoB(flight_i)),:,cell_i)+flight_i*2)
        title(['Cell # ' num2str(cell_i)]);
        ylabel('# Spikes');
        yt = get(gca,'YTick');
        %set(gca,'YTick',yt,'YTickLabel',yt/10);
        set(gca,'xticklabel',{[]});
        xlim([1 length(traceDataAll{1}(1,:))]);
        %plot heatmap for each pre-flight A to B
        p7 = subplot(5,2,9);
        imagesc(traceDataAll{1}(cell_i,:))
        colormap('hot');
        xt = get(gca, 'XTick');
        set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
        xlabel('time (s)');
        ylabel('flight trials');
        
        %plot neuron traces for each pre-flight B to A
        p6 = subplot(5,2,[6 8]); 
        plot(1:length(traceDataAll{2}(1,:)),traceDataAll{2}(cell_i,:));
        hold on        
        plot(1:length(traceDataAll{2}(1,:)),smooth(traceDataAll{2}(cell_i,:),meanSmooth)); %plot the smoothed trace on top
        %normTraceRawAllBtoA(flight_i,:) = [traceDataAll{2}(find([flightPaths.clusterIndex{:}]==flightTransitions.BtoA(flight_i))-length(traceDataAll{1}(:,1)),:,cell_i)];
        %plot(1:length(traceDataAll{2}(1,:,1)),traceDataAll{2}(find([flightPaths.clusterIndex{:}]==flightTransitions.BtoA(flight_i))-length(traceDataAll{1}(:,1)),:,cell_i)+flight_i*2)
        title(['Cell # ' num2str(cell_i)]);
        ylabel('# Spikes');
        yt = get(gca,'YTick');
        %set(gca,'YTick',yt,'YTickLabel',yt/10);
        set(gca,'xticklabel',{[]});
        xlim([1 length(traceDataAll{2}(1,:))]);
         %plot heatmap for each pre-flight B to A       
        p8 = subplot(5,2,10);
        imagesc(traceDataAll{2}(cell_i,:))
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


