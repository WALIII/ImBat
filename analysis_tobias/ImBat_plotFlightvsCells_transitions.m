function ImBat_plotFlightvsCells_transitions(snakeTrace,flightTransitions,flightPaths)
%plot the pre/post/and flight traces for each cell by their transitions

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

%for each flight in the subgroup AtoB
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
        plot(1:length(snakeTrace.normTraceRawPre{1}(1,:,1)),snakeTrace.normTraceRawPre{1}(find([flightPaths.clusterIndex{:}]==flightTransitions.AtoB(flight_i)),:,cell_i))
        hold on
        title(['Cell # ' num2str(cell_i)]);
        ylabel('Norm df/f');
        yt = get(gca,'YTick');
        %set(gca,'YTick',yt,'YTickLabel',yt/10);
        xt = get(gca, 'XTick');
        set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
        xlabel('time (s)');
        
        p14 = subplot(3,6,14); %plot neuron traces for each pre-flight A to B
        plot(1:length(snakeTrace.normTraceRawFlight{1}(1,:,1)),snakeTrace.normTraceRawFlight{1}(find([flightPaths.clusterIndex{:}]==flightTransitions.AtoB(flight_i)),:,cell_i))
        hold on
        title(['Cell # ' num2str(cell_i)]);
        ylabel('Norm df/f');
        yt = get(gca,'YTick');
        %set(gca,'YTick',yt,'YTickLabel',yt/10);
        xt = get(gca, 'XTick');
        set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
        xlabel('time (s)');
        
        p15 = subplot(3,6,15); %plot neuron traces for each pre-flight A to B
        plot(1:length(snakeTrace.normTraceRawPost{1}(1,:,1)),snakeTrace.normTraceRawPost{1}(find([flightPaths.clusterIndex{:}]==flightTransitions.AtoB(flight_i)),:,cell_i))
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
        plot(1:length(snakeTrace.normTraceRawPre{2}(1,:,1)),snakeTrace.normTraceRawPre{2}(find([flightPaths.clusterIndex{:}]==flightTransitions.BtoA(flight_i))-length(snakeTrace.normTraceRawPost{1}(:,1)),:,cell_i))
        hold on
        title(['Cell # ' num2str(cell_i)]);
        ylabel('Norm df/f');
        yt = get(gca,'YTick');
        %set(gca,'YTick',yt,'YTickLabel',yt/10);
        xt = get(gca, 'XTick');
        set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
        xlabel('time (s)');
        
        p17 = subplot(3,6,17); %plot neuron traces for each during flight B to A
        plot(1:length(snakeTrace.normTraceRawFlight{2}(1,:,1)),snakeTrace.normTraceRawFlight{2}(find([flightPaths.clusterIndex{:}]==flightTransitions.BtoA(flight_i))-length(snakeTrace.normTraceRawPost{1}(:,1)),:,cell_i))
        hold on
        title(['Cell # ' num2str(cell_i)]);
        ylabel('Norm df/f');
        yt = get(gca,'YTick');
        %set(gca,'YTick',yt,'YTickLabel',yt/10);
        xt = get(gca, 'XTick');
        set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
        xlabel('time (s)');
        
        p18 = subplot(3,6,18); %plot neuron traces for each post-flight B to A
        plot(1:length(snakeTrace.normTraceRawPost{2}(1,:,1)),snakeTrace.normTraceRawPost{2}(find([flightPaths.clusterIndex{:}]==flightTransitions.BtoA(flight_i))-length(snakeTrace.normTraceRawPost{1}(:,1)),:,cell_i))
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
