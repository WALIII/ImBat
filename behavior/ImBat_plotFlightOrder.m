function [plotFlightTimes] = ImBat_plotFlightOrder(flightPaths,varargin)

saveFlag = 0;
batName = [];
dateSesh = [];
sessionType = [];

nparams=length(varargin);
for i=1:2:nparams
    switch lower(varargin{i})
        case 'saveflag'
            saveFlag = varargin{i+1};
        case 'batname'
            batName = varargin{i+1};
        case 'datesesh'
            dateSesh = varargin{i+1};
        case 'sessiontype'
            sessionType = varargin{i+1};
    end
end

    jj = jet(flightPaths.nClusters); %make color vector
    plotFlightTimes = figure('units','normalized','outerposition',[0 0 1 0.5]);
    for clust_i = 1:flightPaths.nClusters
        %make height vector for the plot
        ht = ones(1,length(flightPaths.clusterIndex{clust_i}));
        %plot each flight in the color of the trajectory
        for flight_i = 1:length(flightPaths.clusterIndex{clust_i})
            legendData{clust_i} = plot(flightPaths.flight_starts_idx(flightPaths.clusterIndex{clust_i}(flight_i)),ht,'*','Color',jj(clust_i,:));
                ylim([0 1]);
            xline(flightPaths.flight_starts_idx(flightPaths.clusterIndex{clust_i}(flight_i)),'-','Color',jj(clust_i,:),'LineWidth',3);
            hold on;
        end
        %make dummy points for the legend so the colors match the trajectories
        l(clust_i) = plot([NaN,NaN],'*', 'color', jj(clust_i,:));
        
    end
    title(['Flight times of each trajectory: ' batName ' ' dateSesh ' ' sessionType]);
    xt = get(gca, 'XTick');
    set(gca,'XTick',xt,'XTickLabel',round(xt/120,1));
    xlabel('time (s)');
    set(gca,'yticklabel',{[]});
    legend([l(1) l(2),l(3)],{'traj 1','traj 2','traj 3'});

    
   if saveFlag ==1
        saveas(plotFlightTimes, [pwd '/' batName '_' dateSesh '_' sessionType '_plotFlightTimes.svg']);
        saveas(plotFlightTimes, [pwd '/' batName '_' dateSesh '_' sessionType '_plotFlightTimes.tif']);
        savefig(plotFlightTimes, [pwd '/' batName '_' dateSesh '_' sessionType '_plotFlightTimes.fig']);
    end
