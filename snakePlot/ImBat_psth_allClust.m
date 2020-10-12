function [psth_plots] = ImBat_psth_allClust(flightPaths,snakeTrace_cRaw,varargin)
% User inputs overrides
saveFlag = 0;

nparams=length(varargin);
for i=1:2:nparams
    switch lower(varargin{i})
        case 'saveflag'
            saveFlag = varargin{i+1};
        case 'analysisfolder'
            analysis_Folder = varargin{i+1};
        case 'batname'
            batName=varargin{i+1};
        case 'datesesh'
            dateSesh = varargin{i+1};
        case 'sessiontype'
            sessionType = varargin{i+1};
    end
end

%% plot
% for each cell, for each cluster, plot flight paths, velocity, mean/std, and each trial
nCells = size(snakeTrace_cRaw.smoothTraceRawFlight{1},3); %number of cells in this day
nClusts = size(snakeTrace_cRaw.smoothTraceRawFlight,2); %number of clusters in this day

for cell_i = 1:nCells
    psth_each_cell = figure('units','normalized','outerposition',[0 0.1 0.8 0.8]);
    sgtitle(['Flight Aligned Activity per cluster: ' batName ' ' dateSesh ' Cell# ' num2str(cell_i)]);
    colClust = jet(nClusts); %color each cluster differently)
    for clust_i = 1:nClusts
        %plot flight paths of each cluster
        nFlights = size(snakeTrace_cRaw.smoothTraceRawFlight{clust_i},1); %number of flights in this cluster
        colFlight = jet(nFlights); %color each cluster differently)
        for flight_i = 1:nFlights
            plotPaths = subplot(10,nClusts,[clust_i nClusts+clust_i]);
            plot3(flightPaths.pos(1,:,flightPaths.clusterIndex{clust_i}(flight_i)),flightPaths.pos(2,:,flightPaths.clusterIndex{clust_i}(flight_i)),flightPaths.pos(3,:,flightPaths.clusterIndex{clust_i}(flight_i)),...
                '-','LineWidth',1,'Color', colClust(clust_i,:));
            hold on;
            view(0,90);
            xlim([-3 3])
            ylim([-3 3])
            title(['Clust ' num2str(clust_i)]);
            if clust_i == 1
                xlabel('m'); ylabel('m');
            else
                set(gca,'xticklabel',[],'yticklabel',[]);
            end
            
            %plot velocity of each cluster
            plotVelocity = subplot(10,nClusts,[clust_i+2*nClusts clust_i+3*nClusts]);
            plot(1:length(snakeTrace_cRaw.smoothSpeedRawFlight{clust_i}(flight_i,:)),snakeTrace_cRaw.smoothSpeedRawFlight{clust_i}(flight_i,:));
            hold on;
            if clust_i == 1
                ylabel('Velocity (m/s)');
                yt = get(gca,'YTick');
                xt = get(gca,'XTick');
                set(gca,'xticklabel',[]);
                %set(gca,'YTick',yt,'YTickLabel',yt,'xTickLabel',round(xt/120,1));
            else
                set(gca,'xticklabel',[],'yticklabel',[]);
            end
            ylim([0 4.5]);
            xlim([0 length(snakeTrace_cRaw.smoothSpeedRawFlight{clust_i}(flight_i,:))]);
            
            %plot each flight activity
            plotIndivTraces = subplot(10,nClusts,[clust_i+6*nClusts:nClusts:clust_i+9*nClusts]);
            plot(1:length(snakeTrace_cRaw.normTraceRawFlight{clust_i}(flight_i,:,cell_i)),(snakeTrace_cRaw.normTraceRawFlight{clust_i}(flight_i,:,cell_i))+flight_i*6,'Color',colClust(clust_i,:));
            hold on;
            if clust_i == 1
                ylabel('Flight trials (df/f)');
            end
            xlabel('Time (s)');
            xt = get(gca, 'XTick');
            set(gca,'XTick',xt,'XTickLabel',round(xt/30,1));
            xlim([0 length(snakeTrace_cRaw.normTraceRawFlight{clust_i}(flight_i,:,cell_i))]);
        end
        %plot mean and stdev of activity per cluster
        plotMeanTraces = subplot(10,nClusts,[clust_i+4*nClusts clust_i+5*nClusts]);
        boundedline(1:length(snakeTrace_cRaw.normTraceFlight{clust_i}(cell_i,:)),snakeTrace_cRaw.normTraceFlight{clust_i}(cell_i,:),snakeTrace_cRaw.sdTraceFlight{clust_i}(cell_i,:),'alpha','cmap',colClust(clust_i,:));
        hold on;
        ylabel('Mean stdev df/f');
        set(gca,'xticklabel',[]);
        xlim([0 length(snakeTrace_cRaw.normTraceFlight{clust_i}(cell_i,:))]);
    end
    
    %% save
    % Save 'place cells' as jpg and fig files..
    saveas(psth_each_cell,[pwd filesep batName '_' dateSesh '_' sessionType '_psth_clusts' num2str(cell_i) '.tif']);
    savefig(psth_each_cell,[pwd filesep batName '_' dateSesh '_' sessionType '_psth_clusts' num2str(cell_i) '.fig']);
    %saveas(psth_each_cell,[pwd filesep batName '_' dateSesh '_' sessionType '_psth_clusts' num2str(cell_i) '.svg']);
    close gcf;
end