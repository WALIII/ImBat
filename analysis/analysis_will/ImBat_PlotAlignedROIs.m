function ImBat_PlotAlignedROIs(FlightAlignedROI,ROI_Data,flightPaths,varargin);
% ImBat_PlotAlignedROIs

% Plot individual ROIs aligned to flights from FlightAlignedROI, Rde lines
% indicate dates

% To generate aligned ROIs, first run:
%[FlightAlignedROI] = ImBat_Align_FC(CombinedROI,flightPaths,clst);


% WAL3
% 01/10/2021
nparams=length(varargin);
manIn = 0;

% User input
for i=1:2:nparams
	switch lower(varargin{i})
		case 'cells2use' % downsample video
             manIn =1;
			cells2use=varargin{i+1};
		case 'deinterlace' % de-interlace 60fps video
			deinterlace =varargin{i+1};
		case 'max_proj' % to make a max projectio
			max_proj=varargin{i+1};
           mkdir('processed/MAX');
    end
end



figure();
hold on;

clst = FlightAlignedROI.clust_number;

stop_time = 1;
sub2plot = 1;


% Plot all flights for 5 days, hilighting the clustered flights
A = flightPaths.tracjectoriesRaw*1000;

for i = 1:sub2plot%max(flightPaths.day); % number of days
    subplot(1,sub2plot,i);
    start_time = stop_time;
    stop_time = start_time+length(ROI_Data{i}.Alignment.out.Flights_Repaired(:,1))-1;
    % plot all flights for one day:
    bound = start_time:stop_time;
    
    hold on;
    % plot all flights
    plot1 =  plot3(A(1,bound),A(2,bound),A(3,bound),'k');
    %plot1.Color(4) = 0.4;
    % plot clust flights
    dayFlights = flightPaths.day(flightPaths.clusterIndex{clst});
    todaysFlight = find(dayFlights==i);
    for ii = 1:size(todaysFlight,1)
        plot3(squeeze(FlightAlignedROI.ClustFlight(1:FlightAlignedROI.ROI_OFF,1,todaysFlight(ii))),squeeze(FlightAlignedROI.ClustFlight(1:FlightAlignedROI.ROI_OFF,2,todaysFlight(ii))),squeeze(FlightAlignedROI.ClustFlight(1:FlightAlignedROI.ROI_OFF,3,todaysFlight(ii))),'LineWidth',2,'color','r')
    end
    grid on;
    view( -37.5000,30)
    axis off
    
end

% Stats on the clustered fights


% %% ROI Analysis
CutCells = FlightAlignedROI.C_raw;
ROI_ON = FlightAlignedROI.ROI_ON;


% get dates:
transition_points = find((diff(FlightAlignedROI.CutCells_date)>=1));  %?
transition_points = [1 transition_points size(CutCells,3)];


% get bound to plot
bound2plot = 1:500;
counter = 1;
col = hsv(size(transition_points,2)+1);
figure();
if manIn ==1;
else
    cells2use = 1:size(CutCells,1);
end
for i = 1:length(cells2use)
    i = cells2use(i);
    subplot(10,1,1:2);% plot shadding for each cell, over lay for each day
hold on;
    for ii = 1:size(transition_points,2)-1
        
        adata = zscore(squeeze(CutCells(i,bound2plot,transition_points(ii):transition_points(ii+1))),[],1)';
        
        L = size(adata,2);
        se = std(adata)/2;%/10;%sqrt(length(adata));
        mn = nanmean(adata);
        h = fill([1:L L:-1:1],[mn-se fliplr(mn+se)],col(ii,:)); alpha(0.5);
        plot(mn,'Color',col(ii,:));
    end
    hold off
    subplot(10,1,3:10);
    
    adata = zscore(squeeze(CutCells(i,bound2plot,:)),[],1)';
%     adata = adata-median(adata(:,500:600)')';
    imagesc(adata);
    % draw date lines:
    colormap(parula.^2)
    hold on;
    for ii = 1:size(transition_points,2)
        line([1,size(CutCells,2)], [transition_points(ii),transition_points(ii)], 'Color','r','LineWidth',2);
    end
    % Get axis handle
    ax = gca;
    % Set where ticks will be
    ax.XTick = [ROI_ON-60 ROI_ON-30 ROI_ON ROI_ON+30 ROI_ON+60  ROI_ON+90 ROI_ON+120 ROI_ON+150 ROI_ON+180];
    % Set TickLabels;
    ylabel('filghts')
    xlabel('time from takeoff');
    ax.XTickLabel = {'-2','-1','0','1','2','3','4','5','6'};
    title([ 'ROI ' num2str(i)]);
    pause();
    clf('reset');
    clear adata
end






