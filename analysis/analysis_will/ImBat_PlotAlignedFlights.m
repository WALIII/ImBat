function ImBat_PlotAlignedFlights(ROI_Data,flightPaths,FlightAlignedROI,days2use)

figure(); 
stop_time = 1;
A = flightPaths.tracjectoriesRaw*1000;

clst = FlightAlignedROI.clust_number;
counter = 1;
for i = 1:length(days2use);%max(flightPaths.day); % number of days
    
    day2use = days2use(i);
    subplot(1,length(days2use),counter);
    start_time = stop_time;
    stop_time = start_time+length(ROI_Data{day2use}.Alignment.out.Flights_Repaired(:,1))-1;
    % plot all flights for one day:
    bound = start_time:stop_time;
    
    hold on;
    % plot all flights
    plot1 =  plot3(A(1,bound),A(2,bound),A(3,bound),'color',[0 0 0 0.5]);
    %plot1.Color(4) = 0.4;
    % plot clust flights
    dayFlights = flightPaths.day(flightPaths.clusterIndex{clst});
    todaysFlight = find(dayFlights==day2use);
    for ii = 1:size(todaysFlight,1)
        plot3(squeeze(FlightAlignedROI.ClustFlight(1:FlightAlignedROI.ROI_OFF,1,todaysFlight(ii))),squeeze(FlightAlignedROI.ClustFlight(1:FlightAlignedROI.ROI_OFF,2,todaysFlight(ii))),squeeze(FlightAlignedROI.ClustFlight(1:FlightAlignedROI.ROI_OFF,3,todaysFlight(ii))),'LineWidth',2,'color','r')
    end
    grid on;
    view( -37.5000,30)
   xlim([-2500 2500]);
    ylim([-2500 2500]);
     zlim([0 2500]);
    axis off
    counter = counter+1;
      title(ROI_Data{day2use}.date);

end


