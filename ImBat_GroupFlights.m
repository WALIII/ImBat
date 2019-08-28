function out = ImBat_GroupFlights(ROI_Data);

clear allFlightTime allFlights
nClust = 10;


% plot all flights
col = jet(size(ROI_Data,2)); 
figure();
hold on;
for i = 1: size(ROI_Data,2);
    A = ROI_Data{1, i}.Alignment.out.flights;
    B = ROI_Data{1, i}.Alignment.out.Location_time;
    C = ones(size(B))*i;
    plot3(A(:,1),A(:,2),A(:,3),'Color',col(i,:));
    
    if i ==1;
        AllFlights = A;
        AllFlightsTime = B;
        DayIndex = C;
    else
        AllFlightsTemp = A;
        AllFlightsTimeTemp = B;
        DayIndexTemp = C;
        AllFlights = cat(1, AllFlights, AllFlightsTemp);
        AllFlightsTime = cat(1, AllFlightsTime, AllFlightsTimeTemp);
        DayIndex = cat(1,DayIndex,DayIndexTemp);

    end
    clear A AllFlightsTemp AllFlightsTimeTemp
end

grid on; 
colormap(jet)
colorbar;


% Segregate flights:

[out] =  ImBat_SegTrajectories(AllFlights,AllFlightsTime,'nclusters',nClust,'day_index',DayIndex);





% Plot all flight occurances:


for i = 1:size(ROI_Data,2)
    qq(i) = size(out.day(out.day ==i),1);
end
    
figure(); 
hold on;
for i = 1: 10;
    temp = zeros(1,size(ROI_Data,2));
 ax(i) =  subplot(5,2,i);
h = histogram(out.day(out.ClusterIndex{i}),'BinMethod','integers');
temp(round(h.BinEdges(2:end)-1)) = h.BinCounts;
histDat(:,i) = temp;
clear temp;
end
linkaxes(ax, 'xy');

figure(); bar(histDat(:,1:6),'stacked'); 
legend('flightpath 1','flightpath 2','flightpath 3','flightpath 4','flightpath 5','flightpath 6','flightpath 7')
title('Distribution of most common sterotyped flight paths');
xlabel('days')
ylabel('Number of flights');
