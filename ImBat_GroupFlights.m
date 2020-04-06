function out = ImBat_GroupFlights(ROI_Data,varargin);

clear allFlightTime allFlights

do_mtf =0;
% Manual inputs
vin=varargin;
for i=1:length(vin)
    if isequal(vin{i},'mtf') % manually inputing a sort order
        MTF=vin{i+1};
        do_mtf = 1;
    end
end



% Align to track file
if do_mtf ==1;
    disp('aligning to master track file');
    for i = 1: size(ROI_Data,2);
        d = ROI_Data{i}.date(end-5:end);
        try
            ROI_Data{1, i}.Alignment.out.flights = ImBat_Align_Tracking(ROI_Data{1, i}.Alignment.out.flights,MTF,d);
        catch
            disp([' no track data for ',d]);
        end
    end
    disp('finished alignment');
end

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

[out] =  ImBat_SegTrajectories(AllFlights,AllFlightsTime,'nclusters',10,'day_index',DayIndex);





% Plot all flight occurances:


for i = 1:size(ROI_Data,2)
    qq(i) = size(out.day(out.day ==i),1);
end

figure();
hold on;
for i = 1:10
    temp = zeros(1,size(ROI_Data,2));
    ax(i) =  subplot(5,2,i);
    h = histogram(out.day(out.ClusterIndex{i}),'BinMethod','integers');
    temp(round(h.BinEdges(2:end)-1)) = h.BinCounts;
    
    h2 = histogram(out.day,'BinMethod','integers');
    temp2(round(h2.BinEdges(2:end)-1)) = h2.BinCounts;
    histDat(:,i) = temp; % flights of this type
    histDat2(:,1) = temp2; % all flights
    clear temp temp2;
end
linkaxes(ax, 'xy');

hist2plot = 9;% top flights
colors = hsv(hist2plot+4);
remains = histDat2 - sum(histDat(:,1:hist2plot)')'; % get remainder
figure();
b = bar(cat(2,histDat(:,1:hist2plot),remains),'stacked');
colormap(lines(winter));
legend('flightpath 1','flightpath 2','flightpath 3','flightpath 4','flightpath 5','flightpath 6','flightpath 7','flightpath 8','flightpath 9','All other flights')
title('Distribution of most common sterotyped flight paths');
xlabel('days')
ylabel('Number of flights');

for K = 1 : length(b); b(K).FaceColor = colors(K,:).'; end
