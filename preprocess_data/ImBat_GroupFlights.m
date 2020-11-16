function [flightPaths] = ImBat_GroupFlights(ROI_Data,varargin);
% Group flights across days, now w/ Angelo's function
% updated 10/20/2020

% WAL3

% Default params
n_splines = 20;
dist_met = 1.2; %1.2 % lower is more selective ( more clusters)

pause(2);
do_mtf =0;
% Manual inputs
% vin=varargin;
% for i=1:length(vin)
%     if isequal(vin{i},'mtf') % manually inputing a sort order
%         MTF=vin{i+1};
%         do_mtf = 1;
%     end
%
% end

% User inputs overrides
nparams=length(varargin);
for i=1:2:nparams
    switch lower(varargin{i})
        case 'mtf'
            MTF=varargin{i+1};
            do_mtf = 1;
        case 'dist'
            dist_met=varargin{i+1};
            disp(['WARNING: distance metric set to: ', num2str(dist_met), ' default is 1.2']);
        case 'n_splines'
            n_splines=varargin{i+1};
            disp(['WARNING: splines metric set to: ', num2str(dist_met), ' default is 20']);
    end
end


% Align to track file
if do_mtf ==1;
    disp('aligning to master track file');
    for i = 1: size(ROI_Data,2);
        d = ROI_Data{i}.date(end-5:end);
        try
            ROI_Data{1, i}.Alignment.out.flights = ImBat_Align_Tracking(ROI_Data{1, i}.Alignment.out.Flights_Repaired,MTF,d);
        catch
            disp([' no track data for ',d]);
            ROI_Data{1, i}.Alignment.out.flights= ROI_Data{1, i}.Alignment.out.Flights_Repaired;
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
    
    D = ROI_Data{1, i}.Alignment.out.video_times(1:end-1); % align timestamps
    % trim the end, otherwise the flights will be longer or shorter than the
    % calcium..
    if max(D)>max(B); disp('adding extra timepoint to flight data');
        A = cat(1,A,A(end,:));
        B = cat(1,B,max(D));
    else
        disp(' Calcium is shorter than flights, trimming data'); % cut flight data down
        % find closest
        [minValue_1,closestIndex_1] = min(abs(D(end)-B));
        % cut off trailing data
        A(closestIndex_1-1:end,:) = [];
        B(closestIndex_1-1:end,:) = [];
        A = cat(1,A,A(end,:));
        B = cat(1,B,max(D));
        
        % add closest timepoint
    end
    
    if i ==1;
        AllFlights = A;
        AllFlightsTime = B;
        AllFlightsMasterTime = B; % accumulating time as an index
        DayIndex = C;
    else
        AllFlightsTemp = A;
        AllFlightsTimeTemp = B;
        DayIndexTemp = C;
        AllFlights = cat(1, AllFlights, AllFlightsTemp);
        AllFlightsTime = cat(1, AllFlightsTime, AllFlightsTimeTemp);
        AllFlightsMasterTime = cat(1, AllFlightsMasterTime, AllFlightsTimeTemp+max(AllFlightsMasterTime));
        
        DayIndex = cat(1,DayIndex,DayIndexTemp);
        
    end
    clear A AllFlightsTemp AllFlightsTimeTemp
end

grid on;
colormap(jet)
colorbar;


% Segregate flights:

%[out] =  ImBat_SegTrajectories(AllFlights,AllFlightsTime,'nclusters',8,'day_index',DayIndex);
Fs = ROI_Data{1, 1}.ROIs.results.metadata.cnmfe.Fs;
[flightPaths] = ImBat_flightsAngelo(AllFlights,AllFlightsTime,'fs',Fs,'n_splines',n_splines,'dist',dist_met,'day_index',DayIndex);

flightPaths.AllFlights = AllFlights;
flightPaths.AllFlightsTime = AllFlightsTime;
flightPaths.AllFlightsMasterTime = AllFlightsMasterTime;

% Plot all flight occurances:


% for i = 1:size(ROI_Data,2)
%     qq(i) = size(out.day(out.day ==i),1);
% end

flight_clust_size = size(flightPaths.clusterIndex,2);
if flight_clust_size>9;
    num2plot = 10;
else
    num2plot = flight_clust_size;
end
figure();
hold on;
for i = 1:num2plot
    temp = zeros(1,size(ROI_Data,2));
    ax(i) =  subplot(5,2,i);
    h = histogram(flightPaths.day(flightPaths.clusterIndex{i}),'BinMethod','integers');
    temp(round(h.BinEdges(2:end)-1)) = h.BinCounts;
    
    h2 = histogram(flightPaths.day,'BinMethod','integers');
    temp2(round(h2.BinEdges(2:end)-1)) = h2.BinCounts;
    histDat(:,i) = temp; % flights of this type
    histDat2(:,1) = temp2; % all flights
    clear temp temp2;
end
linkaxes(ax, 'xy');

hist2plot = num2plot-1;% top flights
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
