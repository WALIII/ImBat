function [dark_cluster, out] = ImBat_AfterDark(flightPaths)

% Find flights that happened in the dark:


% Get flightpaths of type:

FF = flightPaths.flight_starts_idx;


%figure(); plot(flightPaths.AllFlightsTime)

% plot all flights
figure();
histogram(flightPaths.AllFlightsTime(FF)/60,100);

% plot unclustered
for i = 1:4
    cluster = i;
    flightPaths.clusterIndex{cluster};
    FF2 = flightPaths.flight_starts_idx(flightPaths.clusterIndex{cluster});
    
    figure();
    hold on;
    histout{i} = histogram(flightPaths.AllFlightsTime(FF2)/60,60);
    xlabel('time (min)');
    ylabel('count');
    title(['Cluster ',num2str(i)]);
    plot([20 20],[0 10],'--r')
    plot([40 40],[0 10],'--r')
    clear FF2
    xlim([0 60]);
end




% Export Data
for i = 1:4; % for cluster
    FF2 = flightPaths.flight_starts_idx(flightPaths.clusterIndex{i});
    dc = flightPaths.AllFlightsTime(FF2)/60;
    out.Light_cluster{i}.FirstLight = flightPaths.clusterIndex{i}(find(dc<20));
    out.Light_cluster{i}.Darkness = flightPaths.clusterIndex{i}(find(dc>20 & dc<40));
    out.Light_cluster{i}.LastLight = flightPaths.clusterIndex{i}(find(dc>40));
    
    dark_cluster{i} = find(dc>20 & dc<40); % would use out.Light_cluster{i}.Darkness' to index into the flight order... 
    clear FF2 dc;
end
FF3 = flightPaths.flight_starts_idx;
dc = flightPaths.AllFlightsTime(FF3)/60;

% for all flights
out.Light_all.FirstLight = find(dc<20);
out.Light_all.Darkness = find(dc>20 & dc<40);
out.Light_all.LastLight = find(dc>40);
out.HistData = histout;




% some basic plotting
%% Plot unclustered Ratio
figure();
hold on;
%plot(medfilt1(histout{1}.Values./(histout{1}.Values+histout{2}.Values+histout{3}.Values),1,[],2))
plot(smooth(histout{1}.Values./(histout{1}.Values+histout{2}.Values+histout{3}.Values+histout{4}.Values),7),'--r')

xlabel('time (min)');
ylabel('proportion of unique flights');
title(['proportion of unique flights to top flights']);
plot([20 20],[0 1],'-b')
plot([40 40],[0 1],'-b')

%% Plot Flights

figure();
hold on;
A = flightPaths.tracjectoriesRaw*1000;

% Pre-Light
subplot(1,3,1); hold on;
Ind2use = out.Light_all.FirstLight;
for iii = 1:length(Ind2use)  
    bound = flightPaths.flight_starts_idx(Ind2use(iii)):flightPaths.flight_ends_idx(Ind2use(iii));
    plot1 =  plot3(A(1,bound),A(2,bound),A(3,bound),'color',[0 0 0 0.5]); % plot all flights
end
title('Lights ON');
% Darkness
subplot(1,3,2); hold on;
Ind2use = out.Light_all.Darkness;
for iii = 1:length(Ind2use)
    bound = flightPaths.flight_starts_idx(Ind2use(iii)):flightPaths.flight_ends_idx(Ind2use(iii));
    plot1 =  plot3(A(1,bound),A(2,bound),A(3,bound),'color',[0 0 0 0.5]); % plot all flights
end
title('Lights Off');

% Last Light
subplot(1,3,3); hold on;
Ind2use = out.Light_all.LastLight;
for iii = 1:length(Ind2use)
    bound = flightPaths.flight_starts_idx(Ind2use(iii)):flightPaths.flight_ends_idx(Ind2use(iii));
    plot1 =  plot3(A(1,bound),A(2,bound),A(3,bound),'color',[0 0 0 0.5]); % plot all flights
end
title('Lights Return on');

%% Now, Plot the same figure, but overlay the clustered flights:

%% %% %% Plot Flights

figure();
cluster = 3;
col = [0 0 1 1];
hold on;
A = flightPaths.tracjectoriesRaw*1000;

% Pre-Light
subplot(1,3,1); hold on;
Ind2use = out.Light_all.FirstLight;
Ind2use2 = out.Light_cluster{cluster}.FirstLight  ;

for iii = 1:length(Ind2use)  
    bound = flightPaths.flight_starts_idx(Ind2use(iii)):flightPaths.flight_ends_idx(Ind2use(iii));
    plot1 =  plot3(A(1,bound),A(2,bound),A(3,bound),'color',[0 0 0 0.5]); % plot all flights
end
for ii = 1: length(Ind2use2);
    bound2 = flightPaths.flight_starts_idx(Ind2use2(ii)):flightPaths.flight_ends_idx(Ind2use2(ii));
    plot2 =  plot3(A(1,bound2),A(2,bound2),A(3,bound2),'color',col,'LineWidth',2); % plot all flights
end
title('Lights ON');

% Darkness
subplot(1,3,2); hold on;
Ind2use = out.Light_all.Darkness;
Ind2use2 = out.Light_cluster{cluster}.Darkness    ;
for iii = 1:length(Ind2use)  
    bound = flightPaths.flight_starts_idx(Ind2use(iii)):flightPaths.flight_ends_idx(Ind2use(iii));
    plot1 =  plot3(A(1,bound),A(2,bound),A(3,bound),'color',[0 0 0 0.5]); % plot all flights
end
for ii = 1: length(Ind2use2);
    bound2 = flightPaths.flight_starts_idx(Ind2use2(ii)):flightPaths.flight_ends_idx(Ind2use2(ii));
    plot2 =  plot3(A(1,bound2),A(2,bound2),A(3,bound2),'color',col,'LineWidth',2); % plot all flights
end
title('Lights Off');

% Last Light
subplot(1,3,3); hold on;
Ind2use = out.Light_all.LastLight;
Ind2use2 = out.Light_cluster{cluster}.LastLight    ;
for iii = 1:length(Ind2use)  
    bound = flightPaths.flight_starts_idx(Ind2use(iii)):flightPaths.flight_ends_idx(Ind2use(iii));
    plot1 =  plot3(A(1,bound),A(2,bound),A(3,bound),'color',[0 0 0 0.5]); % plot all flights
end
for ii = 1: length(Ind2use2);
    bound2 = flightPaths.flight_starts_idx(Ind2use2(ii)):flightPaths.flight_ends_idx(Ind2use2(ii));
    plot2 =  plot3(A(1,bound2),A(2,bound2),A(3,bound2),'color',col,'LineWidth',2); % plot all flights
end
title('Lights Return on');




