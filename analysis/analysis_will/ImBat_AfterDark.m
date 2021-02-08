function [dark_cluster, out] = ImBat_AfterDark(flightPaths,varargin)
% Find flights that happened in the dark:

% Assumption is the lights are turned off exactly 20min, then back on at
% 40min

% Get flightpaths of type:
FF = flightPaths.flight_starts_idx;
FlightPaths = 2;

% User Inputs:
nparams=length(varargin);
if mod(nparams,2)>0
    error('Parameters must be specified as parameter/value pairs');
end

for i=1:2:nparams
    switch lower(varargin{i})
        case 'flightpaths'
            FlightPaths=varargin{i+1};
    end
end
%figure(); plot(flightPaths.AllFlightsTime)

% plot all flights
figure();
histogram(flightPaths.AllFlightsTime(FF)/60,100);

% plot unclustered
for i = 1:size(FlightPaths,2)+1
    cluster = i;
    flightPaths.clusterIndex{cluster};
    FF2 = flightPaths.flight_starts_idx(flightPaths.clusterIndex{cluster});
    
    figure();
    hold on;
    histout{i} = histogram(flightPaths.AllFlightsTime(FF2)/60,'NumBins',60,'BinLimits',[1,60]);
    xlabel('time (min)');
    ylabel('count');
    title(['Cluster ',num2str(i)]);
    plot([20 20],[0 10],'--r')
    plot([40 40],[0 10],'--r')
    clear FF2
    xlim([0 60]);
end




% Export Data
for i = 1:size(FlightPaths,2)+1; % for cluster
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
    plot1 =  plot3(A(1,bound),A(2,bound),A(3,bound),'color',[0 0 0 0.1]); % plot all flights
end
axis off
title('Lights ON');
% Darkness
subplot(1,3,2); hold on;
Ind2use = out.Light_all.Darkness;
for iii = 1:length(Ind2use)
    bound = flightPaths.flight_starts_idx(Ind2use(iii)):flightPaths.flight_ends_idx(Ind2use(iii));
    plot1 =  plot3(A(1,bound),A(2,bound),A(3,bound),'color',[0 0 0 0.1]); % plot all flights
end
title('Lights Off');
axis off
% Last Light
subplot(1,3,3); hold on;
Ind2use = out.Light_all.LastLight;
for iii = 1:length(Ind2use)
    bound = flightPaths.flight_starts_idx(Ind2use(iii)):flightPaths.flight_ends_idx(Ind2use(iii));
    plot1 =  plot3(A(1,bound),A(2,bound),A(3,bound),'color',[0 0 0 0.1]); % plot all flights
end
title('Lights Return on');
axis off
%% Now, Plot the same figure, but overlay the clustered flights:

%% %% %% Plot Flights


figure();
col = {[1 0 0 0.2],[0 0 1 0.2],[0 1 0 0.2],[1 0 1 0.2],[1 1 0 0.2],[0 1 1 0.2]};

for i = 1:length(FlightPaths);
cluster = FlightPaths(i);
hold on;
A = flightPaths.tracjectoriesRaw*1000;

% Pre-Light
subplot(1,3,1); hold on;
Ind2use = out.Light_all.FirstLight;
Ind2use2 = out.Light_cluster{cluster}.FirstLight  ;
axis off
if i ==1;
for iii = 1:length(Ind2use)  
    bound = flightPaths.flight_starts_idx(Ind2use(iii)):flightPaths.flight_ends_idx(Ind2use(iii));
    plot1 =  plot3(A(1,bound),A(2,bound),A(3,bound),'color',[0 0 0 0.1]); % plot all flights

end
end
for ii = 1: length(Ind2use2);
    bound2 = flightPaths.flight_starts_idx(Ind2use2(ii)):flightPaths.flight_ends_idx(Ind2use2(ii));
    plot2 =  plot3(A(1,bound2),A(2,bound2),A(3,bound2),'color',col{i},'LineWidth',2); % plot all flights
end
title('Lights ON');
axis off

% Darkness
subplot(1,3,2); hold on;
Ind2use = out.Light_all.Darkness;
Ind2use2 = out.Light_cluster{cluster}.Darkness;
if i ==1;
for iii = 1:length(Ind2use)  
    bound = flightPaths.flight_starts_idx(Ind2use(iii)):flightPaths.flight_ends_idx(Ind2use(iii));
    plot1 =  plot3(A(1,bound),A(2,bound),A(3,bound),'color',[0 0 0 0.1]); % plot all flights
end
end
for ii = 1: length(Ind2use2);
    bound2 = flightPaths.flight_starts_idx(Ind2use2(ii)):flightPaths.flight_ends_idx(Ind2use2(ii));
    plot2 =  plot3(A(1,bound2),A(2,bound2),A(3,bound2),'color',col{i},'LineWidth',2); % plot all flights
end
title('Lights Off');
axis off

% Last Light
subplot(1,3,3); hold on;
Ind2use = out.Light_all.LastLight;
Ind2use2 = out.Light_cluster{cluster}.LastLight;
if i ==1;
for iii = 1:length(Ind2use)  
    bound = flightPaths.flight_starts_idx(Ind2use(iii)):flightPaths.flight_ends_idx(Ind2use(iii));
    plot1 =  plot3(A(1,bound),A(2,bound),A(3,bound),'color',[0 0 0 0.1]); % plot all flights
end
end
for ii = 1: length(Ind2use2);
    bound2 = flightPaths.flight_starts_idx(Ind2use2(ii)):flightPaths.flight_ends_idx(Ind2use2(ii));
    plot2 =  plot3(A(1,bound2),A(2,bound2),A(3,bound2),'color',col{i},'LineWidth',2); % plot all flights
end
title('Lights Return on');
end
axis off

% reate custom color:
col2use(1,:) = [0 0 0];
col2use(2,:) = [1 0 0];
col2use(3,:) = [0 0 1];
col2use(4,:) = [1 0 1];
col2use(5,:) = [0 1 0];
col2use(6,:) = [0 1 1];
col2use(7,:) = [1 1 0];



clear a
for i = 1:size(FlightPaths,2)+1;
a(:,i) = histout{i}.Values/max(flightPaths.day);
end
figure();
hold on;
b = bar(a,'stacked');
    plot([20 20],[0 30/max(flightPaths.day)],'--k','LineWidth',2)
    plot([40 40],[0 30/max(flightPaths.day)],'--k','LineWidth',2)
title(' Number of flights vs time of day')
xlabel('time in session ( min)');
ylabel('avg flights/min');
legend('Unique','FlightPath 1','FlightPath 2','FlightPath 3','FlightPath 4');


for K = 1 : length(b);
    if K ==1;
       b(K).FaceAlpha = 0.5;
       b(K).EdgeAlpha = 0.5;
    else
       b(K).FaceAlpha = 0.8;
       b(K).EdgeAlpha = 0.8;
    b(K).FaceColor = col2use(K,:).'; 
    b(K).EdgeColor = col2use(K,:).'; 
    end
end


