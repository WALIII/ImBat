function ImBat_Quantify_Flights(flightPaths);

% TO DO:

% 1. ratio of clustered flights to all flights
% 2. Ratio of top clusters to all flights
% 3. Ratio of top clusters to all clustered flights

clst2use  = size(flightPaths.clusterIndex,2);

for clst = 1:clst2use;
    % get all flights of a certian cluster:
dayFlights = flightPaths.day(flightPaths.clusterIndex{clst});
% store as a function of days 
for day = 1:max(flightPaths.day);
AllFlights(day,clst) = size(find(dayFlights==day),1);
end
end

figure(); 
hold on;
plot(AllFlights(:,1),'--r','LineWidth',2);
plot(AllFlights(:,2:end),'b');
legend('unclustered','clustered Flights')
title('Total Flights ');
xlabel('day');
ylabel('Total # of Flights');


sumFlights = sum(AllFlights')';
figure(); 
hold on;
plot((AllFlights(:,1)./sumFlights)*100,'--r','LineWidth',2);
plot((AllFlights(:,2:end)./sumFlights)*100,'b');
legend('unclustered','clustered Flights')
title('Flight Ratio');
xlabel('day');
ylabel('% of Flights');


% Percent of clustered flights
sumClustFlights = sum(AllFlights(:,2:end)')';
figure(); 
hold on;
for ii = 2:size(AllFlights,2)
    if ii <5 % top 3 flights
plot((AllFlights(:,ii)./sumClustFlights)*100,'g','Linewidth',2);
    else
plot((AllFlights(:,ii)./sumClustFlights)*100,'b');
    end
end

legend('Top3 Flights','','','Other FLights')
title('Ratio of clustered Flights ');
xlabel('day');
ylabel('% of Flights');


figure(); 
hold on;
plot(AllFlights(:,1)./sumFlights,'--r','LineWidth',2);
plot(sum(AllFlights(:,2:end)')./sumFlights','g','LineWidth',2);
legend('unclustered','clustered Flights')
title('Clustered vs Unclustered Flights ');
xlabel('day');
ylabel('% of Flights');
