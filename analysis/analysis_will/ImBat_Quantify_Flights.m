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
plot(AllFlights(:,2:end),'color',[0 0 1 0.2]);
legend('unclustered','clustered Flights')
title('Total Flights ');
xlabel('day');
ylabel('Total # of Flights');


sumFlights = sum(AllFlights')';
figure(); 
hold on;
plot((AllFlights(:,1)./sumFlights)*100,'--r','LineWidth',2);
plot((AllFlights(:,2:end)./sumFlights)*100,'color',[0 0 1 0.2]);
legend('unclustered','clustered Flights')
title('Flight Ratio');
xlabel('day');
ylabel('% of Flights');


% Percent of clustered flights
sumClustFlights = sum(AllFlights(:,2:end)')';
figure(); 
hold on;
col = [0 1 0];
for ii = 2:size(AllFlights,2)
    if ii <5 % top 3 flights
plot((AllFlights(:,ii)./sumClustFlights)*100,'color',col,'Linewidth',2);
col = col+ [0.4 -0.4  00];
    else
plot((AllFlights(:,ii)./sumClustFlights)*100,'color',[0 0 1 0.2]);
    end
end

legend('Top3 Flights','','','Other Flights')
title('Ratio of clustered Flights ');
xlabel('day');
ylabel('% of Flights');


figure(); 
hold on;
plot((AllFlights(:,1)./sumFlights)*100,'--r','LineWidth',2);
plot((sum(AllFlights(:,2:end)')./sumFlights')*100,'g','LineWidth',2);
legend('unclustered','clustered Flights')
title('Clustered vs Unclustered Flights ');
xlabel('day');
ylabel('% of Flights');



figure(); 
hold on;
plot(sumFlights,'color',[0 0 1 0.2],'LineWidth',2)
plot( AllFlights(:,1),'--r','LineWidth',2);
plot(sumClustFlights,'g','LineWidth',2);
legend('all Flights', 'unique Flights','clustered Flights')
title('Number of Flights per group ');
xlabel('day');
ylabel('# of Flights');
