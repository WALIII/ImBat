function output = ImBat_Quantify_Flights(flightPaths);

% TO DO:

% 1. ratio of clustered flights to all flights
% 2. Ratio of top clusters to all flights
% 3. Ratio of top clusters to all clustered flights

DateString = flightPaths.Dates{1, 1};
formatIn = 'yymmdd';

compare_date = datenum(DateString,formatIn);

for i = 1:size(flightPaths.Dates,2)
    DateString = flightPaths.Dates{i};
    temp_date = datenum(DateString,formatIn);
    diff_date(i) = temp_date-compare_date;
end
diff_date = diff_date+1;
% Pre-allocate

clst2use  = size(flightPaths.clusterIndex,2);

for clst = 1:clst2use;
    % get all flights of a certian cluster:
dayFlights = flightPaths.day(flightPaths.clusterIndex{clst});
% store as a function of days 
for day = 1:max(flightPaths.day);
AllFlights(diff_date(day),clst) = size(find(dayFlights==day),1);
AllFlights_List(day,clst) = size(find(dayFlights==day),1); % for exporting 

end
end

figure(); 
hold on;
plot(AllFlights(:,1),'-r','LineWidth',2);
plot(AllFlights(:,2:end),'color',[0 0 1 0.2]);
legend('unclustered','clustered Flights')
title('Total Flights ');
xlabel('day');
ylabel('Total # of Flights');


sumFlights = sum(AllFlights')';
figure(); 
hold on;
plot((AllFlights(:,1)./sumFlights)*100,'-r','LineWidth',2);
plot((AllFlights(:,2:end)./sumFlights)*100,'color',[0 0 1 0.2]);
legend('unclustered','clustered Flights')
title('Flight Ratio');
xlabel('day');
ylabel('% of Flights');


% Percent of clustered flights
sumClustFlights = sum(AllFlights(:,2:end)')';
figure(); 
subplot(1,2,1)
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

subplot(1,2,2)
% Percent of clustered flights
hold on;
col = [0 1 0];
for ii = 2:size(AllFlights,2)
    if ii <5 % top 3 flights
plot(AllFlights(:,ii),'color',col,'Linewidth',2);
col = col+ [0.4 -0.4  00];
    else
plot((AllFlights(:,ii)),'color',[0 0 1 0.2]);
    end
end

legend('Top3 Flights','','','Other Flights')
title('Number of clustered Flights ');
xlabel('day');
ylabel('# Flights');







figure(); 
hold on;
plot((AllFlights(:,1)./sumFlights)*100,'r','LineWidth',2);
plot((sum(AllFlights(:,2:end)')./sumFlights')*100,'g','LineWidth',2);
legend('unclustered','clustered Flights')
title('Clustered vs Unclustered Flights ');
xlabel('day');
ylabel('% of Flights');



figure(); 
hold on;
plot(sumFlights,'color',[0 0 1 0.2],'LineWidth',2)
plot(AllFlights(:,1),'-r','LineWidth',2);
plot(sumClustFlights,'g','LineWidth',2);
legend('all Flights', 'unique Flights','clustered Flights')
title('Number of Flights per group ');
xlabel('day');
ylabel('# of Flights');

for i = 1: 4;
days2use{i} = find(AllFlights_List(:,i)>1);
end

output.days2use = days2use; % for max projection creation
[FL_ID FirstThreeClusters]  = ImBat_FlightNum_OverDays(flightPaths);
output.FL_ID = FL_ID;
output.FirstThreeClusters = FirstThreeClusters;


function [FL_ID FirstThreeClusters] = ImBat_FlightNum_OverDays(flightPaths);
%  ImBat_FlightNum_OverDays.m

% Time of day:


% ratio of flights
for i = 1:max(flightPaths.day); % days
    todaysFlightsID = find(flightPaths.day==i);
    
todaysFlights = flightPaths.id(todaysFlightsID);

for ii = 1: max(todaysFlights);
    FL_ID(i,ii) = size(find(todaysFlights == ii),1);
end
end

figure(); 
hold on;
plot(sum(FL_ID,2));
title('Num of flight');
ylabel('# of Flights');
xlabel('Sessions')


UniqueFlights = FL_ID(:,1);
ClusterdFlights = sum(FL_ID(:,2:end),2);
FirstThreeClusters = FL_ID(:,1:4)./sum(FL_ID,2);

figure();
plotSpread(FirstThreeClusters);
title('Ratio of top flights');


