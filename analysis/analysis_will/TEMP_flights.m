
function TEMP_flights(flightPaths,FlightPaths2plot,Days2use);




figure();
col = distinguishable_colors(length(unique(flightPaths.id)));
alpha2use = 0.4;

% for all flights






hold on;
A = flightPaths.tracjectoriesRaw*1000;
counter = 1;

% plot all flights first:

for day2use =  Days2use;

    % plot all flights
     %subplot(length(Day2use),1,day2use); 

     subplot(1,length(Days2use),counter); 
     counter = counter+1;
        
for clust2use = FlightPaths2plot;
    % Hilight stable flights
hold on;
Ind2use =  find(flightPaths.day == day2use );
%Ind2use2 = find(flightPaths.day == day2use & flightPaths.id == clust2use);
axis off
for iii = 1:length(Ind2use)  
    bound = flightPaths.flight_starts_idx(Ind2use(iii)):flightPaths.flight_ends_idx(Ind2use(iii));
    plot1 =  plot3(A(1,bound),A(2,bound),A(3,bound),'color',[0 0 0 alpha2use]); % plot all flights
end
% for ii = 1: length(Ind2use2);
%     bound2 = flightPaths.flight_starts_idx(Ind2use2(ii)):flightPaths.flight_ends_idx(Ind2use2(ii));
%     plot2 =  plot3(A(1,bound2),A(2,bound2),A(3,bound2),'color',col{clust2use-1},'LineWidth',2); % plot all flights
% end
   % view( -37.5000,30)

end
   xlim([-2800 2800]);
    ylim([-2800 2800]);
     zlim([0 2500]);
end

axis off

figure(); 
counter = 1;

for day2use =  Days2use;

    % plot all flights
     %subplot(length(Day2use),1,day2use); 

     subplot(1,length(Days2use),counter); 
     counter = counter+1;
        
for clust2use = FlightPaths2plot;
    % Hilight stable flights
hold on;
title(strcat("Clusterable flights in color. Day(s): "," ",num2str(Days2use)));
Ind2use =  find(flightPaths.day == day2use );
Ind2use2 = find(flightPaths.day == day2use & flightPaths.id == clust2use);
axis off
for iii = 1:length(Ind2use)  
    bound = flightPaths.flight_starts_idx(Ind2use(iii)):flightPaths.flight_ends_idx(Ind2use(iii));
    plot1 =  plot3(A(1,bound),A(2,bound),A(3,bound),'color',[0 0 0 alpha2use]); % plot all flights
end
for ii = 1: length(Ind2use2);
    bound2 = flightPaths.flight_starts_idx(Ind2use2(ii)):flightPaths.flight_ends_idx(Ind2use2(ii));
    plot2 =  plot3(A(1,bound2),A(2,bound2),A(3,bound2),'color',col(clust2use,:),'LineWidth',1.5); % plot all flights
end
   % view( -37.5000,30)

end
   xlim([-2800 2800]);
    ylim([-2800 2800]);
     zlim([0 2500]);
end

axis off

