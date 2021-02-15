function [Clust2save, AllFLights2save] = ImBat_Session_FlightRate(flightPaths)

% Plot the rate of flights over the session for figure 1

% WAL3
% 2/12/2021

% Get flightpaths of type:
FlightPaths = 1:4;



% plot unclustered
% for day2use = 1:max(flightPaths.day);
for i = 1:size(FlightPaths,2)+1
    cluster = i;
    %         idx =  find (flightPaths.id == cluster & flightPaths.day == day2use);
    idx =  find (flightPaths.id == cluster);
    
    FF2 = flightPaths.flight_starts_idx(idx);
    
    hold on;
    histout = histogram(flightPaths.AllFlightsTime(FF2)/60,'NumBins',60,'BinLimits',[1,60]);
    % Clust2save(cluster,day2use,:) = histout.Values;
    Clust2save(cluster,:) = histout.Values;
    
end

FF2 = flightPaths.flight_starts_idx(:);
histout_all = histogram(flightPaths.AllFlightsTime(FF2)/60,'NumBins',60,'BinLimits',[1,60]);
AllFLights2save(:) = histout_all.Values;


figure();
 plot(Clust2save(3,:)./AllFLights2save);
