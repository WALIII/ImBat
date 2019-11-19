function [clustOrder] = transition_flights

clustOrder = zeros(length(flightPaths.flight_starts_idx(1,:),1));
for n = 1:length(flightPaths.flight_starts_idx(1,:))
    for clust_i = 1:length(flightPaths.clusterIndex)
        if isempty(clustOrder(n))
            clustOrder(n) = find(flightPaths.clusterIndex{clust_i}==n); 
        end
    end  
end