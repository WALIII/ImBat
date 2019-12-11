function [flightTransitions] = ImBat_transition_flights(flightPaths)

%this is to merge specific clusters if 2 or more of them are the same but
%were separated by the flight k-means
flightPaths.clusterIndex{2}=cat(2,flightPaths.clusterIndex{2},flightPaths.clusterIndex{3});
flightPaths.clusterIndex{3} =[];
flightPaths.clusterIndex= flightPaths.clusterIndex(~cellfun('isempty',flightPaths.clusterIndex));

%initialize the cluster order
clustOrder = zeros(length(flightPaths.flight_starts_idx(1,:)),1);
for n = 1:length(flightPaths.flight_starts_idx(1,:))
    for clust_i = 1:length(flightPaths.clusterIndex)
            clustTemp = find(flightPaths.clusterIndex{clust_i}==n); 
            if ~isempty(clustTemp)
                clustOrder(n) = clust_i;
            end
    end  
end

AtoB = [];
AtonotB = [];
BtoA = [];
BtonotA = [];
notAtonotB = [];
for n = 1:length(flightPaths.flight_starts_idx(1,:))-1
    if clustOrder(n) == 1 && clustOrder(n+1) == 2 || clustOrder(n) == 1 && clustOrder(n+1) == 3
        AtoB = [AtoB n];
    elseif clustOrder(n) == 1 && clustOrder(n+1) ~= 2 || clustOrder(n) == 1 && clustOrder(n+1) ~= 3
        AtonotB = [AtonotB n];     
    elseif clustOrder(n) == 2 && clustOrder(n+1) == 1 || clustOrder(n) == 3 && clustOrder(n+1) == 1
        BtoA = [BtoA n];
    elseif clustOrder(n) == 2 && clustOrder(n+1) ~= 1 || clustOrder(n) == 3 && clustOrder(n+1) ~= 1
        BtonotA = [BtonotA n];
    else 
        notAtonotB = [notAtonotB n];
    end
end

for traj_i = 1:length(flightPaths.clusterIndex)
    
end

flightTransitions.clustOrder = clustOrder;
flightTransitions.AtoB = AtoB; 
flightTransitions.AtonotB = AtonotB;
flightTransitions.BtoA = BtoA;
flightTransitions.BtonotA = BtonotA;
flightTransitions.notAtonotB = notAtonotB;