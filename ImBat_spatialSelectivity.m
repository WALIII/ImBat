function ImBat_spatialSelectivity(cellData,trackData,flightPaths,flightFeedersStartStop)

binSizeX = 300; %size of bin for rate code along the x axis (mm)

for cn = 1:length(flightPaths.clusterIndex) %for each cluster
    for n = 1:length(cellData.S(:,1)) %for each neuron from cellData
        for ci = 1:length(flightPaths.clusterIndex{cn}) %for each flight in that cluster
            %for every bin along the x axis
            if flightPaths.flight_ends_xyz(flightPaths.clusterIndex{cn}(ci)) < flightPaths.flight_starts_xyz(flightPaths.clusterIndex{cn}(ci))
                for  flightPaths.flight_ends_xyz(flightPaths.clusterIndex{cn}(ci)):binSizeX:flightPaths.flight_starts_xyz(flightPaths.clusterIndex{cn}(ci))
                    
                    
                end
            elseif flightPaths.flight_ends_xyz(flightPaths.clusterIndex{cn}(ci)) > flightPaths.flight_starts_xyz(flightPaths.clusterIndex{cn}(ci))
                for  flightPaths.flight_starts_xyz(flightPaths.clusterIndex{cn}(ci)):binSizeX:flightPaths.flight_ends_xyz(flightPaths.clusterIndex{cn}(ci))
                    
                    
                end
                
            end
            
            
        end
        
    end
end




