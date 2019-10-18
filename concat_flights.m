function [rateMap1D] = concat_flights(flightPaths, alignment,rateMap1D)
cn = 1; %cluster number to concatenate
n = 1; %neuron number

%make cells for each cluster
allFlightsConcat = cell(1,length(flightPaths.clusterIndex));
allSpikePosConcat = cell(1,length(flightPaths.clusterIndex));
%for each cluster of flights
for cn = 1:2%length(flightPaths.clusterIndex)
    %concatenate all flights
    for f = 1:length(flightPaths.clusterIndex{cn})
        allFlightsConcat{cn} = vertcat(allFlightsConcat{cn},alignment.flights(flightPaths.flight_starts_idx(flightPaths.clusterIndex{cn}(f)):flightPaths.flight_ends_idx(flightPaths.clusterIndex{cn}(f)),:));
    end
    %allFlightsConcat{cn} = cell(1,length(rateMap1D.spikePos{cn})); %make cells for each neuron
    
    
    allSpikePosConcat{cn} = cell(1,length(rateMap1D.spikePos{cn})); %make cells for each neuron
    for n = 1:10%length(rateMap1D.spikePos{cn})
        %concatenate all spike positions
        for i = 1:length(flightPaths.clusterIndex{cn})
            allSpikePosConcat{cn}{n} = vertcat(allSpikePosConcat{cn}{n},rateMap1D.spikePos{cn}{n}{i});
        end
    end
    
end

rateMap1D.allFlightsConcat = allFlightsConcat;
rateMap1D.allSpikePosConcat = allSpikePosConcat;