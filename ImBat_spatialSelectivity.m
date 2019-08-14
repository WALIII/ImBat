function ImBat_spatialSelectivity(cellData,alignment,flightPaths,flightFeedersStartStop)

binSizeX = 200; %size of bin for rate code along the x axis (mm)
offset = 0.1; % account for slow calcium estimation ~move locations back 100ms in time... This is the knob to turn for 'prospective' coding...
topROI = 30;
vidFS = 120; %sampling frequency of video to calculate firing rate at each bin


for cn = 1:length(flightPaths.clusterIndex) %for each cluster
    for n = 1:length(cellData.S(:,1)) %for each neuron from cellData
        [~,xy] = find(cellData.S(n,:)>8);  % get time neuron is active
        Spike_times = alignment.out.video_times(xy)-offset; % convert this to 'spike time'
        for i = 1:size(Spike_times,1)
            try
                % Find the closest 'Location time' in x axis to the 'Spike time'
                [minValue(:,i),closestIndex(:,i)] = min(abs(alignment.out.Location_time-Spike_times(i)));
                LX(i) = alignment.out.flights(closestIndex(:,i),1);
            catch % we need this if Spiketime occurs before/after the location tracking was on..
                continue
            end
        end
        for ci = 1:length(flightPaths.clusterIndex{cn}) %for each flight in that cluster
            %reverse the negatives if flying backwards across room
            if flightPaths.flight_ends_xyz(flightPaths.clusterIndex{cn}(ci)) < flightPaths.flight_starts_xyz(flightPaths.clusterIndex{cn}(ci))
                flightPaths.flight_ends_xyz(flightPaths.clusterIndex{cn}(ci)) = -flightPaths.flight_ends_xyz(flightPaths.clusterIndex{cn}(ci));
                flightPaths.flight_starts_xyz(flightPaths.clusterIndex{cn}(ci)) = -flightPaths.flight_starts_xyz(flightPaths.clusterIndex{cn}(ci));
            end
            %find the start & end positions of flights in the video data
            [minValueStart(:),closestIndexStart(:)] = min(abs(alignment.out.Location_time(flightPaths.flight_starts_idx(flightPaths.clusterIndex{cn}(ci)))-alignment.out.video_times));
            [minValueStart(:),closestIndexEnd(:)] = min(abs(alignment.out.Location_time(flightPaths.flight_ends_idx(flightPaths.clusterIndex{cn}(ci)))-alignment.out.video_times));
            
            %bins of the x position, spiking position, and bin limits
            xBin = alignment.out.flights(flightPaths.flight_starts_idx(flightPaths.clusterIndex{cn}(ci))):binSizeX:alignment.out.flights(flightPaths.flight_ends_idx(flightPaths.clusterIndex{cn}(ci)));
            xPos = alignment.out.flights(flightPaths.flight_starts_idx(flightPaths.clusterIndex{cn}(ci)):flightPaths.flight_ends_idx(flightPaths.clusterIndex{cn}(ci)));
            xSpike = LX; %spike positions in x-axis
            yBinSpike = histc(xSpike,xBin);
            yBinPos = histc(xPos,xBin);
            firingRateBin = 1000*yBinSpike/(yBinPos/vidFS); %calculate firing rate of
        
            %for every bin along the x axis
            
            
                for  flightPaths.flight_ends_xyz(flightPaths.clusterIndex{cn}(ci)):binSizeX:flightPaths.flight_starts_xyz(flightPaths.clusterIndex{cn}(ci))

                  
                    
                    
                    
                end
            elseif flightPaths.flight_ends_xyz(flightPaths.clusterIndex{cn}(ci)) > flightPaths.flight_starts_xyz(flightPaths.clusterIndex{cn}(ci))
                for  flightPaths.flight_starts_xyz(flightPaths.clusterIndex{cn}(ci)):binSizeX:flightPaths.flight_ends_xyz(flightPaths.clusterIndex{cn}(ci))
                    
                    
                end
                
            end
            
            
        end
        
    end
end




