function [rateMap1D] = ImBat_spatialSelectivity(cellData,alignment,flightPaths)

xBinSize = 200; %size of bin for rate code along the x axis (mm)
offset = 0.1; % account for slow calcium estimation ~move locations back 100ms in time... This is the knob to turn for 'prospective' coding...
topROI = 30;
vidFS = 120; %sampling frequency of video to calculate firing rate at each bin
activityThresh = 5;
%occupancyThresh = 0.2; %min number of seconds that bat must occupy a space for a 3d rate map
w = gausswin(23,1.5);

for cn = 1:3%length(flightPaths.clusterIndex) %for each cluster
    %reverse the negatives if flying backwards across room
    for i = 1:length(flightPaths.clusterIndex{cn})
    if flightPaths.flight_ends_xyz(flightPaths.clusterIndex{cn}(i)) < flightPaths.flight_starts_xyz(flightPaths.clusterIndex{cn}(i))
        flightPaths.flight_ends_xyz(flightPaths.clusterIndex{cn}(i)) = -flightPaths.flight_ends_xyz(flightPaths.clusterIndex{cn}(i));
        flightPaths.flight_starts_xyz(flightPaths.clusterIndex{cn}(i)) = -flightPaths.flight_starts_xyz(flightPaths.clusterIndex{cn}(i));
    end
    end
    %min and max x position for bin vector
    xBinMin = min(flightPaths.flight_starts_xyz(flightPaths.clusterIndex{cn}));
    xBinMax = max(flightPaths.flight_ends_xyz(flightPaths.clusterIndex{cn}));
    xBin{cn} = xBinMin:xBinSize:xBinMax;        
    
    %for each neuron from cellData, extract spiking position and calculate
    %1d rate map
    for n = 1:40%length(cellData.S(:,1)) 
        [~,xy] = find(cellData.S(n,:)>activityThresh);  % get time neuron is active
        Spike_times = alignment.video_times(xy)-offset; % convert this to 'spike time' in video clock
        for i = 1:size(Spike_times,1)
            try
                % Find the closest 'Location time' in x axis to the 'Spike time'
                [minValue(:,i),closestIndex(:,i)] = min(abs(alignment.Location_time-Spike_times(i)));
                LX(i) = alignment.flights(closestIndex(:,i),1);
            catch % we need this if Spiketime occurs before/after the location tracking was on..
                continue
            end
        end
        
        yBinSpike{cn}{n} = zeros(length(flightPaths.clusterIndex{cn}),length(xBin{cn})-1); %set up array size for the firing rates
        yBinPos{cn}{n} = zeros(length(flightPaths.clusterIndex{cn}),length(xBin{cn})-1); %set up array size for the firing rates
        yBinPosSec{cn}{n} = zeros(length(flightPaths.clusterIndex{cn}),length(xBin{cn})-1); %set up array size for the firing rates
        firingRate{cn}{n} = zeros(length(flightPaths.clusterIndex{cn}),length(xBin{cn})-1); %set up array size for the firing rates
        %for each flight in that cluster, perform 1D firing rate calculation
        for ci = 1:length(flightPaths.clusterIndex{cn})
            %bins of the x position, spiking position, and bin limits
            xPos{cn}{n}{ci} = alignment.flights(flightPaths.flight_starts_idx(flightPaths.clusterIndex{cn}(ci)):flightPaths.flight_ends_idx(flightPaths.clusterIndex{cn}(ci)));
            xSpike{cn}{n}{ci} = LX; %spike positions in x-axis
            yBinSpike{cn}{n}(ci,:) = histcounts(xSpike{cn}{n}{ci},xBin{cn}); %sum all spikes in each bin
            yBinPos{cn}{n}(ci,:) = histcounts(xPos{cn}{n}{ci},xBin{cn}); %sum all occupancy times in each bin
            yBinPosSec{cn}{n}(ci,:) = yBinPos{cn}{n}(ci,:)/vidFS;
            aa = yBinPosSec{cn}{n}(ci,:);
            aa(yBinPosSec{cn}{n}(ci,:) == 0)=nan; %turn 0 into nan so you can divide for rate
            yBinPosSec{cn}{n}(ci,:) = aa;
            firingRate{cn}{n}(ci,:) = yBinSpike{cn}{n}(ci,:)./yBinPosSec{cn}{n}(ci,:); %calculate firing rate of each bin           
        
            
        end
        %sum up all spike and position occupancies to calculate firing rate of overall sum
        %of all flights in each cluster per neuron
        yBinSpikeSum{cn}(n,:) = sum(firingRate{cn}{n},1);
        yBinPosSum{cn}(n,:) = sum(yBinPos{cn}{n},1);
        yBinPosSecSum{cn}(n,:) = sum(yBinPosSec{cn}{n},1);
        firingRateSum{cn}(n,:) = yBinSpikeSum{cn}(n,:)./yBinPosSecSum{cn}(n,:);
        %switch nan out for 0 for the smoothing
        bb = firingRateSum{cn}(n,:); 
        bb(isnan(firingRateSum{cn}(n,:))) = 0;
        firingRateSum{cn}(n,:) = bb;
        %smooth the summed firing rate with guassian filter
        firingRateSmooth{cn}(n,:) = filter(w,1,firingRateSum{cn}(n,:));
        bb = firingRateSmooth{cn}(n,:);
            bb(firingRateSmooth{cn}(n,:) == 0)=nan; %turn 0 into nan so you can divide for rate
            firingRateSmooth{cn}(n,:) = bb;
        
    end
        
end



rateMap1D.yBinSpike = yBinSpike;
rateMap1D.yBinPos = yBinPos;
rateMap1D.yBinPosSec = yBinPosSec;
rateMap1D.firingRate = firingRate;
rateMap1D.xBin = xBin; 
rateMap1D.yBinSpikeSum = yBinSpikeSum;
rateMap1D.yBinPosSum = yBinPosSum;
rateMap1D.yBinPosSecSum = yBinPosSecSum; 
rateMap1D.firingRateSum = firingRateSum;
rateMap1D.firingRateSmooth = firingRateSmooth;
rateMap1D.xPos = xPos;
rateMap1D.xSpike = xSpike;
