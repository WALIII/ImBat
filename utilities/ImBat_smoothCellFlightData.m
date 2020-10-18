function ImBat_smoothCellFlightData(flightPathsAll)

    %smooth and zscore the velocity
    smoothVelocity = zscore(smooth(flightPathsAll.batSpeed,100));

    %image_preprocessing, smooth and average #cells spike responses
    data(i).image_data.numCells = numCells;
    data(i).image_data.avgFiring = zscore(smooth(mean((full(data(i).image_data.S(1:numCells,:))),1),50));%-min(data(p).image_data.S(1:numCells,:));
