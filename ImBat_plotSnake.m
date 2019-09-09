function ImBat_plotSnake

for i = 1:length(results.C(:,1))
    for ii = 1:length(flightPaths.flight_starts_idx)
        trace(ii) = results.C(i,flightPaths.flight_starts_idx(ii):flightPaths.flight_ends_idx(ii))
   end
    %meanTrace =  
end