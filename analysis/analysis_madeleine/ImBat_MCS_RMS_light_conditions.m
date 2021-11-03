function ImBat_MCS_RMS_light_conditions(flightPaths,dark_flights_idx,lite_flights_idx)

% 4. Plot the filtered RMS of each flight, look at this between light and
% dark.

[ss,rr] = sort(flightPaths.flight_starts_idx);
RR = flightPaths.id(rr);

flight_echos = [];
for i=2:length(flightPaths.id)
    F = flightPaths.flight_starts_idx_sorted(i);
    Fe = flightPaths.flight_ends_idx_sorted(i);
    echos = flightPaths.EcholocationVect_padded(F:Fe);
    echos(isnan(echos))=0;
    pad = NaN(1100-length(echos),1);
    flight_echos(:,i) = [echos;pad];
end
figure(); hold on; imagesc(flight_echos'); 
title("Echolocations per flight. Aligned to takeoff (green). Red lines are start and stop of dark period.");
ylabel("Flight (trial)");
xlabel("Time (ms)");
line([20 20], [1 160], 'Color', 'g','LineWidth',2);
line([1 1100], [dark_flights_idx(1) dark_flights_idx(1)], 'Color', 'r','LineWidth',2);
line([1 1100], [dark_flights_idx(end) dark_flights_idx(end)], 'Color', 'm','LineWidth',2);

    



end