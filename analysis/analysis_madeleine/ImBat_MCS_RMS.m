function ImBat_MCS_RMS(flightPaths)

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
    pad = NaN(size(flightPaths.vel,2)-length(echos),1);
    flight_echos(:,i) = [echos;pad];
end
figure(); hold on; imagesc(flight_echos'); 
title("Echolocations per flight.");
ylabel("Flight #");
xlabel("Time");

end