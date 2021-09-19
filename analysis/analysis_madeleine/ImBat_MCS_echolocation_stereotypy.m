function ImBat_MCS_echolocation_stereotypy(flightPaths34,cluster)

% Calculate the mean amplitde of the peaks summed across all flights
% normalized by the number of flights.
[fssorted,ss] = sort(flightPaths34.flight_starts_idx)
fpid = flightPaths34.id(ss);
fpid_cluster = find(fpid==cluster);
ff=100;

clear flights_stacked
for i=1:size(fpid_cluster,1)
    start = flightPaths34.flight_starts_idx_sorted(fpid_cluster(i))-ff;
    stop = flightPaths34.flight_ends_idx_sorted(fpid_cluster(i));
    flight_dur(i) = stop-start;
    padding = zeros(800-size(flightPaths34.EcholocationVect_padded(start:stop),1),1);
    flights_stacked(:,i) = [flightPaths34.EcholocationVect_padded(start:stop);padding];
end

for i=1:size(flights_stacked,2)
    voc = flights_stacked(:,i);
    voc(isnan(voc))=0;
    flights_stacked(:,i) = voc;
end
    
FSS = sum(flights_stacked,2);

figure(); hold on; subplot(2,1,1);imagesc(FSS'); colorbar;
title(strcat("Average echolocation stereotypy of cluster ", num2str(cluster)));
xline(ff,'r','LineWidth',2);
subplot(2,1,2);
hist(flight_dur);

figure();
plot(FSS);