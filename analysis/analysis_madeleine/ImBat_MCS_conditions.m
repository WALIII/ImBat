function ImBat_MCS_conditions(flightPaths,cluster_to_plot)

% Function to ask 
% 1. Are there more echolocations on the flights in the dark in general?
% 2. ^ Normalized by flight length
% 3. How do the echolocations look on the stereotyped flights that persist across
% light/dark?

% For a given stereotyped flight cluster, plot the flights in the dark versus in the flight

plot_counter = 0;
figure(); hold on; title(strcat("Echolocations on Cluster ",num2str(cluster_to_plot)," flights"));
flight_subset = find(flightPaths.id==cluster_to_plot);
for i=1:size(flightPaths.id,1)
    if flightPaths.id(i)==cluster_to_plot
        plot_counter = plot_counter+1;
        timeline_idx = find(flightPaths.flight_starts_idx_sorted==flightPaths.flight_starts_idx(i));
    	subplot(round(sqrt(size(flight_subset,1)))+1,round(sqrt(size(flight_subset,1))),plot_counter);
        hold on; 
        title(num2str(timeline_idx));
        plot3(flightPaths.pos(1,:,i),flightPaths.pos(2,:,i),flightPaths.pos(3,:,i));
        plot3(flightPaths.trajectoriesSpline(1,flightPaths.flight_starts_idx(i):flightPaths.flight_ends_idx(i)),flightPaths.trajectoriesSpline(2,flightPaths.flight_starts_idx(i):flightPaths.flight_ends_idx(i)),flightPaths.trajectoriesSpline(3,flightPaths.flight_starts_idx(i):flightPaths.flight_ends_idx(i)),':');
        fpe = flightPaths.EcholocationVect_padded(flightPaths.flight_starts_idx(i):flightPaths.flight_ends_idx(i));
        fpe(fpe==0)=NaN;
        fpe_idxs = [flightPaths.flight_starts_idx(i):flightPaths.flight_ends_idx(i)];
        for j=1:size(fpe,1)
            if fpe(j)==0.5000
                plot3(flightPaths.trajectoriesSpline(1,fpe_idxs(j)),flightPaths.trajectoriesSpline(2,fpe_idxs(j)),flightPaths.trajectoriesSpline(3,fpe_idxs(j)),'*r');
            end
        end  
    end
end