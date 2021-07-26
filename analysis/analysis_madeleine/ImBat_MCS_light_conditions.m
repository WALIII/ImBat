function ImBat_MCS_light_conditions(flightPaths)

% Function to ask 
% 1. Are there more echolocations on the flights in the dark in general?
% 2. ^ Normalized by flight length
% 3. How do the echolocations look on the stereotyped flights that persist across
% light/dark?

lights_off_range_C = [151200,273600]; lights_off_range_L = [144000,288000];
lights_off_range = lights_off_range_L;
dark_flights_idx = []; lite_flights_idx = []; q_flights_idx = [];
for i=1:length(flightPaths.flight_starts_idx_sorted)
    F = flightPaths.flight_starts_idx_sorted(i);
    Fe = flightPaths.flight_ends_idx_sorted(i);
    if all(F >= lights_off_range(1) & F <= lights_off_range(2)) & all(Fe >= lights_off_range(1) & Fe <= lights_off_range(2))
        dark_flights_idx = [dark_flights_idx,i];
    else
        lite_flights_idx = [lite_flights_idx,i];
    end
end

% Get number of echolocations for each flight in the light
for i=1:length(lite_flights_idx)
    f_idx = flightPaths.flight_starts_idx_sorted(lite_flights_idx(i));
    fe_idx = flightPaths.flight_ends_idx_sorted(lite_flights_idx(i));
    num_echos_light(i) = length(find(flightPaths.EcholocationVect_padded(f_idx:fe_idx)==0.5));
    num_echos_light_norm(i) = length(find(flightPaths.EcholocationVect_padded(f_idx:fe_idx)==0.5))/length(f_idx:fe_idx);
end

% Get number of echolocations for each flight in the "dark"
for i=1:length(dark_flights_idx)
    f_idx = flightPaths.flight_starts_idx_sorted(dark_flights_idx(i));
    fe_idx = flightPaths.flight_ends_idx_sorted(dark_flights_idx(i));
    num_echos_dark(i) = length(find(flightPaths.EcholocationVect_padded(f_idx:fe_idx)==0.5));
    num_echos_dark_norm(i) = length(find(flightPaths.EcholocationVect_padded(f_idx:fe_idx)==0.5))/length(f_idx:fe_idx);
end

% Plot the distributions
figure(); hold on; title("# Echolocations/flight in light vs low-light");
xlabel("# echolocations per flight");
ylabel("Frequency");
legend('Light','Dark');
histogram(num_echos_light,'DisplayName','Light'); histogram(num_echos_dark,'DisplayName','Dark');

% Plot the distributions (normalized by flight duration)
figure(); hold on; title("# Echolocations/flight in light vs low-light normalized");
xlabel("# echolocations per flight normalized by flight duration");
ylabel("Frequency");
legend('Light','Dark');
histogram(num_echos_light_norm,'DisplayName','Light'); histogram(num_echos_dark_norm,'DisplayName','Dark');

% For a given stereotyped flight cluster, plot the flights in the dark versus in the flight

cluster_to_plot = 2;
plot_counter = 0;
figure(); hold on; title(strcat("Echolocations on Cluster ",num2str(cluster_to_plot)," flights"));
flight_subset = find(flightPaths.id==cluster_to_plot);
for i=1:size(flightPaths.id,1)
    if flightPaths.id(i)==cluster_to_plot
        plot_counter = plot_counter+1;
        timeline_idx = find(flightPaths.flight_starts_idx_sorted==flightPaths.flight_starts_idx(i));
    	subplot(round(sqrt(size(flight_subset,1)))+1,round(sqrt(size(flight_subset,1))),plot_counter);
        hold on; 
        if ismember(timeline_idx,dark_flights_idx)
            title(strcat("DARK ",num2str(timeline_idx)));
        else
            title(num2str(timeline_idx));
        end
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


    
        



end