function [outlier_velocity_indexes] = ImBat_MCS_find_velocity_outlier_flights(flightPaths34,outlier_flight_indexes)

flight_start_times = flightPaths34.AllFlightsMasterTime(flightPaths34.flight_starts_idx(:));
[tts,tts_idx]  = sort(flight_start_times(:));

flight_starts_xyz_unsorted = squeeze(flightPaths34.flight_starts_xyz);
flight_starts_xyz = flight_starts_xyz_unsorted(:,tts_idx);

% figure();
% scatter(flight_starts_xyz_unsorted(1,:),flight_starts_xyz_unsorted(2,:),flight_starts_xyz_unsorted(3,:));
% figure();
% scatter(flight_starts_xyz(1,:),flight_starts_xyz(2,:),flight_starts_xyz(3,:));

flight_pos_xyz = flightPaths34.pos(:,:,tts_idx);

fd = flightPaths34.day(tts_idx);

flight_vel_unsorted = squeeze(flightPaths34.vel);
flight_vel = flight_vel_unsorted(:,tts_idx);
% 
% figure();
% scatter([1:size(flight_vel_unsorted,2)],flight_vel_unsorted(1,:));
% figure();
% scatter([1:size(flight_vel,2)],flight_vel(1,:));

% IGNORE ABOVE< COPIED

% Plot first velocity datapoints of the spatial outlier points
bin_outlier_flight_indexes = find(outlier_flight_indexes==1);

vel_vec_sorted = []; vel_vec_unsorted = [];
for i=1:size(outlier_flight_indexes,2)
    if ismember(i,bin_outlier_flight_indexes)
        vel_vec_unsorted = [vel_vec_unsorted,flight_vel_unsorted(1,i)];
        vel_vec_sorted = [vel_vec_sorted,flight_vel(1,i)];        
    end
end
VV = zeros(size(outlier_flight_indexes,2),1);
VV(bin_outlier_flight_indexes) = vel_vec_unsorted;
disp("Starting velocities of outlier points:");
disp(vel_vec_unsorted);
disp("Outlier point indexes")
disp(bin_outlier_flight_indexes);
figure(); hold on; title("All flights' starting velocities. Spatial outliers in red.");
scatter([1:size(outlier_flight_indexes,2)],flight_vel_unsorted(1,:));
scatter([1:size(outlier_flight_indexes,2)],VV');

% Examine the points that are below motion threshold
vel_unsorted_firstrow = flight_vel_unsorted(1,:);
motion_threshold = 0.4;
velocity_outliers = find(vel_unsorted_firstrow > motion_threshold);

% Scatter plot the start points with the velocity outliers in red
figure();
hold on;
title("Velocity outliers (>0.4m/s) in red. Spatial outliers in green.")
set(gca,'XLim',[-2.9 2.9],'YLim',[-2 2.5],'ZLim',[0 2.2])
for i=1:size(flight_starts_xyz_unsorted,2)
    if ismember(i,velocity_outliers)
        scatter3(flight_starts_xyz(1,i),flight_starts_xyz(2,i),flight_starts_xyz(3,i),[], [1 0 0],'filled');;
    elseif ismember(i,bin_outlier_flight_indexes)
        scatter3(flight_starts_xyz(1,i),flight_starts_xyz(2,i),flight_starts_xyz(3,i),[], [0 1 0],'filled');;        
    else
        scatter3(flight_starts_xyz(1,i),flight_starts_xyz(2,i),flight_starts_xyz(3,i),[],[0 0 1]);
    end  
end

outlier_velocity_indexes = velocity_outliers;

% Calculate the overlap of velocity and position outliers.
percent_similar = 0;
for i=1:size(velocity_outliers,2)
    if ismember(velocity_outliers(i),bin_outlier_flight_indexes)
        percent_similar = percent_similar+1;
    end
end
percent_similar2 = percent_similar/size(unique([velocity_outliers,bin_outlier_flight_indexes]),2);

disp(strcat(num2str(percent_similar2*100),"% overlap between velocity and spatial outliers."));

end