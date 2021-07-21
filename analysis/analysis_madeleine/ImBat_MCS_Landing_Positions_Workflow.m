% Script to identify spatially weird and velocity-weird outliers in a
% subset of the data, cluster, and plot them.

flight_start_times = flightPaths34.AllFlightsMasterTime(flightPaths34.flight_starts_idx(:));
[tts,tts_idx]  = sort(flight_start_times(:));

flight_starts_xyz_unsorted = squeeze(flightPaths34.flight_starts_xyz);
flight_starts_xyz = flight_starts_xyz_unsorted(:,tts_idx);

figure();
scatter(flight_starts_xyz_unsorted(1,:),flight_starts_xyz_unsorted(2,:),flight_starts_xyz_unsorted(3,:));
figure();
scatter(flight_starts_xyz(1,:),flight_starts_xyz(2,:),flight_starts_xyz(3,:));

flight_pos_xyz = flightPaths34.pos(:,:,tts_idx);

fd = flightPaths34.day(tts_idx);

flight_vel_unsorted = squeeze(flightPaths34.vel);
flight_vel = flight_vel_unsorted(:,tts_idx);

figure();
scatter([1:size(flight_vel_unsorted,2)],flight_vel_unsorted(1,:));
figure();
scatter([1:size(flight_vel,2)],flight_vel(1,:));

% Identify the spatially-weird outliers (NOT SORTED BY TIME)
[outlier_flight_indexes]=ImBat_MCS_find_outlier_flights(flightPaths34);
% For now, exclude these points

% Identify the velocity-weird outliers
[outlier_velocity_indexes] = ImBat_MCS_find_velocity_outlier_flights(flightPaths34,outlier_flight_indexes);

% Cluster the non-outliers
bin_outlier_flight_indexes = find(outlier_flight_indexes==1)
non_outlier_indexes= [1:size(tts_idx,1)];
non_outlier_indexes(bin_outlier_flight_indexes) = [];

% Remove the indexes of the outlier flights (from the SORTED flight start
% xys position vector)
flight_starts_xyz_pruned = flight_starts_xyz;
flight_starts_xyz_pruned(:,bin_outlier_flight_indexes) = [];

clear Xvar Yvar Zvar
max_clusters=100;
for j=1:max_clusters
    C = kmeans(flight_starts_xyz_pruned',j);
    for i=1:max(C)
        Xvar(j,i) = var(F(1,find(C==i)));
        Yvar(j,i) = var(F(2,find(C==i)));
        Zvar(j,i) = var(F(3,find(C==i)));
    end
end
Xvar(Xvar==0) = NaN;
Yvar(Yvar==0) = NaN;
Zvar(Zvar==0) = NaN;
figure(); hold on;
title("Median variance of X, Y, and Z coordinates within a cluster for 1:100 clusters");
plot([1:max_clusters],nanmedian(Xvar,2));
plot([1:max_clusters],nanmedian(Yvar,2));
plot([1:max_clusters],nanmedian(Zvar,2));

% It seems like 8 clusters is nice 
k=8;
%C = kmeans(Fxy',k);
KC = kmeans(flight_starts_xyz_pruned',k);
figure();
hold on;
colormap = jet(k);
for i=1:max(KC)
    scatter3(flight_starts_xyz_pruned(1,find(KC==i)),flight_starts_xyz_pruned(2,find(KC==i)),flight_starts_xyz_pruned(3,find(KC==i)),[], colormap(i,:),'filled','MarkerEdgeColor','k');
end




