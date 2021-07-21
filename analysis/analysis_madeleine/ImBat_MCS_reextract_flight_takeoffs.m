function [takeoff_locations_cpu_re] = ImBat_MCS_reextract_flight_takeoffs(flightPaths34,outlier_flight_indexes,F,P)

% Function to re-extract the flight takeoff positions from the flight
% position vector to see if there are points later or sooner that make that
% flight includable in the takeoff.
    
% Create vector of all flight cluster IDs 
true_times = flightPaths34.AllFlightsMasterTime(flightPaths34.flight_starts_idx(:));
[tts,tts_idx]  = sort(true_times(:));

[tts_ext,tts_idx_ext] = sort(flightPaths34.AllFlightsMasterTime(:));

% Create vector of all flight days
fd = flightPaths34.day(tts_idx);
xyz_positions = flightPaths34.pos(:,:,tts_idx);
xyz_flight_starts_unsorted = squeeze(flightPaths34.flight_starts_xyz);
xyz_flight_starts = xyz_flight_starts_unsorted(:,tts_idx);
vel_unsorted = squeeze(flightPaths34.vel);
vel = vel_unsorted(:,tts_idx);
trajectories = flightPaths34.trajectoriesContinous(:,tts_idx_ext);

% Plot the outlier points in red and all start points in green
figure(); hold on;
plot3(trajectories(1,:),trajectories(2,:),trajectories(3,:),'Color',[0.7,0.7,0.7]);
values_outlier_flight_indexes = F(:,outlier_flight_indexes);
bin_outlier_flight_indexes = find(outlier_flight_indexes == 1);

for i=1:size(xyz_flight_starts,2)
    scatter3(xyz_flight_starts(1,:),xyz_flight_starts(2,:),xyz_flight_starts(3,:),[],[0 1 0]);
end
for i=1:sum(outlier_flight_indexes)
    scatter3(values_outlier_flight_indexes(1,i),values_outlier_flight_indexes(2,i),values_outlier_flight_indexes(3,i),[],[1 0 0])
end


% Plot velocity of points
vel_vec_sorted = []; vel_vec_unsorted = [];
for i=1:size(outlier_flight_indexes,2)
    if ismember(i,bin_outlier_flight_indexes)
        vel_vec_unsorted = [vel_vec_unsorted,vel_unsorted(1,i)];
        vel_vec_sorted = [vel_vec_sorted,vel(1,i)];        
    end
end
VV = zeros(size(outlier_flight_indexes,2),1);
VV(bin_outlier_flight_indexes) = vel_vec_unsorted;
disp("Starting velocities of outlier points:");
disp(vel_vec_unsorted);
disp("Outlier point indexes")
disp(bin_outlier_flight_indexes);
figure(); hold on; title("All flights' starting velocities");
temp = flightPaths34vel(1,:);%mean(flightPaths34vel(1,:));
scatter([1:size(outlier_flight_indexes,2)],temp);
scatter([1:size(outlier_flight_indexes,2)],VV');

% Examine the points that are below motion threshold
vel_unsorted_firstrow = flight_vel_unsorted(1,:);
motion_threshold = 0.4
velocity_outliers = find(vel_unsorted_firstrow > motion_threshold);

% Scatter plot the start points with the velocity outliers in red
figure();
hold on;
title("Velocity outliers")
set(gca,'XLim',[-2.9 2.9],'YLim',[-2 2.5],'ZLim',[0 2.2])
for i=1:size(xyz_flight_starts,2)
    if ismember(i,velocity_outliers)
        scatter3(xyz_flight_starts(1,i),xyz_flight_starts(2,i),xyz_flight_starts(3,i),[], [1 0 0],'filled');;
    elseif ismember(i,bin_outlier_flight_indexes)
        scatter3(xyz_flight_starts(1,i),xyz_flight_starts(2,i),xyz_flight_starts(3,i),[], [0 1 0],'filled');;        
    else
        scatter3(xyz_flight_starts(1,i),xyz_flight_starts(2,i),xyz_flight_starts(3,i),[],[0 0 1]);
    end  
end


% Find the points -200s from the takeoff
%cmap = jet(sum(outlier_flight_indexes));
%figure(); hold on; 
for i=1:sum(outlier_flight_indexes)
    OF = bin_outlier_flight_indexes(i);
    inOF = flightPaths34.flight_starts_idx(OF)-120; %

    temp = flightPaths34.AllFlightsMasterTime(inOF);
    
    scatter(flightPaths34.AllFlightsMasterTime(inOF),[ones(size(inOF))+i;zeros(size(inOF))+i],'LineWidth',3,'color',cmap(i,:));

end

unsorted_takeoff_times_shifted = flightPaths34.flight_starts_idx(bin_outlier_flight_indexes)-120;
unsorted_takeoff_times = flightPaths34.flight_starts_idx(bin_outlier_flight_indexes);

true_times_a = flightPaths34.AllFlightsMasterTime(unsorted_takeoff_times_shifted);
true_times_b = flightPaths34.AllFlightsMasterTime(unsorted_takeoff_times);

[tts_a,tts_idx_a]  = sort(true_times_a(:));
[tts_b,tts_idx_b]  = sort(true_times_b(:));







end
