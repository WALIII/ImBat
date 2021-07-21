function [outlier_flight_indexes] = ImBat_MCS_find_outlier_flights(flightPaths34)

% Function that identifies the outlier flights, in order to change the
% start time or exclude that flight from the clustering.

flight_start_times = flightPaths34.AllFlightsMasterTime(flightPaths34.flight_starts_idx(:));
[tts,tts_idx]  = sort(flight_start_times(:));

flight_starts_xyz_unsorted = squeeze(flightPaths34.flight_starts_xyz);
flight_starts_xyz = flight_starts_xyz_unsorted(:,tts_idx);

xv = [-2.4,2.05,2.05,-2.4,-2.4];
yv = [2.2,2.2,-1.6,-1.6,2.2];

xq = [flight_starts_xyz(1,:)];
yq = [flight_starts_xyz(2,:)];

[in,on] = inpolygon(xq,yq,xv,yv);

figure(); hold on
title("Spatial outliers are red crosses. Outline is where bat is suspected to actually be on the wall");
plot(xv,yv) % polygon
plot(xq(in&~on),yq(in&~on),'r+') % points strictly inside
plot(xq(on),yq(on),'k*') % points on edge
plot(xq(~in),yq(~in),'bo') % points outside
hold off

outlier_flight_indexes=in;

disp(strcat(num2str(sum(in))," flights identified as not starting on the wall or feeder"));

end