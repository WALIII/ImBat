function [Location, Location_interp] = ImBat_formatTracking(Markers);
% ImBat_format_tracking.m

% Take in marker data from 'cortex' software, and conver it to a 3D position in space.

% WAL3
% D05102019


range = size( Markers,1);
reform_tracking =0;

for i = 1:range;
    
    if Markers(i,1,1) == 0 && Markers(i,1,2) == 0 && Markers(i,1,3) == 0
        Location(i,:) = NaN;
    else
        for ii = 1:3
            Location(i,ii) = squeeze(Markers(i,1,ii));
        end
    end
end

if max(max(max(Markers))) ==0
    disp(' no valid tracking data...');
    reform_tracking =1;
else
    st = min(find(Location(:,1)>1));
    for ii = 1:3;
        Location_interp(:,ii) = fillmissing(Location(:,ii)','linear');
    end
    % remove early vals
    try
        Location_interp(1:st,:) = ones(size(Location_interp(1:st,:))).*Location_interp(st+1,:);
    catch
        Location_interp(1:st,:) = ones(size(Location_interp(1:st,:))).*Location_interp(st,:);
    end
end


% Exclude the NaNs
for i = 1:range
    if i ==1;
        Location(i,1:3) = 0;
    else
        if isnan(Location(i,1));
            Location(i,:) = Location(i-1,:);
        end
    end
end


% Plot the results:
% comet3(Location(:,1),Location(:,2),Location(:,3));
if reform_tracking ==1;
    Location_interp = Location;
end

% Smooth Location data
 for i = 1:3
     Location3(:,i) = movmean(Location(:,i),30);
 end

