function [Location] = ImBat_formatTracking(Markers);
% ImBat_format_tracking.m

% Take in marker data from 'cortex' software, and conver it to a 3D position in space.

% WAL3
% D05102019


range = size( Markers,1);

for i = 1:range;

if Markers(i,1,1) == 0 && Markers(i,1,2) == 0 && Markers(i,1,3) == 0
    Location(i,:) = NaN;
else
    for ii = 1:3
        Location(i,ii) = squeeze(Markers(i,1,ii));
    end
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
comet3(Location(:,1),Location(:,2),Location(:,3));
