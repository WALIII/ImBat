fileList = dir('*.gpx');
colors = jet(length(fileList));
webmap('openstreetmap')
for track_i = 1:length(fileList)
tracks(track_i) = gpxread(fileList(track_i).name, 'Index', 1:2);
wmline(tracks(track_i), 'Color', colors(track_i,:))
%fprintf(fileList(track_i).name)
pause
end
%wmline(tracks, 'Color', colors(track_i,:))

%webmap('openstreetmap')
%colors = {'cyan'};
%wmline(tracks, 'Color', colors)